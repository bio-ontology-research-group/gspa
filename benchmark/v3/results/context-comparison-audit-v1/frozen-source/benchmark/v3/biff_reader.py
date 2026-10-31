"""Minimal read-only reader for legacy Excel (.xls, OLE2 + BIFF8) workbooks, standard library only.

Written because the primary-source supplement of Wintermute and Silver 2010 is a binary .xls and the pinned
benchmark environment has no xlrd.  Reads cell values (strings and numbers) of every sheet; formats, formulas
(beyond their cached value), charts and dates are not interpreted.  Anything unexpected raises ``BiffError``
instead of returning a partial table.
"""

from __future__ import annotations

import struct
from typing import Any

OLE_MAGIC = b"\xd0\xcf\x11\xe0\xa1\xb1\x1a\xe1"
ENDOFCHAIN, FREESECT = 0xFFFFFFFE, 0xFFFFFFFF


class BiffError(ValueError):
    pass


def _ole_stream(data: bytes, wanted: tuple[str, ...] = ("Workbook", "Book")) -> bytes:
    if data[:8] != OLE_MAGIC:
        raise BiffError("not an OLE2 compound file")
    ssz = 1 << struct.unpack_from("<H", data, 30)[0]
    msz = 1 << struct.unpack_from("<H", data, 32)[0]
    dir_start, cutoff, mini_start = (struct.unpack_from("<I", data, o)[0] for o in (48, 56, 60))
    difat_start, n_difat = (struct.unpack_from("<I", data, o)[0] for o in (68, 72))

    def sector(i: int) -> bytes:
        off = 512 + i * ssz
        if off + ssz > len(data):
            raise BiffError("sector outside file")
        return data[off:off + ssz]

    difat = list(struct.unpack_from("<109I", data, 76))
    nxt = difat_start
    for _ in range(n_difat):
        if nxt in (ENDOFCHAIN, FREESECT):
            break
        block = struct.unpack(f"<{ssz // 4}I", sector(nxt))
        difat += list(block[:-1])
        nxt = block[-1]
    fat: list[int] = []
    for s in difat:
        if s in (FREESECT, ENDOFCHAIN):
            continue
        fat += struct.unpack(f"<{ssz // 4}I", sector(s))

    def chain(start: int, table: list[int]) -> list[int]:
        out, seen = [], set()
        while start not in (ENDOFCHAIN, FREESECT):
            if start in seen or start >= len(table):
                raise BiffError("corrupt sector chain")
            seen.add(start)
            out.append(start)
            start = table[start]
        return out

    directory = b"".join(sector(s) for s in chain(dir_start, fat))
    entries = []
    for i in range(len(directory) // 128):
        e = directory[i * 128:(i + 1) * 128]
        nlen = struct.unpack_from("<H", e, 64)[0]
        name = e[:max(nlen - 2, 0)].decode("utf-16le", "replace")
        entries.append((name, e[66], struct.unpack_from("<I", e, 116)[0], struct.unpack_from("<I", e, 120)[0]))
    root = next((e for e in entries if e[1] == 5), None)
    target = next((e for e in entries if e[1] == 2 and e[0] in wanted), None)
    if root is None or target is None:
        raise BiffError("no Workbook stream in compound file")
    _, _, start, size = target
    if size >= cutoff:
        return b"".join(sector(s) for s in chain(start, fat))[:size]
    mini_fat: list[int] = []
    for s in chain(mini_start, fat):
        mini_fat += struct.unpack(f"<{ssz // 4}I", sector(s))
    mini_stream = b"".join(sector(s) for s in chain(root[2], fat))
    buf = b"".join(mini_stream[s * msz:(s + 1) * msz] for s in chain(start, mini_fat))
    return buf[:size]


def _records(stream: bytes, pos: int = 0):
    while pos + 4 <= len(stream):
        rtype, ln = struct.unpack_from("<HH", stream, pos)
        yield rtype, stream[pos + 4:pos + 4 + ln], pos
        pos += 4 + ln


def _rk(v: int) -> float:
    if v & 2:
        x = float(struct.unpack("<i", struct.pack("<I", v & 0xFFFFFFFC))[0] >> 2)
    else:
        x = struct.unpack("<d", struct.pack("<Q", (v & 0xFFFFFFFC) << 32))[0]
    return x / 100 if v & 1 else x


class _Chunks:
    """SST payload with CONTINUE boundaries, so strings that straddle a boundary decode correctly."""

    def __init__(self, parts: list[bytes]):
        self.parts, self.i, self.o = parts, 0, 0

    def _roll(self) -> bool:
        while self.i < len(self.parts) and self.o >= len(self.parts[self.i]):
            self.i, self.o = self.i + 1, 0
            return True
        return False

    def take(self, n: int) -> bytes:
        out = b""
        while n > 0:
            if self.i >= len(self.parts):
                raise BiffError("SST ended inside a string")
            avail = len(self.parts[self.i]) - self.o
            if avail == 0:
                self.i, self.o = self.i + 1, 0
                continue
            k = min(n, avail)
            out += self.parts[self.i][self.o:self.o + k]
            self.o += k
            n -= k
        return out

    def string(self) -> str:
        cch, flags = struct.unpack("<HB", self.take(3))
        runs = struct.unpack("<H", self.take(2))[0] if flags & 8 else 0
        ext = struct.unpack("<I", self.take(4))[0] if flags & 4 else 0
        wide, chars = bool(flags & 1), []
        while cch > 0:
            if self.o >= len(self.parts[self.i]):          # continuation chunk begins with a new encoding flag
                self.i, self.o = self.i + 1, 0
                if self.i >= len(self.parts):
                    raise BiffError("SST ended inside a string")
                wide = bool(self.parts[self.i][0] & 1)
                self.o = 1
                continue
            width = 2 if wide else 1
            k = min(cch, (len(self.parts[self.i]) - self.o) // width)
            raw = self.take(k * width)
            chars.append(raw.decode("utf-16le" if wide else "latin-1"))
            cch -= k
        self.take(4 * runs + ext)
        return "".join(chars)


def read_xls(data: bytes) -> dict[str, list[list[Any]]]:
    """sheet name -> rows (lists, None for empty cells)."""
    wb = _ole_stream(data)
    sheets: list[tuple[str, int]] = []
    sst: list[str] = []
    recs = list(_records(wb))
    i = 0
    while i < len(recs):
        rtype, body, _ = recs[i]
        if rtype == 0x0085:                                   # BOUNDSHEET
            off = struct.unpack_from("<I", body, 0)[0]
            cch, flags = body[6], body[7]
            nm = body[8:8 + cch * (2 if flags & 1 else 1)].decode("utf-16le" if flags & 1 else "latin-1")
            sheets.append((nm, off))
        elif rtype == 0x00FC:                                 # SST (+ CONTINUE)
            parts = [body[8:]]
            total_unique = struct.unpack_from("<II", body, 0)[1]
            while i + 1 < len(recs) and recs[i + 1][0] == 0x003C:
                i += 1
                parts.append(recs[i][1])
            ch = _Chunks(parts)
            sst = [ch.string() for _ in range(total_unique)]
        i += 1
    if not sheets:
        raise BiffError("no sheets found")
    out: dict[str, list[list[Any]]] = {}
    for name, off in sheets:
        cells: dict[tuple[int, int], Any] = {}
        pending, depth = None, 0
        for rtype, body, _pos in _records(wb, off):
            if rtype == 0x0809:                               # BOF; nested BOF/EOF pairs are charts or dialogs
                depth += 1
                continue
            if rtype == 0x000A:                               # EOF
                depth -= 1
                if depth <= 0:
                    break
                continue
            if depth != 1:
                continue
            if rtype == 0x00FD:
                r, c, _x, k = struct.unpack_from("<HHHI", body)
                cells[r, c] = sst[k]
            elif rtype == 0x0203:
                r, c, _x, v = struct.unpack_from("<HHHd", body)
                cells[r, c] = v
            elif rtype == 0x027E:
                r, c, _x, v = struct.unpack_from("<HHHI", body)
                cells[r, c] = _rk(v)
            elif rtype == 0x00BD:
                r, c0 = struct.unpack_from("<HH", body)
                n = (len(body) - 6) // 6
                for j in range(n):
                    cells[r, c0 + j] = _rk(struct.unpack_from("<I", body, 4 + 6 * j + 2)[0])
            elif rtype == 0x0204:
                r, c, _x, cch = struct.unpack_from("<HHHH", body)
                flags = body[8]
                cells[r, c] = body[9:9 + cch * (2 if flags & 1 else 1)].decode("utf-16le" if flags & 1 else "latin-1")
            elif rtype == 0x0006:
                r, c = struct.unpack_from("<HH", body)
                if body[12:14] == b"\xff\xff" and body[6] == 0:
                    pending = (r, c)
                elif body[12:14] != b"\xff\xff":
                    cells[r, c] = struct.unpack_from("<d", body, 6)[0]
            elif rtype == 0x0207 and pending:
                cch, flags = struct.unpack_from("<HB", body)
                cells[pending] = body[3:3 + cch * (2 if flags & 1 else 1)].decode("utf-16le" if flags & 1 else "latin-1")
                pending = None
        if not cells:
            out[name] = []
            continue
        nr, nc = max(r for r, _ in cells) + 1, max(c for _, c in cells) + 1
        out[name] = [[cells.get((r, c)) for c in range(nc)] for r in range(nr)]
    return out
