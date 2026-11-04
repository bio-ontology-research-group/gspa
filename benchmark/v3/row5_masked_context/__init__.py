"""Row-5 masked-locus reaction-context benchmark package."""
import hashlib

if not hasattr(hashlib, 'file_digest'):
    def _file_digest(fileobj, digest):
        h = hashlib.new(digest) if isinstance(digest, str) else digest()
        while True:
            chunk = fileobj.read(65536)
            if not chunk:
                break
            h.update(chunk)
        return h

    hashlib.file_digest = _file_digest  # type: ignore[attr-defined]
