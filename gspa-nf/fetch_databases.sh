#!/usr/bin/env bash
#
# fetch_databases.sh — pull GSPA neural-predictor model artefacts based on
# the versioned manifest at gspa-nf/database_manifest.tsv.
#
# Usage:
#   fetch_databases.sh --predictor proteinfer [--predictor truth] \
#                      --dest ~/gspa-db [--manifest URL_OR_PATH] [--no-extract]
#
# Behaviour per matching manifest row:
#   1. mkdir -p $DEST/$predictor/$version
#   2. curl -L -C - $url -o file (resumable)
#   3. sha256sum -c (skipped if sha256 == REPLACE_WITH_SHA256, with a warning)
#   4. tar xzf if .tar.gz
#   5. gunzip if .gz (and not .tar.gz)
#   6. print the matching `databases.config` export lines
#
# Requires: bash, curl, awk, sha256sum, tar.

set -euo pipefail

PREDICTORS=()
DEST=""
MANIFEST="$(dirname "$0")/database_manifest.tsv"
EXTRACT=1

# Default URL prefix used by manifest rows (override with env GSPA_DB_BASE).
: "${GSPA_DB_BASE:=https://gspa.bio2vec.net/db}"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --predictor)  PREDICTORS+=("$2"); shift 2 ;;
        --dest)       DEST="$2"; shift 2 ;;
        --manifest)   MANIFEST="$2"; shift 2 ;;
        --no-extract) EXTRACT=0; shift ;;
        -h|--help)
            sed -n '1,/^set/p' "$0" | sed 's/^# \{0,1\}//'
            exit 0
            ;;
        *) echo "unknown arg: $1" >&2; exit 2 ;;
    esac
done

[[ -n "$DEST" ]] || { echo "--dest is required" >&2; exit 2; }
[[ ${#PREDICTORS[@]} -gt 0 ]] || { echo "at least one --predictor is required" >&2; exit 2; }

mkdir -p "$DEST"

# Load manifest into memory; expand ${GSPA_DB_BASE}.
if [[ "$MANIFEST" =~ ^https?:// ]]; then
    MANIFEST_LOCAL="$DEST/.manifest.tsv"
    curl -fL "$MANIFEST" -o "$MANIFEST_LOCAL"
    MANIFEST="$MANIFEST_LOCAL"
fi
[[ -s "$MANIFEST" ]] || { echo "manifest empty: $MANIFEST" >&2; exit 2; }

echo "manifest: $MANIFEST"
echo "dest:     $DEST"
echo "GSPA_DB_BASE: $GSPA_DB_BASE"
echo

EXPORT_LINES=()

# Skip header; for each row whose predictor matches one of the requested,
# fetch + verify + extract.
while IFS=$'\t' read -r predictor artefact version url sha256 size_bytes required notes; do
    [[ "$predictor" == "predictor" ]] && continue
    [[ -z "$predictor" ]] && continue

    match=0
    for p in "${PREDICTORS[@]}"; do
        if [[ "$p" == "$predictor" ]]; then match=1; break; fi
    done
    [[ "$match" -eq 0 ]] && continue

    # Expand ${GSPA_DB_BASE}
    url_expanded="${url//\$\{GSPA_DB_BASE\}/$GSPA_DB_BASE}"
    out_dir="$DEST/$predictor/$version"
    out_file="$out_dir/$(basename "$url_expanded")"
    mkdir -p "$out_dir"

    echo "==> $predictor / $artefact ($version)"
    echo "    url:  $url_expanded"
    echo "    dest: $out_file"

    # Resumable download
    curl -fL -C - "$url_expanded" -o "$out_file"

    # sha256 check (skipped if placeholder)
    if [[ "$sha256" == "REPLACE_WITH_SHA256" || -z "$sha256" ]]; then
        echo "    WARN: manifest has no sha256 — skipping checksum verify" >&2
    else
        echo "$sha256  $out_file" | sha256sum -c -
    fi

    # Extract / decompress
    if [[ "$EXTRACT" -eq 1 ]]; then
        case "$out_file" in
            *.tar.gz|*.tgz)
                tar xzf "$out_file" -C "$out_dir"
                ;;
            *.gz)
                gunzip -kf "$out_file"
                ;;
        esac
    fi

    # Build databases.config export hint
    final="$out_file"
    case "$out_file" in
        *.tar.gz|*.tgz) final="$out_dir" ;;
        *.gz)           final="${out_file%.gz}" ;;
    esac
    case "$predictor/$artefact" in
        esm2-deepgoplus/head_ckpt) EXPORT_LINES+=("    esm2_dgp_ckpt        = '$final'") ;;
        esm2-deepgoplus/go_terms)  EXPORT_LINES+=("    esm2_dgp_terms       = '$final'") ;;
        esm2-centroid/centroid_db) EXPORT_LINES+=("    esm2_centroid_db     = '$final'") ;;
        proteinfer/saved_model)    EXPORT_LINES+=("    proteinfer_model_dir = '$final'") ;;
        clean/weights)             EXPORT_LINES+=("    clean_model_dir      = '$final'") ;;
        truth/swissprot_go_ec)     EXPORT_LINES+=("    swissprot_go_ec      = '$final'") ;;
        truth/smoke_fixture)       EXPORT_LINES+=("    truth_smoke          = '$final'") ;;
        go_aspect_map/go_obo)      EXPORT_LINES+=("    go_obo               = '$final'") ;;
    esac
done < "$MANIFEST"

echo
echo "==== Add the following to databases.config (under params { … }) ===="
for line in "${EXPORT_LINES[@]}"; do
    echo "$line"
done
echo "===================================================================="
