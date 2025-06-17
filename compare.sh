#!/bin/bash
# compare.sh: compare two versions of prodigal
# ./compare.sh <old_binary> <input_file>

if [[ $# -lt 2 ]]; then
    echo "Usage: $0 <old_binary> <input_file> [additional opts]"
    echo "Example: $0 ./prodigal.windows.exe anthus_aco.fas -t test/anthus_aco.train"
    exit 1
fi

OLD_BINARY="$1"
INPUT_FILE="$2"
OUT_PREFIX="test/$(basename "$INPUT_FILE")"
shift 2

if ! command -v "$OLD_BINARY" >/dev/null; then
    echo "Error: Old binary '$OLD_BINARY' is not executable."
    exit 1
fi
if [[ ! -f $INPUT_FILE ]]; then
    echo "Error: Input file '$INPUT_FILE' does not exist."
    exit 1
fi

mkdir -p "test"

if ! "$OLD_BINARY" -i "$INPUT_FILE" -f gbk -o "$OUT_PREFIX".old.gbk -d "$OUT_PREFIX".old.ffn -a "$OUT_PREFIX".old.faa "$@"; then
    echo "Error: Old binary failed to run."
    exit 1
fi

if ! ./prodigal -i "$INPUT_FILE" -f gbk -o "$OUT_PREFIX".new.gbk -d "$OUT_PREFIX".new.ffn -a "$OUT_PREFIX".new.faa "$@"; then
    echo "Error: New binary failed to run."
    exit 1
fi

# normalize line endings
sed -i -e 's/\r//g' "$OUT_PREFIX".{old,new}.*

if ! diff -u "$OUT_PREFIX".old.gbk "$OUT_PREFIX".new.gbk > "$OUT_PREFIX".gbk.diff; then
    echo "GBK files differ."
else
    echo "GBK files are identical."
fi

if ! diff -u "$OUT_PREFIX".old.ffn "$OUT_PREFIX".new.ffn > "$OUT_PREFIX".ffn.diff; then
    echo "FFN files differ."
else
    echo "FFN files are identical."
fi

if ! diff -u "$OUT_PREFIX".old.faa "$OUT_PREFIX".new.faa > "$OUT_PREFIX".faa.diff; then
    echo "FAA files differ."
else
    echo "FAA files are identical."
fi
