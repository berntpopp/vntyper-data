"""Derive the adVNTR MUC1 model from a chromosome 1 FASTA and pinned coordinates.

The model used to ship as an opaque binary seed (`seeds/vntr_db_advntr.zip`) that
nobody could regenerate, and the hg38 copy was wrong: it carried hg19's `repeats` and
`pattern` verbatim, describing 840 bp of an array that GRCh38 carries at 3,525 bp.
See berntpopp/vntyper-data#1 and hassansaei/VNtyper#268.

Two things this script establishes:

1. **The model is derivable.** Everything in the row comes from chr1 plus the pinned
   constants below -- flanks, segments, pattern, coordinates.
2. **The derivation is trustworthy**, because it reproduces the shipped hg19 database
   content byte-for-byte. hg19 is correct for its own array, so it is a free regression
   test of the unit-boundary rule. See `--verify`.

The unit boundary is the 13 bp anchor `GGC[CT]GAGGTGACA`. Cutting the hg19 array at
every occurrence reproduces its 14 stored segments exactly, including the irregular
78, 48 and 54 bp ones. Note the honest limitation: hg19 round-trips at mismatch
tolerance 0, 1, 2 and 3 alike, because its array contains no degenerate anchors. The
round-trip validates the bounds, flanks, `ref_start` and serialisation -- it is NOT
evidence for the mismatch tolerance. Only the benchmark can settle that.

`advntr addmodel` is not used: it runs `train_classifier_threshold`, documented as
taking hours on a human genome, and the shipped model has `scaled_score = 0.0`, i.e.
that training was never applied. Writing the row directly reproduces current behaviour
without the multi-hour step. `scaled_score` stays 0.0 deliberately.
"""
from __future__ import print_function

import argparse
import hashlib
import os
import re
import sqlite3
import struct
import sys


# --- Pinned inputs -----------------------------------------------------------------
#
# The array bounds were recovered from the flanks stored in the shipped databases
# (left_flanking ends at ref_start; right_flanking begins one base after the array) and
# are verified against chr1 by --verify-provenance. Two independent lines of support:
# on BOTH assemblies the array start coincides exactly with an anchor occurrence, and
# the first and last 24 bp of the array are identical across the two assemblies.
# config.json's `vntr_array_coords` are deliberately NOT used: they are self-described
# as heuristic, drive coverage normalisation (#222), and fall 16/22 bp INSIDE the first
# unit on hg19/hg38 respectively.

ANCHOR = re.compile(r'GGC[CT]GAGGTGACA')
ANCHOR_CANONICAL = ('GGCCGAGGTGACA', 'GGCTGAGGTGACA')

VID = 25561
GENE = 'MUC1'
ANNOTATION = 'Coding'
FLANK = 500
SCALED_SCORE = 0.0
CHROM = 'chr1'

# 1-based inclusive array bounds.
ARRAYS = {
    'hg19': (155160984, 155161823),   # 840 bp
    'hg38': (155188508, 155192032),   # 3525 bp
}

# The two shipped databases differ in SQL whitespace, page size and the
# `nonoverlapping` literal. Reproduced per assembly so the hg19 round-trip is exact.
# NOTE: hg38 stores '1', which advntr/models.py evaluates as non_overlapping=False
# (it tests `== 'True'`). That is a real latent defect, deliberately NOT corrected
# here -- this change alters one thing at a time. Filed separately.
LEGACY = {
    'hg19': {
        'page_size': 1024,
        'nonoverlapping': 'True',
        'schema': ('CREATE TABLE vntrs(id INTEGER PRIMARY KEY, nonoverlapping TEXT,\n'
                   '                       chromosome TEXT, ref_start INTEGER, gene_name TEXT, annotation TEXT, '
                   'pattern TEXT, left_flanking TEXT, right_flanking TEXT, repeats TEXT, scaled_score REAL DEFAULT 0)'),
    },
    'hg38': {
        'page_size': 4096,
        'nonoverlapping': '1',
        'schema': ('CREATE TABLE vntrs(id INTEGER PRIMARY KEY, nonoverlapping TEXT, chromosome TEXT, '
                   'ref_start INTEGER, gene_name TEXT,\n'
                   '    annotation TEXT, pattern TEXT, left_flanking TEXT, right_flanking TEXT, repeats TEXT, '
                   'scaled_score REAL default 0)'),
    },
}

# The canonical 60 bp repeat unit. Identical in both shipped databases.
PATTERN = 'GGCCGAGGTGACACCGTGGGCTGGGGGGGCGGTGGAGCCCGGGGCCGGCCTGGTGTCCGG'

# The 9 header bytes a clean rebuild cannot reproduce, and why. Every one is a property
# of the writing SQLite library, not of the model: all content pages are byte-identical.
HEADER_DELTAS = (
    (24, 28, 'file change counter'),
    (43, 44, 'schema cookie'),
    (92, 96, 'version-valid-for'),
    (96, 100, 'SQLite version stamp'),
)


def read_chr1(path, start, end):
    """Return the 1-based inclusive [start, end] slice of the single record in path.

    Accepts either a whole-chromosome FASTA or a pre-cut region whose header carries
    `chrom:start-end`, so the derivation can run against either.
    """
    header = None
    chunks = []
    with open(path) as handle:
        for line in handle:
            if line.startswith('>'):
                if header is not None:
                    break
                header = line[1:].strip()
            else:
                chunks.append(line.strip())
    sequence = ''.join(chunks).upper()
    offset = 1
    match = re.match(r'\S*?:(\d+)-(\d+)', header or '')
    if match:
        offset = int(match.group(1))
    lo = start - offset
    hi = end - offset + 1
    if lo < 0 or hi > len(sequence):
        raise ValueError('%s does not cover %d-%d' % (path, start, end))
    return sequence[lo:hi]


def segment(array, max_mismatch=0, min_gap=20):
    """Cut the array at every anchor occurrence; return the segments in order.

    max_mismatch relaxes the anchor to a Hamming distance over the 13-mer, which
    recovers unit boundaries whose anchor has mutated. min_gap stops a near-miss
    adjacent to a real anchor from cutting twice.
    """
    if max_mismatch == 0:
        cuts = [m.start() for m in ANCHOR.finditer(array)]
    else:
        cuts = []
        for i in range(len(array) - 12):
            window = array[i:i + 13]
            best = min(sum(1 for a, b in zip(window, canon) if a != b)
                       for canon in ANCHOR_CANONICAL)
            if best <= max_mismatch:
                cuts.append(i)
    kept = []
    for cut in cuts:
        if not kept or cut - kept[-1] >= min_gap:
            kept.append(cut)
    if not kept or kept[0] != 0:
        kept = [0] + kept
    return [array[a:b] for a, b in zip(kept, kept[1:] + [len(array)])]


def build_row(chr1_path, assembly, mode, max_mismatch, chr1_hg19=None):
    """Assemble the single database row for an assembly."""
    start, end = ARRAYS[assembly]
    array = read_chr1(chr1_path, start, end)
    left = read_chr1(chr1_path, start - FLANK, start - 1)
    right = read_chr1(chr1_path, end + 1, end + FLANK)

    if mode == 'legacy':
        # Reproduce what the assembly's own array yields. For hg19 this is the shipped
        # content; for hg38 the shipped content is hg19's and is NOT reproducible from
        # GRCh38 sequence -- which is the defect.
        segments = segment(array, max_mismatch=0)
    elif mode == 'inherited':
        # The shipped hg38 model carries hg19's segments verbatim -- that is the copy
        # bug. Deriving them from the hg19 array rather than copying a blob keeps the
        # artefact reproducible while leaving the HMM byte-identical to what has been
        # benchmarked. Only the coordinates and ref_end come from this assembly.
        #
        # This is deliberately NOT a claim that these units describe GRCh38. They do
        # not. It records the historical content in derivable form so a model can be
        # regenerated, while the fetch window is corrected by ref_end.
        if not chr1_hg19:
            raise ValueError('--segmentation inherited needs --chr1-hg19: the units come '
                             'from the GRCh37 array, and reading GRCh37 coordinates out '
                             'of a GRCh38 FASTA would silently yield different sequence')
        hg19_start, hg19_end = ARRAYS['hg19']
        hg19_array = read_chr1(chr1_hg19, hg19_start, hg19_end)
        segments = segment(hg19_array, max_mismatch=0)
        return {
            'id': VID,
            'nonoverlapping': LEGACY[assembly]['nonoverlapping'],
            'chromosome': CHROM,
            'ref_start': start - 1,
            'gene_name': GENE,
            'annotation': ANNOTATION,
            'pattern': PATTERN,
            'left_flanking': left,
            'right_flanking': right,
            'repeats': ','.join(segments),
            'scaled_score': SCALED_SCORE,
            'ref_end': end,
            'n_segments': len(segments),
            'n_distinct': len(set(segments)),
            'max_segment': max(len(s) for s in segments),
        }
    else:
        segments = segment(array, max_mismatch=max_mismatch)

    assert ''.join(segments) == array, 'segments must tile the array exactly'
    return {
        'id': VID,
        'nonoverlapping': LEGACY[assembly]['nonoverlapping'],
        'chromosome': CHROM,
        'ref_start': start - 1,          # 0-based
        'gene_name': GENE,
        'annotation': ANNOTATION,
        'pattern': PATTERN,
        'left_flanking': left,
        'right_flanking': right,
        'repeats': ','.join(segments),
        'scaled_score': SCALED_SCORE,
        'ref_end': end,                  # 0-based exclusive == 1-based inclusive end
        'n_segments': len(segments),
        'n_distinct': len(set(segments)),
        'max_segment': max(len(s) for s in segments),
    }


def write_db(path, row, assembly, schema_version):
    """Write the single-row database.

    schema_version 'legacy' reproduces the shipped table exactly. 'v2' writes a
    `vntrs_v2` table carrying ref_end -- a separate table so that an adVNTR predating
    ref_end fails loudly instead of silently ignoring the column and recreating the
    truncated window.
    """
    if os.path.exists(path):
        os.remove(path)
    db = sqlite3.connect(path)
    db.execute('PRAGMA page_size=%d' % LEGACY[assembly]['page_size'])
    columns = ['id', 'nonoverlapping', 'chromosome', 'ref_start', 'gene_name',
               'annotation', 'pattern', 'left_flanking', 'right_flanking', 'repeats',
               'scaled_score']
    if schema_version == 'legacy':
        db.execute(LEGACY[assembly]['schema'])
        table = 'vntrs'
    else:
        table = 'vntrs_v2'
        db.execute('CREATE TABLE vntrs_v2(id INTEGER PRIMARY KEY, nonoverlapping TEXT, '
                   'chromosome TEXT, ref_start INTEGER, gene_name TEXT, annotation TEXT, '
                   'pattern TEXT, left_flanking TEXT, right_flanking TEXT, repeats TEXT, '
                   'scaled_score REAL DEFAULT 0, ref_end INTEGER)')
        columns.append('ref_end')
    db.execute('INSERT INTO %s (%s) VALUES (%s)'
               % (table, ','.join(columns), ','.join(['?'] * len(columns))),
               [row[c] for c in columns])
    db.commit()
    db.close()


def compare(derived_path, shipped_path):
    """Compare a derived database with a shipped one; return (ok, report lines)."""
    a = open(shipped_path, 'rb').read()
    b = open(derived_path, 'rb').read()
    lines = []
    if len(a) != len(b):
        return False, ['size differs: shipped %d, derived %d' % (len(a), len(b))]

    diffs = [i for i in range(len(a)) if a[i:i + 1] != b[i:i + 1]]
    allowed = set()
    for lo, hi, _name in HEADER_DELTAS:
        allowed.update(range(lo, hi))
    unexpected = [i for i in diffs if i not in allowed]

    lines.append('differing bytes: %d' % len(diffs))
    for lo, hi, name in HEADER_DELTAS:
        sa = struct.unpack('>I', a[lo:hi].ljust(4, b'\0'))[0]
        sb = struct.unpack('>I', b[lo:hi].ljust(4, b'\0'))[0]
        if sa != sb:
            lines.append('  offs %3d-%-3d %-22s shipped=%-10d derived=%d'
                         % (lo, hi - 1, name, sa, sb))
    lines.append('unexpected differing bytes (outside the SQLite header): %d'
                 % len(unexpected))
    if unexpected:
        lines.append('  offsets: %s' % unexpected[:20])

    # Content equality, independent of file bytes.
    def dump(path):
        con = sqlite3.connect(path)
        schema = list(con.execute("SELECT sql FROM sqlite_master WHERE type='table'"))[0][0]
        cols = [r[1] for r in con.execute('PRAGMA table_info(vntrs)')]
        rows = list(con.execute('SELECT %s FROM vntrs' % ','.join(cols)))
        return schema, cols, rows

    same_content = dump(derived_path) == dump(shipped_path)
    lines.append('logical content identical (schema SQL + every column value): %s'
                 % same_content)
    lines.append('content digest shipped=%s derived=%s'
                 % (hashlib.sha256(repr(dump(shipped_path)).encode()).hexdigest()[:16],
                    hashlib.sha256(repr(dump(derived_path)).encode()).hexdigest()[:16]))
    return (not unexpected) and same_content, lines


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--chr1', required=True,
                        help='chromosome 1 FASTA (whole chromosome or a region covering the locus)')
    parser.add_argument('--assembly', required=True, choices=sorted(ARRAYS))
    parser.add_argument('--out', help='database to write')
    parser.add_argument('--schema', default='v2', choices=('legacy', 'v2'))
    parser.add_argument('--segmentation', default='anchor',
                        choices=('legacy', 'anchor', 'inherited'),
                        help="'legacy' pins the exact-anchor rule on this assembly's array; "
                             "'anchor' honours --max-mismatch; 'inherited' reproduces the "
                             "historical hg19-derived units (needs --chr1-hg19)")
    parser.add_argument('--chr1-hg19', help='GRCh37 chr1 FASTA, required by --segmentation inherited')
    parser.add_argument('--max-mismatch', type=int, default=0,
                        help='anchor Hamming tolerance (0 = exact)')
    parser.add_argument('--verify', metavar='SHIPPED',
                        help='compare the derived database against a shipped one and report')
    args = parser.parse_args(argv)

    row = build_row(args.chr1, args.assembly, args.segmentation, args.max_mismatch,
                    chr1_hg19=getattr(args, 'chr1_hg19', None))
    print('%s: array %d-%d (%d bp), %d segments, %d distinct, longest %d bp'
          % (args.assembly, ARRAYS[args.assembly][0], ARRAYS[args.assembly][1],
             ARRAYS[args.assembly][1] - ARRAYS[args.assembly][0] + 1,
             row['n_segments'], row['n_distinct'], row['max_segment']))

    out = args.out
    if args.verify and not out:
        out = args.verify + '.derived'
    if not out:
        parser.error('--out is required unless --verify is given')
    write_db(out, row, args.assembly, args.schema)
    print('wrote %s (schema=%s)' % (out, args.schema))

    if args.verify:
        ok, lines = compare(out, args.verify)
        for line in lines:
            print('  ' + line)
        print('ROUND TRIP: %s' % ('PASS' if ok else 'FAIL'))
        return 0 if ok else 1
    return 0


if __name__ == '__main__':
    sys.exit(main())
