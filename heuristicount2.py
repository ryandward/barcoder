#!/usr/bin/env python3
"""
heuristicount2 - flank-anchored, mismatch-tolerant barcode/sgRNA counter for pooled screens.

Successor to heuristicount.py. What changed and why (see heuristicount.py for the original):

  * FLANK-ANCHORED localization (per read) instead of a single fixed offset.
    The original sliced every read at one auto-detected coordinate, so staggered /
    phased-primer libraries (deliberately variable spacer length, common in Illumina
    amplicon prep to boost cluster diversity) were undercounted by ~1/N. v2 finds the
    constant flanking primer sequence in each read and reads the barcode beside it,
    so the barcode can sit at any offset.

  * MISMATCH TOLERANCE in the barcode body (default 1, configurable with -m).
    The original required a byte-for-byte match of flank+barcode+flank, silently
    dropping every read that carried a sequencing error in that ~50 bp window
    (~half of reads at a realistic single-substitution rate). v2 recovers them via
    bounded Hamming search; a read within range of two or more library barcodes is
    flagged AMBIGUOUS and skipped, never miscounted.

  * NO "undocumented" bin. The original emitted flank-matching-but-unknown sequences
    as discovered barcodes; with mismatch recovery in place these are just sequencing
    noise, so they are summarized (as the unmatched fraction), not enumerated.

  * COMPLETE, headered, deterministic output: one row per library barcode, sorted,
    INCLUDING zero-count guides (a dropped-out guide is the primary CRISPRi signal),
    ready to paste into a MAGeCK/edgeR count matrix. Optional guide-ID column.

  * Native .tsv/.csv library loading via --column (the shipped Example_Libraries TSVs
    work directly), FASTA, or one-sequence-per-line; loud failure on malformed input.

  * Streaming multiprocessing (lazy chunk pull, static data shipped once per worker),
    a summary table whose value column actually renders under non-TTY width, and
    pipeline-safe exit codes. Uses a minimal rich logger to avoid the babel/locale
    crash the shared Logger triggers under a POSIX 'C' locale.
"""

import argparse
import csv
import gzip
import os
import platform
import sys
from collections import Counter, defaultdict
from contextlib import nullcontext
from datetime import datetime
from itertools import combinations, islice, product
from multiprocessing import Pool

import zstandard as zstd
import rich.table
from rich.console import Console
from rich.table import Table

__version__ = "2.0.0"

ISSUES_URL = "https://github.com/ryandward/barcoder/issues"
_console = Console(stderr=True)


# --------------------------------------------------------------------------- #
# Minimal logging (rich, no babel -> no locale crash)
# --------------------------------------------------------------------------- #
def info(msg):
    _console.print(f"[bold blue]INFO[/]    {msg}")


def warn(msg):
    _console.print(f"[bold yellow]WARN[/]    {msg}")


def error(msg):
    _console.print(f"[bold red]ERROR[/]   {msg}")


# --------------------------------------------------------------------------- #
# Sequence helpers
# --------------------------------------------------------------------------- #
_COMPLEMENT = str.maketrans("ACGTN", "TGCAN")


def rev_comp(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


def hamming_neighbors(seq, max_mm):
    """Yield every sequence within Hamming distance 1..max_mm of seq (ACGT substitutions)."""
    if max_mm <= 0:
        return
    n = len(seq)
    bases = "ACGT"
    for d in range(1, max_mm + 1):
        for combo in combinations(range(n), d):
            choices = [[b for b in bases if b != seq[p]] for p in combo]
            for repl in product(*choices):
                s = list(seq)
                for p, c in zip(combo, repl):
                    s[p] = c
                yield "".join(s)


def low_quality_positions(barcode, qual, min_q):
    """Indices of barcode bases that are untrustworthy: an 'N', or Phred < min_q.

    These are the ONLY positions where a deviation from a library barcode may be
    forgiven. High-confidence bases must match exactly, so guides that differ only by
    a designed (intentional) mismatch stay distinct as long as that base is called
    confidently.
    """
    if qual is not None:
        return [i for i, ch in enumerate(barcode)
                if ch == "N" or (ord(qual[i]) - 33) < min_q]
    return [i for i, ch in enumerate(barcode) if ch == "N"]


def assign(barcode, qual, lib_set, min_q, max_low):
    """Quality-gated mapping of an observed barcode to a library barcode.

    Exact match wins outright. Otherwise the only deviations considered are at
    low-confidence positions (see low_quality_positions): those bases are treated as
    wildcards and every fill is looked up. Returns (library_barcode, 'exact'|'recovered'),
    (None, 'ambiguous') if >1 library barcode is reachable that way, or
    (None, 'unmatched') if none is.
    """
    if barcode in lib_set:
        return barcode, "exact"
    if min_q <= 0:
        return None, "unmatched"  # quality gating off -> exact only
    low = low_quality_positions(barcode, qual, min_q)
    if not low or len(low) > max_low:
        return None, "unmatched"
    chars = list(barcode)
    hit = None
    for repl in product("ACGT", repeat=len(low)):
        for p, c in zip(low, repl):
            chars[p] = c
        cand = "".join(chars)
        if cand in lib_set:
            if hit is None:
                hit = cand
            elif cand != hit:
                return None, "ambiguous"
    if hit is not None:
        return hit, "recovered"
    return None, "unmatched"


# --------------------------------------------------------------------------- #
# File I/O
# --------------------------------------------------------------------------- #
def _open(path, mode="rt"):
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    if path.endswith(".zst"):
        return zstd.open(path, mode)
    return open(path, mode)


def _strip_compression(path):
    return path[:-3] if path.endswith(".gz") else (path[:-4] if path.endswith(".zst") else path)


def iter_seqs(path):
    """Yield (sequence, quality) from a FASTQ or .reads file (optionally compressed).

    quality is the raw Phred+33 string for FASTQ, or None for .reads (no quality).
    """
    base = _strip_compression(path)
    if base.endswith((".fastq", ".fq")):
        ftype = "fastq"
    elif base.endswith(".reads"):
        ftype = "reads"
    else:
        raise ValueError(
            f'"{path}" is not a supported reads file (.fastq, .fq, .reads, optionally .gz/.zst).'
        )
    name = os.path.basename(path)
    with _open(path, "rt") as f:
        if ftype == "fastq":
            while True:
                header = f.readline()
                if not header:
                    break
                seq = f.readline()
                plus = f.readline()
                qual = f.readline()
                if not seq or not plus or not qual:
                    warn(f"Truncated final FASTQ record in {name}; ignoring it.")
                    break
                if header[0] != "@" or plus[0] != "+":
                    raise ValueError(
                        f"Malformed FASTQ in {name}: expected '@'/'+' record lines but got "
                        f"{header[:20]!r}/{plus[:20]!r}. The file is misaligned or not FASTQ."
                    )
                yield seq.strip().upper(), qual.strip()
        else:
            for line in f:
                s = line.strip()
                if s:
                    yield s.upper(), None


def iter_chunks(path1, path2, chunk_size):
    """Yield (reads1, reads2) chunks; reads2 is None for single-end.

    For paired-end input the two files must have equal length. If they diverge we
    stop at the common prefix and warn loudly rather than silently dropping the
    trailing unpaired reads (a mismatched-files error).
    """
    g1 = iter_seqs(path1)
    g2 = iter_seqs(path2) if path2 else None
    while True:
        r1 = list(islice(g1, chunk_size))
        if g2 is None:
            if not r1:
                break
            yield r1, None
            continue
        r2 = list(islice(g2, chunk_size))
        if not r1 and not r2:
            break
        if len(r1) != len(r2):
            n = min(len(r1), len(r2))
            warn(
                f"Paired read files have unequal length: stopping after {n} aligned "
                f"pair(s) in this chunk and ignoring {abs(len(r1) - len(r2))} unpaired "
                f"trailing read(s). Check that the R1/R2 files match."
            )
            if n:
                yield r1[:n], r2[:n]
            break
        yield r1, r2


# --------------------------------------------------------------------------- #
# Library loading
# --------------------------------------------------------------------------- #
def load_library(path, column, id_column):
    """Load barcodes from .tsv/.csv (header + --column), FASTA, or one-per-line.

    Returns (spacer_to_id: dict[str, str|None], duplicate_count: int).
    """
    base = _strip_compression(path)
    ext = os.path.splitext(base)[1].lower()
    entries = []  # (id_or_None, seq)

    with _open(path, "rt") as f:
        if ext in (".tsv", ".csv"):
            delim = "\t" if ext == ".tsv" else ","
            reader = csv.DictReader(f, delimiter=delim)
            fields = reader.fieldnames or []
            if column not in fields:
                raise ValueError(
                    f"Column '{column}' not found in {os.path.basename(path)}. "
                    f"Available columns: {fields}. Use --column to pick the barcode column."
                )
            if id_column and id_column not in fields:
                raise ValueError(
                    f"--id-column '{id_column}' not found in {os.path.basename(path)}. "
                    f"Available columns: {fields}."
                )
            for row in reader:
                seq = (row.get(column) or "").strip().upper()
                if not seq:
                    continue
                rid = (row.get(id_column) or "").strip() if id_column else None
                entries.append((rid or None, seq))
        elif ext in (".fasta", ".fa", ".fna"):
            rid = None
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    rid = line[1:].split()[0] if len(line) > 1 else None
                else:
                    entries.append((rid, line.upper()))
        else:
            # plain: one sequence per line; tolerate stray FASTA headers
            for line in f:
                line = line.strip()
                if not line or line.startswith(">"):
                    continue
                entries.append((None, line.upper()))

    valid = set("ACGTN")
    spacer_to_id = {}
    duplicates = 0
    bad = []
    for rid, seq in entries:
        if not seq or any(ch not in valid for ch in seq):
            bad.append(seq)
            continue
        if seq in spacer_to_id:
            duplicates += 1
            if spacer_to_id[seq] is None and rid is not None:
                spacer_to_id[seq] = rid
            continue
        spacer_to_id[seq] = rid

    if bad:
        raise ValueError(
            f"{len(bad)} library entr{'y' if len(bad) == 1 else 'ies'} contain non-ACGTN "
            f"characters (e.g. {bad[0][:60]!r}). If this is a TSV/CSV, select the barcode "
            f"column with --column (default 'spacer')."
        )
    if not spacer_to_id:
        raise ValueError("No barcodes loaded from the library file.")
    lengths = {len(s) for s in spacer_to_id}
    if len(lengths) != 1:
        raise ValueError(
            f"Barcodes have multiple lengths {sorted(lengths)}; all barcodes must be the "
            f"same length.\nIf this disrupts your pipeline, tell us: {ISSUES_URL}"
        )
    return spacer_to_id, duplicates


def count_near_collisions(lib_set, max_mm, bc_len, work_cap=200_000_000):
    """How many library barcodes have another library barcode within max_mm?

    Returns None when the check would be too expensive (so callers can note it was
    skipped rather than report a misleadingly small number).
    """
    if max_mm <= 0:
        return None
    # neighbors per barcode for ACGT substitutions up to max_mm
    per = 0
    for d in range(1, max_mm + 1):
        c = 1
        for i in range(d):
            c = c * (bc_len - i) // (i + 1)
        per += c * (3 ** d)
    if per * len(lib_set) > work_cap:
        return None
    colliding = 0
    for seq in lib_set:
        for nb in hamming_neighbors(seq, max_mm):
            if nb in lib_set:
                colliding += 1
                break
    return colliding


# --------------------------------------------------------------------------- #
# Flank detection (sampling)
# --------------------------------------------------------------------------- #
def best_flank(counts, max_flank):
    """Pick the longest reliable constant flank from a Counter of suffix/prefix variants.

    Descends from the longest length; stops shortening once the per-length max count
    plateaus (the constant-primer region), via a 3x gate. (Validated behaviour carried
    over from heuristicount.py's extract_best_flank.)
    """
    chosen = None
    for fl_len in range(max_flank, 0, -1):
        candidates = [s for s in counts if len(s) == fl_len]
        if not candidates:
            continue
        top = max(candidates, key=lambda x: counts[x])
        if chosen is None:
            chosen = top
        elif counts[top] > 3 * counts[chosen]:
            chosen = top
    return chosen


def consensus_flanks(pairs, bc_len, max_flank):
    """Given (read, barcode_offset) pairs, derive the constant left/right flanks."""
    L, R = Counter(), Counter()
    for read, off in pairs:
        l_ctx = read[max(0, off - max_flank):off]
        r_ctx = read[off + bc_len:off + bc_len + max_flank]
        for i in range(1, len(l_ctx) + 1):
            L[l_ctx[-i:]] += 1
        for i in range(1, len(r_ctx) + 1):
            R[r_ctx[:i]] += 1
    return best_flank(L, max_flank), best_flank(R, max_flank)


def first_hit_offset(read, bc_len, lib_set, rev_set):
    """Scan a read for the first window matching a library barcode. Returns (offset, orient)."""
    for i in range(len(read) - bc_len + 1):
        kmer = read[i:i + bc_len]
        if kmer in lib_set:
            return i, "forward"
        if kmer in rev_set:
            return i, "reverse"
    return None, None


def detect(file1, file2, lib_set, bc_len, max_flank, sample_size):
    """Sample reads to determine orientation, supported offsets, and flanks per stream."""
    rev_set = {rev_comp(b) for b in lib_set}
    paired = file2 is not None

    g1 = iter_seqs(file1)
    g2 = iter_seqs(file2) if paired else None

    f1_fwd, f1_rev, f2_fwd, f2_rev = [], [], [], []
    f1_off, f2_off = Counter(), Counter()
    seen = 0

    while seen < sample_size:
        rec1 = next(g1, None)
        if rec1 is None:
            break
        rec2 = next(g2, None) if g2 is not None else None
        if g2 is not None and rec2 is None:
            break
        seen += 1
        r1 = rec1[0]
        r2 = rec2[0] if rec2 is not None else None

        off, orient = first_hit_offset(r1, bc_len, lib_set, rev_set)
        if orient == "forward":
            f1_fwd.append((r1, off)); f1_off[off] += 1
        elif orient == "reverse":
            f1_rev.append((r1, off))

        if r2 is not None:
            off2, orient2 = first_hit_offset(r2, bc_len, lib_set, rev_set)
            if orient2 == "forward":
                f2_fwd.append((r2, off2))
            elif orient2 == "reverse":
                f2_rev.append((r2, off2)); f2_off[off2] += 1

    info(f"Sampled {seen:,} reads "
         f"(forward-read hits: {len(f1_fwd):,} fwd / {len(f1_rev):,} rev"
         + (f"; reverse-read hits: {len(f2_fwd):,} fwd / {len(f2_rev):,} rev)" if paired else ")"))

    return {
        "paired": paired,
        "f1_fwd": f1_fwd, "f1_rev": f1_rev, "f2_fwd": f2_fwd, "f2_rev": f2_rev,
    }


def supported_offsets(pairs, max_keep=6, min_frac=0.02):
    """Offsets with non-trivial support, as a fixed-offset fallback for the anchor search."""
    c = Counter(off for _, off in pairs)
    total = sum(c.values()) or 1
    offs = [off for off, n in c.most_common(max_keep) if n / total >= min_frac]
    return offs or [off for off, _ in c.most_common(1)]


# --------------------------------------------------------------------------- #
# Worker (counting pass)
# --------------------------------------------------------------------------- #
class _Cfg:
    __slots__ = (
        "lib_set", "bc_len", "min_q", "max_low", "soft", "paired", "need_swap", "single_reverse",
        "fL", "fLn", "fR", "fRn", "f_off", "f_lo", "f_hi",
        "rL", "rLn", "rR", "rRn", "r_off", "r_lo", "r_hi",
    )


_G = None


def _init_worker(cfg):
    global _G
    _G = cfg


def _candidates(seq, qual, forward):
    """Forward-oriented (barcode, quality) candidates extracted from a read.

    Anchors are searched only within the barcode-start window [lo, hi] (the offsets seen
    during sampling, plus slack), which stops a short anchor from matching a spurious
    earlier/later occurrence while still absorbing staggering. For reverse reads the
    barcode is rev-complemented and its quality reversed so both align to the library.
    """
    g = _G
    bc = g.bc_len
    if forward:
        L, Ln, R, Rn, offs, lo, hi = g.fL, g.fLn, g.fR, g.fRn, g.f_off, g.f_lo, g.f_hi
    else:
        L, Ln, R, Rn, offs, lo, hi = g.rL, g.rLn, g.rR, g.rRn, g.r_off, g.r_lo, g.r_hi

    starts = []
    if L:
        p = seq.find(L, max(0, lo - Ln), max(0, hi - Ln) + Ln)
        if p != -1:
            starts.append(p + Ln)
    if R:
        q = seq.find(R, max(0, lo + bc), max(0, hi + bc) + Rn)
        if q != -1 and q - bc >= 0:
            starts.append(q - bc)
    starts.extend(offs)

    out = []
    seen = set()
    for st in starts:
        if st < 0 or st + bc > len(seq) or st in seen:
            continue
        seen.add(st)
        b = seq[st:st + bc]
        ql = qual[st:st + bc] if qual is not None else None
        if not forward:
            b = rev_comp(b)
            ql = ql[::-1] if ql is not None else None
        out.append((b, ql))
    return out


def _classify(seq, qual, forward):
    g = _G
    saw_ambiguous = False
    for b, ql in _candidates(seq, qual, forward):
        canon, kind = assign(b, ql, g.lib_set, g.min_q, g.max_low)
        if canon is not None:
            return canon, kind
        if kind == "ambiguous":
            saw_ambiguous = True
    return None, ("ambiguous" if saw_ambiguous else "unmatched")


def score_barcode(b, ql, lib_set, min_q, max_low):
    """Library guides reachable from barcode b by varying its low-confidence positions,
    each with its likelihood score restricted to those positions.

    Returns {guide: score}. The exact match (no flips) is the score-1 member. Returns {}
    if no library guide is reachable, or if the barcode has too many low-confidence bases
    to enumerate (more than max_low), in which case only an exact match is honored.
    """
    if b in lib_set and not any(ch == "N" for ch in b):
        # fast path: clean exact hit with no low-quality wildcarding needed below
        low_clean = low_quality_positions(b, ql, min_q)
        if not low_clean:
            return {b: 1.0}
    low = low_quality_positions(b, ql, min_q)
    if not low or len(low) > max_low:
        return {b: 1.0} if b in lib_set else {}
    es = []
    for j in low:
        if b[j] == "N" or ql is None:
            es.append(0.75)  # no usable base call -> uniform over the 4 bases
        else:
            es.append(10.0 ** (-(ord(ql[j]) - 33) / 10.0))
    out = {}
    chars = list(b)
    for repl in product("ACGT", repeat=len(low)):
        for k, j in enumerate(low):
            chars[j] = repl[k]
        cand = "".join(chars)
        if cand in lib_set:
            s = 1.0
            for k, j in enumerate(low):
                e = es[k]
                s *= (1.0 - e) if repl[k] == b[j] else e / 3.0
            out[cand] = out.get(cand, 0.0) + s
    return out


def _localize_soft(seq, qual, forward):
    """{guide: score} from the first localization (anchor/offset) that yields any match."""
    g = _G
    for b, ql in _candidates(seq, qual, forward):
        sc = score_barcode(b, ql, g.lib_set, g.min_q, g.max_low)
        if sc:
            return sc
    return {}


def _bucket_key(sc):
    """Canonical key for a contested observation: sorted (guide, normalized-score) pairs.

    Scores are normalized and rounded so reads with the same likelihood profile share a
    bucket. run_em renormalizes per iteration, so only the rounded ratio matters; 6 places
    keeps that approximation well below sequencing noise while bounding the bucket count.
    """
    z = sum(sc.values()) or 1.0
    return tuple(sorted((guide, round(s / z, 6)) for guide, s in sc.items()))


def _process_soft(chunk):
    g = _G
    r1, r2 = chunk
    exact = Counter()
    contested = Counter()
    stats = Counter()

    def record(sc):
        if not sc:
            stats["unmatched"] += 1
        elif len(sc) == 1:
            (guide, _), = sc.items()
            exact[guide] += 1
            stats["exact"] += 1
        else:
            contested[_bucket_key(sc)] += 1
            stats["contested"] += 1

    if g.paired:
        fwd_reads, rev_reads = (r2, r1) if g.need_swap else (r1, r2)
        for (sa, qa), (sb, qb) in zip(fwd_reads, rev_reads):
            stats["total"] += 1
            c1 = _localize_soft(sa, qa, True)
            c2 = _localize_soft(sb, qb, False)
            if c1 and c2:
                common = set(c1) & set(c2)
                if common:
                    record({gg: c1[gg] * c2[gg] for gg in common})
                else:
                    stats["pair_disagree"] += 1
            elif c1:
                record(c1)
            elif c2:
                record(c2)
            else:
                stats["unmatched"] += 1
    else:
        forward = not g.single_reverse
        for s, q in r1:
            stats["total"] += 1
            record(_localize_soft(s, q, forward))

    return exact, contested, stats


def process_chunk(chunk):
    if _G.soft:
        return _process_soft(chunk)
    counts, stats = _process_hard(chunk)
    return counts, None, stats


def _process_hard(chunk):
    g = _G
    r1, r2 = chunk
    counts = Counter()
    stats = Counter()

    if g.paired:
        fwd_reads, rev_reads = (r2, r1) if g.need_swap else (r1, r2)
        for (sa, qa), (sb, qb) in zip(fwd_reads, rev_reads):
            stats["total"] += 1
            ca, ka = _classify(sa, qa, True)
            cb, kb = _classify(sb, qb, False)
            if ca is not None and cb is not None:
                if ca == cb:
                    counts[ca] += 1
                    stats["exact" if "exact" in (ka, kb) else "recovered"] += 1
                else:
                    stats["pair_disagree"] += 1
            elif ca is not None:
                counts[ca] += 1
                stats[ka] += 1
            elif cb is not None:
                counts[cb] += 1
                stats[kb] += 1
            elif "ambiguous" in (ka, kb):
                stats["ambiguous"] += 1
            else:
                stats["unmatched"] += 1
    else:
        forward = not g.single_reverse
        for s, q in r1:
            stats["total"] += 1
            canon, kind = _classify(s, q, forward)
            if canon is not None:
                counts[canon] += 1
            stats[kind] += 1

    return counts, stats


# --------------------------------------------------------------------------- #
# Orchestration
# --------------------------------------------------------------------------- #
MIN_ANCHOR = 4  # anchors shorter than this are not specific enough; fall back to offsets


def _anchor(flank, anchor_len, side):
    """Take the barcode-adjacent end of a flank as the anchor (L: suffix, R: prefix).

    Returns None for flanks too short to anchor reliably; localization then relies on
    the supported-offset fallback.
    """
    if not flank:
        return None
    a = flank[-anchor_len:] if side == "L" else flank[:anchor_len]
    return a if len(a) >= MIN_ANCHOR else None


def build_config(args, spacer_to_id, bc_len, det):
    paired = det["paired"]
    cfg = _Cfg()
    cfg.lib_set = set(spacer_to_id)
    cfg.bc_len = bc_len
    cfg.min_q = args.min_quality
    cfg.max_low = args.mismatches
    cfg.soft = args.soft
    cfg.paired = paired
    cfg.need_swap = False
    cfg.single_reverse = False
    cfg.fL = cfg.fR = cfg.rL = cfg.rR = None
    cfg.fLn = cfg.fRn = cfg.rLn = cfg.rRn = 0
    cfg.f_off = []
    cfg.r_off = []
    cfg.f_lo = cfg.f_hi = cfg.r_lo = cfg.r_hi = 0

    al = args.anchor_length
    mf = args.max_flank

    if not paired:
        if len(det["f1_fwd"]) >= len(det["f1_rev"]):
            pairs = det["f1_fwd"]
            cfg.single_reverse = False
            L, R = consensus_flanks(pairs, bc_len, mf)
            cfg.fL, cfg.fR = _anchor(L, al, "L"), _anchor(R, al, "R")
            cfg.f_off = supported_offsets(pairs)
            fwd_flanks, rev_flanks = (L, R), (None, None)
        else:
            pairs = det["f1_rev"]
            cfg.single_reverse = True
            L, R = consensus_flanks(pairs, bc_len, mf)
            cfg.rL, cfg.rR = _anchor(L, al, "L"), _anchor(R, al, "R")
            cfg.r_off = supported_offsets(pairs)
            fwd_flanks, rev_flanks = (None, None), (L, R)
    else:
        f1_forward = len(det["f1_fwd"]) >= len(det["f1_rev"])
        cfg.need_swap = not f1_forward
        fwd_pairs = det["f1_fwd"] if f1_forward else det["f2_fwd"]
        rev_pairs = det["f2_rev"] if f1_forward else det["f1_rev"]
        Lf, Rf = consensus_flanks(fwd_pairs, bc_len, mf)
        Lr, Rr = consensus_flanks(rev_pairs, bc_len, mf)
        cfg.fL, cfg.fR = _anchor(Lf, al, "L"), _anchor(Rf, al, "R")
        cfg.rL, cfg.rR = _anchor(Lr, al, "L"), _anchor(Rr, al, "R")
        cfg.f_off = supported_offsets(fwd_pairs)
        cfg.r_off = supported_offsets(rev_pairs)
        fwd_flanks, rev_flanks = (Lf, Rf), (Lr, Rr)

    cfg.fLn = len(cfg.fL) if cfg.fL else 0
    cfg.fRn = len(cfg.fR) if cfg.fR else 0
    cfg.rLn = len(cfg.rL) if cfg.rL else 0
    cfg.rRn = len(cfg.rR) if cfg.rR else 0

    slack = max(3, args.anchor_length // 2)

    def _window(offs):
        if not offs:
            return 0, 0
        return max(0, min(offs) - slack), max(offs) + slack

    cfg.f_lo, cfg.f_hi = _window(cfg.f_off)
    cfg.r_lo, cfg.r_hi = _window(cfg.r_off)

    if not any([cfg.fL, cfg.fR, cfg.rL, cfg.rR, cfg.f_off, cfg.r_off]):
        raise ValueError(
            "Could not localize barcodes from the sample (no flanks or offsets found). "
            "Are the reads from the same library as the barcodes?"
        )

    fwd_active = (cfg.paired) or (not cfg.single_reverse)
    rev_active = (cfg.paired) or cfg.single_reverse
    if fwd_active and not (cfg.fL or cfg.fR) and len(cfg.f_off) > 1:
        warn("No reliable forward flank found; localizing by offset across "
             f"{len(cfg.f_off)} positions - staggered-primer reads may be undercounted.")
    if rev_active and not (cfg.rL or cfg.rR) and len(cfg.r_off) > 1:
        warn("No reliable reverse flank found; localizing by offset across "
             f"{len(cfg.r_off)} positions - staggered-primer reads may be undercounted.")
    return cfg, fwd_flanks, rev_flanks


def run_em(exact, contested, guides, alpha, max_iter=500, tol=1e-9):
    """Apportion contested reads across guides by EM.

    exact:     Counter of whole-count (unambiguous) reads per guide.
    contested: Counter of canonical bucket keys -> multiplicity; each key is a tuple of
               (guide, normalized-score) pairs.
    guides:    full set of library guides (so zero-count guides stay in the output).
    alpha:     symmetric Dirichlet pseudocount used inside the E-step, which keeps the
               estimate interior (so a guide seen only in contested reads is not stuck at
               the spurious pi=0 fixed point).

    Returns (counts dict guide->float, info dict).
    """
    counts = Counter({g: float(exact.get(g, 0)) for g in guides})
    if not contested:
        return counts, {"iters": 0, "contested_reads": 0}

    buckets = []
    involved = set()
    for key, mult in contested.items():
        gs = [guide for guide, _ in key]
        ss = [s for _, s in key]
        buckets.append((gs, ss, mult))
        involved.update(gs)
    contested_reads = sum(m for _, _, m in buckets)

    iters = 0
    delta = 0.0
    for iters in range(1, max_iter + 1):
        contrib = defaultdict(float)
        for gs, ss, mult in buckets:
            ws = [(counts.get(gg, 0.0) + alpha) * s for gg, s in zip(gs, ss)]
            z = sum(ws)
            if z <= 0.0:
                for gg in gs:
                    contrib[gg] += mult / len(gs)
                continue
            for gg, w in zip(gs, ws):
                contrib[gg] += mult * w / z
        delta = 0.0
        for gg in involved:
            new = float(exact.get(gg, 0)) + contrib.get(gg, 0.0)
            delta = max(delta, abs(new - counts[gg]))
            counts[gg] = new
        if delta < tol:
            break

    if iters >= max_iter and delta >= tol:
        warn(f"EM did not converge after {max_iter} iterations (delta={delta:.2e}); "
             f"using current estimates.")
    return counts, {"iters": iters, "contested_reads": contested_reads}


def emit_table(args, spacer_to_id, counts, out):
    name = args.name or "count"
    has_id = args.id_column is not None

    def fmt(seq):
        c = counts.get(seq, 0)
        rc = round(c)
        return str(int(rc)) if abs(c - rc) < 1e-9 else f"{c:.4f}"

    if has_id:
        out.write(f"id\tbarcode\t{name}\n")
        for seq in sorted(spacer_to_id):
            out.write(f"{spacer_to_id[seq] or ''}\t{seq}\t{fmt(seq)}\n")
    else:
        out.write(f"barcode\t{name}\n")
        for seq in sorted(spacer_to_id):
            out.write(f"{seq}\t{fmt(seq)}\n")


def summary_table(args, bc_len, cfg, fwd_flanks, rev_flanks, spacer_to_id,
                  counts, stats, total_reads, duplicates, near_collisions, em_info=None):
    t = Table(
        box=rich.table.box.SIMPLE_HEAVY,
        caption=f"Finished at [u]{datetime.now()}[/u]",
        header_style="bold bright_white",
        border_style="bold bright_white",
        show_header=True,
        expand=False,
    )
    t.add_column(os.path.basename(sys.argv[0]), justify="right", style="white", no_wrap=True)
    t.add_column("Summary", justify="right", no_wrap=True)

    def row(k, v):
        t.add_row(k, str(v))

    t.add_section()
    row("[bold bright_magenta]Input & Config[/]", "")
    row("Barcode library", os.path.basename(args.library))
    row("Forward reads", os.path.basename(args.file1))
    if args.file2:
        row("Reverse reads", os.path.basename(args.file2))
    row("Mode", "soft (EM expected counts)" if args.soft else "hard (quality-gated)")
    if args.min_quality > 0:
        row("Min base quality (Phred)", args.min_quality)
        row("Low-Q bases forgiven", args.mismatches)
    else:
        row("Matching", "exact only (-q 0)")
    if args.soft:
        row("EM prior (pseudocount)", args.em_prior)
    row("Threads", args.threads)
    row("Operating system", platform.system())

    t.add_section()
    row("[bold bright_blue]Heuristics[/]", "")
    row("Barcode length", bc_len)
    if cfg.paired:
        row("Read1 orientation", "reverse (swapped)" if cfg.need_swap else "forward")
    elif cfg.single_reverse:
        row("Orientation", "reverse")
    else:
        row("Orientation", "forward")
    if cfg.fL or cfg.fR:
        row("Forward flanks", f"{cfg.fL or ''}...{cfg.fR or ''}")
    if cfg.rL or cfg.rR:
        row("Reverse flanks", f"{cfg.rL or ''}...{cfg.rR or ''}")
    if cfg.f_off:
        row("Forward offsets", ",".join(map(str, cfg.f_off)))
    if cfg.r_off:
        row("Reverse offsets", ",".join(map(str, cfg.r_off)))

    t.add_section()
    row("[bold bright_green]Library[/]", "")
    row("Library barcodes", f"{len(spacer_to_id):,}")
    if duplicates:
        row("Duplicate spacers collapsed", f"{duplicates:,}")
    if near_collisions:
        row("Barcodes 1 mismatch from another", f"{near_collisions:,}")

    assigned = sum(counts.values())
    seen = sum(1 for v in counts.values() if round(v) > 0)
    t.add_section()
    row("[bold]Total reads[/]", f"{total_reads:,}")
    row("Assigned reads", f"{assigned:,.0f}")
    if args.soft:
        row("  whole (unambiguous)", f"{stats.get('exact', 0):,}")
        row("  contested (EM-split)", f"{stats.get('contested', 0):,}")
        if em_info:
            row("  EM iterations", em_info.get("iters", 0))
    else:
        row("  exact", f"{stats.get('exact', 0):,}")
        row("  recovered (low-Q wildcard)", f"{stats.get('recovered', 0):,}")
        row("Ambiguous (skipped)", f"{stats.get('ambiguous', 0):,}")
    if cfg.paired:
        row("Pair disagreed (skipped)", f"{stats.get('pair_disagree', 0):,}")
    row("Unmatched (noise)", f"{stats.get('unmatched', 0):,}")
    frac = assigned / total_reads if total_reads else 0
    row("[bold]Assigned fraction[/]", f"{frac:.3f}")

    t.add_section()
    row("Guides detected", f"{seen:,} / {len(spacer_to_id):,}")
    row("Guides with zero reads", f"{len(spacer_to_id) - seen:,}")

    t.add_section()
    top = sorted(counts.items(), key=lambda kv: kv[1], reverse=True)[:5]
    row("[bold bright_green]Top barcodes[/]", "")
    for bc, n in top:
        row(bc, f"{n:,.4g}" if args.soft else f"{n:,}")

    _console.print(t)


def run(args):
    info(f"heuristicount2 v{__version__} - initializing...")

    info("Reading barcode library...")
    spacer_to_id, duplicates = load_library(args.library, args.column, args.id_column)
    bc_len = len(next(iter(spacer_to_id)))
    lib_set = set(spacer_to_id)
    info(f"Loaded {len(spacer_to_id):,} unique barcodes of length {bc_len}"
         + (f" ({duplicates:,} duplicate spacers collapsed)" if duplicates else "") + ".")
    if len(spacer_to_id) < 10:
        warn("Fewer than 10 barcodes loaded - detection may be unreliable.")

    if args.min_quality > 0:
        info(f"Quality-gated matching: bases at Phred >= {args.min_quality} must match exactly; "
             f"up to {args.mismatches} sub-threshold base(s) per barcode may be forgiven.")
        if not _strip_compression(args.file1).endswith((".fastq", ".fq")):
            warn("Reads carry no quality scores (.reads): only 'N' bases can be forgiven; "
                 "all other positions require an exact match.")
    else:
        info("Exact matching (--min-quality 0).")

    near_collisions = None
    if args.min_quality > 0 and args.mismatches > 0:
        near_collisions = count_near_collisions(lib_set, 1, bc_len)
        if near_collisions:
            warn(f"{near_collisions:,} library barcode(s) differ from another by a single base "
                 f"(e.g. a designed mismatch series). A read is skipped as ambiguous only when "
                 f"that differing base is itself low-confidence - designed mismatches called at "
                 f"Phred >= {args.min_quality} stay distinct.")

    info("Sampling reads to detect orientation, offsets, and flanks...")
    det = detect(args.file1, args.file2, lib_set, bc_len, args.max_flank, args.sample_size)
    cfg, fwd_flanks, rev_flanks = build_config(args, spacer_to_id, bc_len, det)
    info(f"Forward flanks {cfg.fL or '-'}...{cfg.fR or '-'}"
         + (f" | reverse flanks {cfg.rL or '-'}...{cfg.rR or '-'}" if (cfg.rL or cfg.rR) else "")
         + ".")

    info(f"Counting with {args.threads} worker(s)" + (" [soft/EM]" if args.soft else "") + "...")
    stats = Counter()
    exact = Counter()        # soft: whole-count reads ; hard: unused
    contested = Counter()    # soft: contested buckets
    counts = Counter()       # hard: integer counts
    chunks = iter_chunks(args.file1, args.file2, args.chunk_size)

    def absorb(primary, secondary, s):
        stats.update(s)
        if args.soft:
            exact.update(primary)
            contested.update(secondary)
        else:
            counts.update(primary)

    if args.threads <= 1:
        _init_worker(cfg)
        for chunk in chunks:
            absorb(*process_chunk(chunk))
    else:
        with Pool(args.threads, initializer=_init_worker, initargs=(cfg,)) as pool:
            for res in pool.imap_unordered(process_chunk, chunks, chunksize=1):
                absorb(*res)

    total_reads = stats.get("total", 0)
    info("Collating results...")

    em_info = None
    if args.soft:
        counts, em_info = run_em(exact, contested, lib_set, args.em_prior)
        info(f"EM apportioned {em_info['contested_reads']:,} contested read(s) "
             f"in {em_info['iters']} iteration(s).")

    out_ctx = open(args.output, "w") if args.output else nullcontext(sys.stdout)
    with out_ctx as out:
        emit_table(args, spacer_to_id, counts, out)

    summary_table(args, bc_len, cfg, fwd_flanks, rev_flanks, spacer_to_id,
                  counts, stats, total_reads, duplicates, near_collisions, em_info)


def build_parser():
    p = argparse.ArgumentParser(
        prog="heuristicount2.py",
        description="Flank-anchored, quality-gated barcode/sgRNA counter for pooled screens.",
    )
    p.add_argument("library", help="Barcode library: .tsv/.csv (use --column), FASTA, or one-per-line.")
    p.add_argument("file1", help="Forward reads: FASTQ or .reads (optionally .gz/.zst).")
    p.add_argument("file2", nargs="?", default=None, help="Reverse reads for paired-end (optional).")
    p.add_argument("--column", default="spacer",
                   help="Column holding the barcode sequence in a .tsv/.csv library (default: spacer).")
    p.add_argument("--id-column", default=None,
                   help="Optional column holding a guide ID; adds an 'id' column to the output.")
    p.add_argument("-q", "--min-quality", type=int, default=20,
                   help="Phred cutoff: bases at or above this must match a guide EXACTLY; only "
                        "sub-threshold (low-confidence) bases may be forgiven. 0 = exact matching "
                        "with no forgiveness (default: 20). Keeps designed-mismatch guides distinct.")
    p.add_argument("-m", "--mismatches", type=int, default=2, metavar="MAX_LOW_Q",
                   help="Max number of low-confidence (sub --min-quality) bases that may be "
                        "wildcarded per barcode; reads with more are dropped as unreliable "
                        "(default: 2). High-confidence bases are never forgiven.")
    p.add_argument("--soft", action="store_true",
                   help="Emit EM expected counts: each read with two or more candidate guides "
                        "is apportioned by posterior probability rather than dropped as "
                        "ambiguous. Counts are generally fractional. See heuristicount2_model.md.")
    p.add_argument("--em-prior", type=float, default=0.5, metavar="ALPHA",
                   help="Symmetric Dirichlet pseudocount used inside the soft-mode EM step "
                        "(default: 0.5); keeps the estimate interior so a guide seen only in "
                        "contested reads is not pinned at zero. Only used with --soft.")
    p.add_argument("-t", "--threads", type=int, default=None,
                   help="Worker processes (default: half of available CPUs).")
    p.add_argument("-o", "--output", default=None,
                   help="Write the count table here instead of stdout.")
    p.add_argument("--name", default="count",
                   help="Header name for the count column (e.g. the sample name).")
    p.add_argument("--anchor-length", type=int, default=8,
                   help="Length of the flank anchor used to localize the barcode (default: 8).")
    p.add_argument("--max-flank", type=int, default=12,
                   help="Maximum flank length to detect (default: 12).")
    p.add_argument("--sample-size", type=int, default=100_000,
                   help="Reads to sample for detection (default: 100000).")
    p.add_argument("--chunk-size", type=int, default=65_536,
                   help="Reads per work chunk (default: 65536).")
    p.add_argument("--debug", action="store_true", help="Re-raise tracebacks instead of a tidy error.")
    p.add_argument("--version", action="version", version=f"heuristicount2 {__version__}")
    return p


def main(argv=None):
    args = build_parser().parse_args(argv)

    if args.threads is None:
        try:
            args.threads = max(1, len(os.sched_getaffinity(0)) // 2)
        except AttributeError:
            args.threads = max(1, (os.cpu_count() or 2) // 2)
    if args.mismatches < 0:
        error("--mismatches (max low-quality bases) must be >= 0.")
        return 2
    if args.min_quality < 0:
        error("--min-quality must be >= 0.")
        return 2
    if args.anchor_length < MIN_ANCHOR:
        error(f"--anchor-length must be >= {MIN_ANCHOR} (shorter anchors are not specific enough).")
        return 2
    if args.sample_size < 1:
        error("--sample-size must be >= 1.")
        return 2
    if args.chunk_size < 1:
        error("--chunk-size must be >= 1.")
        return 2
    if args.em_prior < 0:
        error("--em-prior must be >= 0.")
        return 2

    try:
        run(args)
        return 0
    except ValueError as ve:
        if args.debug:
            raise
        error(str(ve))
        return 2
    except FileNotFoundError as fe:
        if args.debug:
            raise
        error(str(fe))
        return 2
    except Exception as e:  # noqa: BLE001
        if args.debug:
            raise
        error(f"Unexpected error: {e}")
        return 1


if __name__ == "__main__":
    sys.exit(main())
