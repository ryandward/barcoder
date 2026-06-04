#!/usr/bin/env python3
"""Regression tests for heuristicount2.py (quality-gated matching model).

Runnable two ways:
    python test_heuristicount2.py     # prints PASS/FAIL, exits non-zero on failure
    pytest test_heuristicount2.py     # standard pytest discovery (test_* functions)

Self-contained: synthetic libraries and reads are generated into a temp dir.
Integration tests invoke the CLI as a subprocess; unit tests import directly.

Matching model under test: a base called at Phred >= --min-quality must match a guide
EXACTLY; only sub-threshold bases may be wildcarded. So a *low-quality* sequencing error
is recovered, but a *high-quality* difference is NOT forgiven (this is what keeps a
designed single-mismatch guide from absorbing reads of its perfect-match sibling).
"""

import os
import random
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPT = os.path.join(HERE, "heuristicount2.py")

BASES = "ACGT"
LAD = "GTACGTACGT"   # 10 bp left constant  -> barcode starts at offset 10
RAD = "TGCATGCATG"   # 10 bp right constant
BCLEN = 20
NBC = 20
PER = 100

HIQ = chr(40 + 33)   # 'I'  Phred 40 (trusted)
LOQ = chr(2 + 33)    # '#'  Phred 2  (low confidence)


def _rc(s):
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def _rnd(rng, n):
    return "".join(rng.choice(BASES) for _ in range(n))


def _sub_at(rng, seq, i):
    alt = rng.choice([b for b in BASES if b != seq[i]])
    return seq[:i] + alt + seq[i + 1:]


def _write_fastq(path, records):
    """records: iterable of (seq, qual) with equal lengths."""
    with open(path, "w") as f:
        for i, (seq, qual) in enumerate(records):
            f.write(f"@r{i}\n{seq}\n+\n{qual}\n")


def _read(seq):
    """Wrap a barcode body in constant flanks; all bases high-quality."""
    full = LAD + seq + RAD
    return full, HIQ * len(full)


def _make_data(d):
    rng = random.Random(42)
    barcodes = set()
    while len(barcodes) < NBC:
        barcodes.add(_rnd(rng, BCLEN))
    barcodes = sorted(barcodes)

    p = {"barcodes": barcodes}
    p["lib"] = os.path.join(d, "lib.fasta")
    with open(p["lib"], "w") as f:
        f.write("\n".join(barcodes) + "\n")

    nomatch = set()
    while len(nomatch) < NBC:
        c = _rnd(rng, BCLEN)
        if c not in barcodes:
            nomatch.add(c)
    p["lib_nomatch"] = os.path.join(d, "lib_nomatch.fasta")
    with open(p["lib_nomatch"], "w") as f:
        f.write("\n".join(sorted(nomatch)) + "\n")

    # clean: all high quality, exact
    p["clean"] = os.path.join(d, "clean.fastq")
    _write_fastq(p["clean"], [_read(bc) for bc in barcodes for _ in range(PER)])

    # 50% of reads carry a single substitution whose base is LOW quality
    def lowq_err():
        for bc in barcodes:
            for _ in range(PER):
                if rng.random() < 0.5:
                    i = rng.randrange(BCLEN)
                    body = _sub_at(rng, bc, i)
                    full = LAD + body + RAD
                    qual = list(HIQ * len(full))
                    qual[len(LAD) + i] = LOQ          # error base flagged low-confidence
                    yield full, "".join(qual)
                else:
                    yield _read(bc)
    p["mut_lowq"] = os.path.join(d, "mut_lowq.fastq")
    _write_fastq(p["mut_lowq"], list(lowq_err()))

    # 50% of reads carry a single substitution but the base is HIGH quality
    def highq_err():
        for bc in barcodes:
            for _ in range(PER):
                if rng.random() < 0.5:
                    i = rng.randrange(BCLEN)
                    full = LAD + _sub_at(rng, bc, i) + RAD
                    yield full, HIQ * len(full)        # error base looks trustworthy
                else:
                    yield _read(bc)
    p["mut_highq"] = os.path.join(d, "mut_highq.fastq")
    _write_fastq(p["mut_highq"], list(highq_err()))

    # staggered primer (0-4 bp pad before the constant region), all high quality
    def stag():
        for bc in barcodes:
            for _ in range(PER):
                pad = _rnd(rng, rng.randint(0, 4))
                full = pad + LAD + bc + RAD
                yield full, HIQ * len(full)
    p["stag"] = os.path.join(d, "stag.fastq")
    _write_fastq(p["stag"], list(stag()))

    # paired-end (clean)
    p["r1"] = os.path.join(d, "R1.fastq")
    p["r2"] = os.path.join(d, "R2.fastq")
    amps = [LAD + bc + RAD for bc in barcodes for _ in range(PER)]
    _write_fastq(p["r1"], [(a, HIQ * len(a)) for a in amps])
    _write_fastq(p["r2"], [(_rc(a), HIQ * len(a)) for a in amps])

    # reverse single-end (clean)
    p["rev"] = os.path.join(d, "rev.fastq")
    _write_fastq(p["rev"], [(_rc(a), HIQ * len(a)) for a in amps])
    return p


def _run(*args):
    r = subprocess.run([sys.executable, SCRIPT, *args], capture_output=True, text=True)
    return r.returncode, r.stdout, r.stderr


def _counts(stdout):
    out = {}
    for line in stdout.splitlines():
        if not line or line.startswith(("barcode\t", "id\t")):
            continue
        parts = line.split("\t")
        out[parts[0]] = float(parts[-1])
    return out


def _designed_pair_data(d):
    """Library with a designed single-mismatch pair X, Y (differ at barcode position 0)
    plus 18 distinct decoys. Reads: 1000 clean X (high quality), 100 X molecules whose
    discriminating base ERRORED to Y's base with LOW quality (the 'evil twin' case), and
    50 clean reads per decoy. True abundance of Y is zero."""
    X = "A" * BCLEN
    Y = "C" + "A" * (BCLEN - 1)
    rng = random.Random(5)
    others = set()
    while len(others) < 18:
        c = _rnd(rng, BCLEN)
        if c not in (X, Y):
            others.add(c)
    lib = [X, Y] + sorted(others)
    libp = os.path.join(d, "pair_lib.fasta")
    with open(libp, "w") as f:
        f.write("\n".join(lib) + "\n")

    recs = [_read(X) for _ in range(1000)]                 # clean X, all high quality
    for _ in range(100):                                   # X miscalled to Y's base, low Q there
        full = LAD + Y + RAD                               # sequence reads as Y...
        q = list(HIQ * len(full))
        q[len(LAD)] = LOQ                                  # ...but the deciding base is low-confidence
        recs.append((full, "".join(q)))
    for g in others:
        recs += [_read(g) for _ in range(50)]
    rp = os.path.join(d, "pair_reads.fastq")
    _write_fastq(rp, recs)
    return libp, rp, X, Y


# --------------------------------------------------------------------------- #
# Integration tests
# --------------------------------------------------------------------------- #
def test_clean_exact():
    with tempfile.TemporaryDirectory() as d:
        p = _make_data(d)
        rc, out, _ = _run(p["lib"], p["clean"], "-t", "2")
        assert rc == 0
        c = _counts(out)
        assert len(c) == NBC and sum(c.values()) == NBC * PER
        assert all(v == PER for v in c.values())


def test_low_quality_error_recovered():
    with tempfile.TemporaryDirectory() as d:
        p = _make_data(d)
        rc, out, _ = _run(p["lib"], p["mut_lowq"], "-q", "20", "-t", "2")
        assert rc == 0
        c = _counts(out)
        # every read (including the ~50% with a *low-quality* substitution) recovered
        assert sum(c.values()) == NBC * PER
        assert all(v == PER for v in c.values())


def test_high_quality_error_NOT_forgiven():
    """The key intentional-mismatch safety property."""
    with tempfile.TemporaryDirectory() as d:
        p = _make_data(d)
        rc, out, _ = _run(p["lib"], p["mut_highq"], "-q", "20", "-t", "2")
        assert rc == 0
        total = sum(_counts(out).values())
        # only the ~50% unmutated reads match; high-quality differences are dropped
        assert NBC * PER * 0.4 < total < NBC * PER * 0.6


def test_exact_mode_q0():
    with tempfile.TemporaryDirectory() as d:
        p = _make_data(d)
        # with -q 0, even the low-quality errors are not forgiven -> ~50%
        rc, out, _ = _run(p["lib"], p["mut_lowq"], "-q", "0", "-t", "2")
        assert rc == 0
        total = sum(_counts(out).values())
        assert NBC * PER * 0.4 < total < NBC * PER * 0.6


def test_staggered_recovered():
    with tempfile.TemporaryDirectory() as d:
        p = _make_data(d)
        rc, out, _ = _run(p["lib"], p["stag"], "-t", "2")
        assert rc == 0
        assert sum(_counts(out).values()) == NBC * PER


def test_nomatch_graceful():
    with tempfile.TemporaryDirectory() as d:
        p = _make_data(d)
        rc, _, err = _run(p["lib_nomatch"], p["clean"])
        assert rc == 2 and "localize" in err.lower()


def test_paired_end():
    with tempfile.TemporaryDirectory() as d:
        p = _make_data(d)
        rc, out, _ = _run(p["lib"], p["r1"], p["r2"], "-t", "2")
        assert rc == 0
        c = _counts(out)
        assert sum(c.values()) == NBC * PER and all(v == PER for v in c.values())


def test_reverse_single_end():
    with tempfile.TemporaryDirectory() as d:
        p = _make_data(d)
        rc, out, err = _run(p["lib"], p["rev"], "-t", "2")
        assert rc == 0
        assert sum(_counts(out).values()) == NBC * PER
        assert "reverse" in err.lower()


def test_zero_count_rows_present():
    with tempfile.TemporaryDirectory() as d:
        p = _make_data(d)
        extra = "AAAAAAAAAACCCCCCCCCC"
        lib2 = os.path.join(d, "lib2.fasta")
        with open(lib2, "w") as f:
            f.write("\n".join(p["barcodes"]) + "\n" + extra + "\n")
        rc, out, _ = _run(lib2, p["clean"], "-t", "2")
        assert rc == 0
        c = _counts(out)
        assert extra in c and c[extra] == 0
        assert len(c) == NBC + 1


def test_malformed_fastq_errors():
    with tempfile.TemporaryDirectory() as d:
        p = _make_data(d)
        bad = os.path.join(d, "bad.fastq")
        with open(bad, "w") as f:
            f.write("@r0\n" + LAD + p["barcodes"][0] + RAD + "\n+\n" + HIQ * 40 + "\n")
            f.write("NOT_A_HEADER\nACGT\n+\nIIII\n")
        rc, _, err = _run(p["lib"], bad)
        assert rc == 2 and "malformed" in err.lower()


def test_unequal_paired_files_warns_not_crashes():
    with tempfile.TemporaryDirectory() as d:
        p = _make_data(d)
        short_r2 = os.path.join(d, "R2_short.fastq")
        with open(p["r2"]) as src, open(short_r2, "w") as dst:
            dst.write("".join(src.readlines()[: 4 * 1500]))
        rc, out, err = _run(p["lib"], p["r1"], short_r2, "-t", "2")
        assert rc == 0
        assert "unequal" in err.lower()
        assert sum(_counts(out).values()) == 1500


def test_tsv_column_loading():
    with tempfile.TemporaryDirectory() as d:
        p = _make_data(d)
        tsv = os.path.join(d, "lib.tsv")
        with open(tsv, "w") as f:
            f.write("spacer\tgene\n")
            for i, b in enumerate(p["barcodes"]):
                f.write(f"{b}\tgene_{i:02d}\n")
        rc, out, _ = _run(tsv, p["clean"], "--column", "spacer", "--id-column", "gene", "-t", "2")
        assert rc == 0
        assert out.splitlines()[0] == "id\tbarcode\tcount"
        assert sum(int(l.split("\t")[-1]) for l in out.splitlines()[1:]) == NBC * PER


# --------------------------------------------------------------------------- #
# Soft (EM) mode
# --------------------------------------------------------------------------- #
def test_soft_clean_equals_hard():
    with tempfile.TemporaryDirectory() as d:
        p = _make_data(d)
        rc, out, _ = _run(p["lib"], p["clean"], "--soft", "-t", "2")
        assert rc == 0
        c = _counts(out)
        assert sum(c.values()) == NBC * PER and all(v == PER for v in c.values())


def test_soft_recovers_contested_to_abundant_guide():
    """Soft mode recovers the 100 low-quality 'evil twin' reads, apportioning them to the
    abundant true source X rather than the absent designed-mismatch sibling Y."""
    with tempfile.TemporaryDirectory() as d:
        lib, reads, X, Y = _designed_pair_data(d)
        rc, out, _ = _run(lib, reads, "--soft", "-q", "20", "-t", "2")
        assert rc == 0
        c = _counts(out)
        assert 1095 < c[X] <= 1101, c[X]      # 1000 clean + ~100 recovered
        assert c[Y] < 5, c[Y]                  # not dumped on the absent sibling


def test_hard_miscounts_contested_to_sibling():
    """Contrast: hard mode trusts the exact (but low-quality) base call and miscounts all
    100 reads as the wrong guide Y. This is what soft mode fixes."""
    with tempfile.TemporaryDirectory() as d:
        lib, reads, X, Y = _designed_pair_data(d)
        rc, out, _ = _run(lib, reads, "-q", "20", "-t", "2")   # hard
        assert rc == 0
        c = _counts(out)
        assert c[X] == 1000
        assert c[Y] == 100        # 100 X-molecules wrongly banked as Y by the exact short-circuit


# --------------------------------------------------------------------------- #
# Unit tests
# --------------------------------------------------------------------------- #
def _import():
    sys.path.insert(0, HERE)
    import heuristicount2 as h
    return h


def test_rev_comp():
    h = _import()
    assert h.rev_comp("ACGTN") == "NACGT"


def test_assign_quality_gated():
    h = _import()
    G = "AAAACCCCGGGGTTTTACGT"          # 20 bp guide
    lib = {G, "TTTTAAAACCCCGGGGACGT", "CCCCGGGGTTTTAAAAACGT"}
    Q = lambda phreds: "".join(chr(p + 33) for p in phreds)
    allhi = Q([40] * len(G))

    # exact, regardless of quality
    assert h.assign(G, allhi, lib, 20, 2) == (G, "exact")

    # single error at a LOW-quality base -> recovered to G
    err = G[:5] + ("T" if G[5] != "T" else "A") + G[6:]
    qlo = list(allhi); qlo[5] = chr(2 + 33)
    assert h.assign(err, "".join(qlo), lib, 20, 2) == (G, "recovered")

    # same single difference but HIGH quality -> NOT forgiven
    assert h.assign(err, allhi, lib, 20, 2) == (None, "unmatched")

    # quality gating off (min_q=0) -> exact only
    assert h.assign(err, "".join(qlo), lib, 0, 2) == (None, "unmatched")


def test_assign_ambiguous_on_low_quality_discriminator():
    h = _import()
    # two guides differ only at position 0 (a designed single-mismatch pair)
    X = "A" * 20
    Y = "C" + "A" * 19
    lib = {X, Y, "G" * 20, "T" * 20}
    allhi = "".join(chr(40 + 33) for _ in range(20))

    # exact reads are unambiguous
    assert h.assign(X, allhi, lib, 20, 2) == (X, "exact")
    assert h.assign(Y, allhi, lib, 20, 2) == (Y, "exact")

    # a read with a *third* base at the discriminating position, called LOW quality,
    # is within reach of BOTH X and Y once that base is wildcarded -> ambiguous
    read = "G" + "A" * 19
    assert read not in lib
    qlo = list(allhi); qlo[0] = chr(2 + 33)
    assert h.assign(read, "".join(qlo), lib, 20, 2) == (None, "ambiguous")
    # ...but the same read with that base called HIGH quality is a trusted, real
    # difference from every guide -> unmatched (never silently merged into X or Y)
    assert h.assign(read, allhi, lib, 20, 2) == (None, "unmatched")


def test_reverse_strand_quality_alignment():
    """A low-quality error on the reverse strand must be recovered, which only works
    if the quality string is reversed to stay aligned with the rev-complemented barcode."""
    h = _import()
    G, bc = "AAAACCCG", 8
    lib = {G, "TTTTAAAA", "CCCCGGGG", "GGGGTTTT"}
    LF, RF = "CACACACA", "GAGAGAGA"
    fwd_amp = LF + G + RF
    rev_read = h.rev_comp(fwd_amp)                  # barcode region = rev_comp(G) at index 8
    j = 8 + 3                                        # corrupt one base inside the rev barcode
    rev_read = rev_read[:j] + ("A" if rev_read[j] != "A" else "T") + rev_read[j + 1:]
    qual = list(chr(40 + 33) for _ in rev_read)
    qual[j] = chr(2 + 33)                            # flag that base low-confidence
    qual = "".join(qual)

    cfg = h._Cfg()
    cfg.lib_set, cfg.bc_len, cfg.min_q, cfg.max_low = lib, bc, 20, 2
    cfg.paired = cfg.need_swap = False
    cfg.single_reverse = True
    cfg.fL = cfg.fR = None
    cfg.fLn = cfg.fRn = 0
    cfg.f_off, cfg.f_lo, cfg.f_hi = [], 0, 0
    cfg.rL, cfg.rLn = h.rev_comp(RF), len(RF)        # left flank of the reverse read
    cfg.rR, cfg.rRn = h.rev_comp(LF), len(LF)
    cfg.r_off, cfg.r_lo, cfg.r_hi = [8], 5, 11
    h._init_worker(cfg)

    canon, kind = h._classify(rev_read, qual, forward=False)
    assert canon == G and kind == "recovered", f"got {canon!r}/{kind}"

    # if the SAME corrupted base were high-quality, it must NOT be forgiven
    canon2, kind2 = h._classify(rev_read, chr(40 + 33) * len(rev_read), forward=False)
    assert canon2 is None, f"high-quality reverse error wrongly forgiven: {canon2!r}"


def test_n_base_treated_as_uncertain():
    h = _import()
    G = "AAAACCCCGGGGTTTTACGT"
    lib = {G, "TTTTAAAACCCCGGGGACGT"}
    with_n = "N" + G[1:]
    allhi = chr(40 + 33) * len(G)
    # an N is uncertain regardless of its quality char -> wildcarded -> recovered
    assert h.assign(with_n, allhi, lib, 20, 2) == (G, "recovered")
    # ...but in exact mode (min_q=0) an N cannot match -> unmatched
    assert h.assign(with_n, allhi, lib, 0, 2) == (None, "unmatched")


def test_score_barcode_unit():
    h = _import()
    X = "A" * 20
    Y = "C" + "A" * 19
    lib = {X, Y, "G" * 20, "T" * 20}
    allhi = chr(40 + 33) * 20
    # clean exact -> single candidate, score 1
    assert h.score_barcode(X, allhi, lib, 20, 2) == {X: 1.0}
    # low-Q at the discriminating position -> contested {X, Y}, X favored (call matches X)
    qlo = list(allhi); qlo[0] = chr(2 + 33)
    sc = h.score_barcode(X, "".join(qlo), lib, 20, 2)
    assert set(sc) == {X, Y} and sc[X] > sc[Y] > 0
    # too many low-Q bases to enumerate -> fall back to exact-only
    q5 = list(allhi)
    for i in range(5):
        q5[i] = chr(2 + 33)
    assert h.score_barcode(X, "".join(q5), lib, 20, 2) == {X: 1.0}


def test_run_em_unit():
    h = _import()
    X = "A" * 20
    Y = "C" + "A" * 19
    exact = {X: 1000}
    sc = {X: 0.21, Y: 0.37}          # one read reading as Y but low-Q at the discriminator
    contested = {h._bucket_key(sc): 100}
    counts, info = h.run_em(exact, contested, {X, Y}, alpha=0.5)
    assert 1099.0 < counts[X] <= 1100.0, counts[X]   # ~all 100 go to abundant X
    assert 0 <= counts[Y] < 1.0, counts[Y]
    assert abs(counts[X] + counts[Y] - 1100) < 1e-6  # mass conserved
    assert info["contested_reads"] == 100


def test_windowed_anchor_ignores_spurious_flank():
    h = _import()
    REAL, DECOY = "AAAACCCC", "GGGGTTTT"
    lib = {REAL, DECOY, "TTTTAAAA", "CCCCGGGG"}
    LF, RF, bc = "CACACACA", "GAGAGAGA", 8
    read = LF + DECOY + LF + REAL + RF      # spurious LF at 0; true barcode at offset 24
    assert read.find(LF) == 0

    cfg = h._Cfg()
    cfg.lib_set, cfg.bc_len, cfg.min_q, cfg.max_low = lib, bc, 20, 2
    cfg.paired = cfg.need_swap = cfg.single_reverse = False
    cfg.fL, cfg.fLn, cfg.fR, cfg.fRn = LF, len(LF), RF, len(RF)
    cfg.f_off, cfg.f_lo, cfg.f_hi = [24], 21, 27
    cfg.rL = cfg.rR = None
    cfg.rLn = cfg.rRn = 0
    cfg.r_off, cfg.r_lo, cfg.r_hi = [], 0, 0
    h._init_worker(cfg)

    canon, kind = h._classify(read, None, forward=True)
    assert canon == REAL and kind == "exact", f"got {canon!r}/{kind}"
    assert all(b != DECOY for b, _ in h._candidates(read, None, forward=True))


# --------------------------------------------------------------------------- #
def main():
    tests = [v for k, v in sorted(globals().items())
             if k.startswith("test_") and callable(v)]
    failed = 0
    for t in tests:
        try:
            t()
            print(f"PASS  {t.__name__}")
        except AssertionError as e:
            failed += 1
            print(f"FAIL  {t.__name__}: {e}")
        except Exception as e:  # noqa: BLE001
            failed += 1
            print(f"ERROR {t.__name__}: {type(e).__name__}: {e}")
    print(f"\n{len(tests) - failed}/{len(tests)} passed")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
