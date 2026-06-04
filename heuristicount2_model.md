# heuristicount2: statistical model for quality-aware barcode counting

This document records the model and the estimator that `heuristicount2.py` uses, with the
derivations and the places where the naive approach is biased. It is self-contained. The
framework is standard finite-mixture / EM estimation (the same machinery RSEM, kallisto,
and salmon use to apportion multi-mapping reads, and that misclassification correction
uses in epidemiology). The specialization to a one-base-apart guide pair, and the closed
forms below, are derived here and were checked both by Monte Carlo and by independent
re-derivation.

## 1. Why this is needed

A pooled library can contain guides that differ by a single base on purpose (a perfect
guide and an intentional single-mismatch variant). Designing them one base apart removes
the safety margin that protects ordinary libraries: in a library where every pair differs
by several bases, no single sequencing error can turn one guide into another, so exact
matching is safe. When two guides differ at one position, a single error at that position
bridges them, and the rate of cross-contamination becomes the per-base error rate rather
than the (vanishing) probability of several coordinated errors.

## 2. Model

Two library sequences X and Y are identical except at one position p, where X has base
alpha and Y has base beta (alpha != beta). True abundances are N_X, N_Y; the mixing
fraction is pi = N_X / (N_X + N_Y). Each molecule is read once; at position p we observe
a base call c with PHRED quality Q, giving error probability

    e = 10^(-Q/10).

Error model (symmetric, uniform miscall): the true base is read correctly with
probability (1 - e), and on error becomes each of the other three bases with probability
e/3. So for a guide g with base g_p at position p,

    P(c | g) = (1 - e)   if c == g_p
             = e/3       otherwise.

PHRED enters the estimator only through e. If a platform's Q is miscalibrated or its
substitutions are not uniform, replace P(c | g) with an empirical substitution matrix
P(c | true, Q); nothing structural below changes, only the numbers.

## 3. The naive hard estimator is biased

Hard counting calls c == alpha as X, c == beta as Y, and discards the other two bases.
By linearity over the N_X molecules of X and the N_Y of Y,

    E[nX_hat] = N_X (1 - e) + N_Y (e/3)
    E[nY_hat] = N_X (e/3)   + N_Y (1 - e)

which is the confusion matrix

    [ E nX_hat ]   [ 1-e    e/3 ] [ N_X ]
    [          ] = [            ] [     ].
    [ E nY_hat ]   [ e/3    1-e ] [ N_Y ]

The observed ratio (with r = N_X / N_Y) is shrunk toward 1:

    r_hat = [ r (1 - e) + e/3 ] / [ (1 - e) + r (e/3) ].

If Y is truly absent (N_Y = 0) you still expect nY_hat = N_X (e/3) phantom reads, so the
ratio of expectations is capped:

    r_hat_max = 3 (1 - e) / e     (about 297 at Q20, about 2997 at Q30).

This is a contamination floor on the depleted guide. Note it bounds the ratio of
expectations, not the realized ratio of a finite sample (which is high variance and can
exceed it by chance when the absent guide's count is tiny). With per-read e_i, the floor
generalizes to 3 (1 - e_bar) / e_bar with e_bar the mean error.

When N_X = N_Y the off-diagonal leakages X -> beta and Y -> alpha occur at the same rate
e/3 and cancel, so the ratio is unbiased and PHRED is irrelevant. This is exactly the
balanced corner. The relative bias on the minor guide grows as (e/3) * r, so it is
negligible for neighbors of similar abundance and dominant only when a low-abundance guide
sits one base from a high-abundance neighbor.

## 4. The ideal estimator: per-read mixture maximum likelihood

The data are a two-component mixture. For reads i = 1..n with calls c_i and error
probabilities e_i,

    L(pi) = prod_i [ pi * P(c_i | X, e_i) + (1 - pi) * P(c_i | Y, e_i) ].

The maximum-likelihood pi is the ideal point estimate: it is consistent and asymptotically
efficient (it attains the Cramer-Rao bound). It is found by EM.

E-step (responsibility, the posterior that read i is X):

    w_i = pi * Lx_i / ( pi * Lx_i + (1 - pi) * Ly_i ),

    where Lx_i = (1 - e_i) if c_i == alpha else e_i/3,
          Ly_i = (1 - e_i) if c_i == beta  else e_i/3.

M-step:

    pi <- (1/n) sum_i w_i.

At convergence N_X_hat = pi_hat * M and N_Y_hat = (1 - pi_hat) * M, where M is the total
read depth (only the fraction pi is identifiable from the calls; the absolute counts need
the depth, which the pipeline supplies as the number of reads).

The responsibility for the three observed cases makes the role of PHRED explicit:

    c == alpha:  w_i = pi (1-e_i) / [ pi (1-e_i) + (1-pi) e_i/3 ]   (-> 1 as e_i -> 0)
    c == beta:   w_i = pi e_i/3   / [ pi e_i/3   + (1-pi)(1-e_i) ]   (-> 0 as e_i -> 0)
    c == other:  w_i = pi                                            (no information)

A high-quality call moves the estimate hard toward one guide; a low-quality call barely
moves it; an off-base call contributes nothing beyond the prior.

### Why the soft estimator is optimal

N_X is the latent count sum_i 1{Z_i = X}. For any per-read contribution a_i that is a
function of the data, the mean squared error E[(sum a_i - sum 1{Z_i = X})^2] is minimized
term by term by the conditional expectation a_i = E[1{Z_i = X} | data_i] = w_i (the
defining property of conditional expectation; Rao-Blackwell). So sum_i w_i is the unique
minimum-MSE estimator, and it is unbiased: E[w_i] = pi. Any hard rule (a_i in {0, 1}, or
discard) is a coarser quantization of w_i and has weakly larger MSE, strictly larger unless
every base is error-free or fully random. This holds with pi known; with pi estimated by
EM, sum_i w_i is a plug-in that is asymptotically MMSE and unbiased.

### Multiple positions and many guides

For a read with per-base error e_j and observed barcode b, the likelihood of guide g
factorizes over positions. The factor prod_j (1 - e_j) is common to all guides and cancels
in the responsibility, so define a score over mismatched positions only:

    s_ig = prod_{ j : b_j != g_j } (e_j / 3) / (1 - e_j),

    w_ig = pi_g * s_ig / sum_h pi_h * s_ih.

A mismatch at a high-quality base multiplies by about e_j/3 (heavy penalty); at a
low-quality base by a factor near 1 (light penalty). The factor (e_j/3)/(1 - e_j) is a
penalty (below 1) for e_j < 3/4, which holds for any realistic quality (Q >= 2). Only
discriminating positions (where candidate guides differ) carry information; positions
where all candidates agree cancel.

## 5. Why a fixed PHRED cutoff is only an approximation

Hard-calling an observed beta as Y is justified only when the posterior P(X | c = beta) is
small, that is pi_X (e/3) << pi_Y (1 - e), i.e. e <~ 3/r, i.e.

    Q >~ 10 * log10(r / 3).

The exact boundary is e = 3/(r + 3); the form above is its small-e limit. The threshold is
positive only for r > 3 (X dominant); for r <= 3 a beta call can be trusted at any quality.
The required quality therefore scales with the abundance ratio, so a single fixed cutoff is
not generally correct. The soft estimator needs no cutoff: it puts the continuous e_i into
the posterior and becomes skeptical of a low-quality minor-guide call in exact proportion
to how outnumbered that guide is.

## 6. What heuristicount2 does

Default mode (hard, quality-gated): a base at Q >= --min-quality must match a guide
exactly; only sub-threshold (or N) bases are treated as wildcards, capped at
--mismatches positions. A read whose low-quality base could equally have come from two
guides is reported as ambiguous and not counted, rather than guessed. This is the
threshold approximation of Section 5 and is correct when neighbors are similar in
abundance.

Soft mode (--soft): the EM estimator of Section 4. Per read, the candidate set is the
library guides reachable by varying the low-quality (and N) positions of the extracted
barcode, each with score

    s_g = prod_{ j in low-quality positions } ( (1 - e_j) if b_j == g_j else e_j/3 ),

which is the Section 4 score restricted to the positions that can vary (high-quality
positions are fixed and cancel). Reads with a single candidate are counted whole; reads
with two or more are accumulated as contested observations. EM then apportions the
contested reads in proportion to the abundances, which are pinned by the whole-count reads.
Output is the expected (generally fractional) count per guide.

### Practical notes from the derivation

1. Identifiability: only the fraction pi is identifiable from the discriminating base.
   Absolute counts come from the total read depth, which the pipeline has.
2. Other-base calls are not discarded by the soft estimator. They are uninformative
   (equal likelihood under all candidates) and cancel, but a correct estimator keeps them.
   Hard counting drops them, a small extra inefficiency.
3. Boundary behavior: pi in {0, 1} are spurious EM fixed points (w == 0 at pi = 0, w == 1
   at pi = 1 for any data). For a guide with no whole-count reads that appears only in
   contested observations, EM must be initialized in the interior or it sticks at zero.
   heuristicount2 uses a small Dirichlet pseudocount (--em-prior, default 0.5) on the
   abundances inside the E-step, which keeps the estimate interior and is the MAP estimator
   under a symmetric Dirichlet prior. This also means a dropped-out guide is recovered to
   near zero asymptotically and is biased slightly upward by the pseudocount in finite
   samples, which is the safe direction (it does not invent a hard floor the way naive
   counting does).
4. Calibration: everything rests on PHRED being P(miscall) with uniform substitution. The
   estimates are only as good as that calibration; swap in an empirical substitution matrix
   if your platform warrants it.
