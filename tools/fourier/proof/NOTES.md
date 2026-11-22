# Track B: toward max_k |Φ_n(k)| ≤ Cρ^n (Tao's question). Status notes

Everything here is in `tools/fourier/proof/`. Labels: **[P]** means proved here (the proof is written
out below), **[P*]** means proved modulo the detail named, **[C]** means conjecture, **[N]** means
numerics only.

## 0. Summary

- **The exp(−cn) bound is not proved.** Below is a rigorous reduction of the oscillatory problem to a
  **positive** problem about a random walk in a fixed arithmetic picture (the gray/white plane of
  2^j mod 3^m). It comes with an exact statement of what is missing, hypothesis (★) in §3, and an
  explanation of why the missing piece is exactly where Tao's argument stops at polynomial decay
  (the "chain" scenario, §5).
- Proved:
  - pair factorization, Lemma 1 [P];
  - martingale reduction, Lemma 2 [P]. Theorem A: (★) ⇒ max|Φ_n| ≤ μ^{εn/2} + Ce^{−cn}, with
    μ = 1 − (1 − cos(π/5))/4 = 0.95225 [P];
  - triangle geometry at threshold 1/5, Lemmas 3–4 and 6 [P], Lemma 5 [P*]. Lemmas 3a, 3b, 3c and
    4 were also brute-force checked with exact integers (m ≤ 16, |j| ≤ 300).
- Computed: the positive majorant F_n (Lemma 1) itself decays exponentially, at rate ≈ 0.034–0.037
  per level for n ≤ 600, against I₀ = 0.055 for Φ. So passing to absolute values costs only about
  1/3 of the exponent. For n ≤ 600 and the low frequencies ξ = 2^k, 0 ≤ k ≤ 1.6n, this gives a float
  (not interval-certified) upper bound, e.g. |Φ_600(k)| ≤ 8.73·10⁻¹¹. This is finite-n evidence
  only.
- Lemmas that fail:
  - (F1) uniform L-block sup-norm contraction, for every L;
  - (F2) a uniform post-exit block lemma H(L,q), for every fixed L: it fails globally, although it
    holds in every window we computed.

## 1. Setup [P]

- Write θ(m,j) = (2^j mod 3^m)/3^m ∈ R/Z for j ∈ Z (negative j means inverses), and ‖·‖ for the
  distance to Z.
- Two exact relations: **θ(m,j+1) = 2θ(m,j)** and **θ(m−1,j) = 3θ(m,j)** (mod 1).
- Unrolling the recursion gives
  Φ_n(k) = E Π_{m=1}^{n} e(θ(m, j_m)), with j_{n+1} = k and j_m = j_{m+1} − a_m, where the a_m are
  iid Geom(1/2) on {1,2,…}.
- The "phase point" at level m is (m, j_m). The walk moves left (m decreases) and down
  (j decreases by a_m).

## 2. Lemma 1: pair factorization [P]

Group levels (m, m−1) and set b = a_m + a_{m−1}. Given b, a_m is uniform on {1,…,b−1}. The phase at
level m−1 and everything below depend only on j_{m+1} − b. So, exactly,

  Φ_m(J) = Σ_{b≥2} 2^{−b} [Σ_{a=1}^{b−1} e(θ(m,J−a))] · e(θ(m−1,J−b)) · Φ_{m−2}(J−b).

The inner sum equals Σ_{r=0}^{b−2} e(2^r x) with x = θ(m, J−b+1). Hence

  |Φ_m(J)| ≤ F_m(J) := Σ_b 2^{−b}(b−1) g_b(θ(m,J−b+1)) F_{m−2}(J−b),

where g_b(x) = |(b−1)^{−1} Σ_{r<b−1} e(2^r x)|, F_0 = 1 and F_1 = |Φ_1| = 1/√3. Note g_3(x) = |cos πx|.

- Pairing can be chosen adaptively (by a stopping rule on the past); the factorization survives.
- With fixed pairing, F is E Π_p g_{b_p}(·) over the pair walk J_{p+1} = J_p − b_p, where
  P(b) = (b−1)2^{−b}.
- Checked: |Φ| ≤ F for every k, n ≤ 15 (`majorant.py`).

## 3. Lemma 2 and Theorem A: martingale reduction [P]

- **Q-walk.** Q_p = (n−2p, J_p − 2), p < P = ⌊n/2⌋, with steps (−2, −b_p). For b_p = 3 the point
  in the g-factor is exactly Q_p.
- **Colours.** Call a point *white* if ‖θ‖ ≥ 1/5 and *gray* otherwise. Let w_p = 1[Q_p white]. It is
  predictable, since it depends on b_0…b_{p−1}. Let V = Σ_p w_p.
- **The bound.** Put λ = cos(π/5) and N = Σ_p 1[b_p=3]·w_p. Then F ≤ E λ^N, because g ≤ 1 always and
  g_3 ≤ λ at white points.
- **The martingale.** P(b=3) = 1/4, so E[λ^{1[b_p=3]w_p} μ^{−w_p} | F_p] = 1 for μ = 1 − (1−λ)/4.
  So λ^N μ^{−V} is a martingale with mean 1, and for every t:

  **|Φ_n(k)| ≤ μ^t + P_k(V < t),   μ = 0.95225.**

**Hypothesis (★)** [C]: there are ε, c, C > 0 such that, uniformly in n and k (k mod 2·3^{n−1}),
P_k(V < εP) ≤ C e^{−cn}.

**Theorem A** [P]: (★) implies max_k |Φ_n(k)| ≤ μ^{εn/2 − 1} + C e^{−cn}. For ε = 0.05 this gives
ρ = μ^{0.025} = 0.99878 per level, provided c ≥ 0.0013.

(★) concerns only the deterministic gray set and a random walk. No oscillation is left in it.
Necessarily c ≤ I₀: the walk that stays in Tao's main triangle (rows 0 ≤ j < m·log₂3) sees almost no
white points.

## 4. Triangle geometry at threshold 1/5

Signed values s ∈ (−½, ½].

- **Lemma 3** [P].
  - **Vertex.** Define a vertex v = (m_v, j_v), with value δ = s(v), as a gray point where
    s(m_v, j_v−1) ≠ δ/2 and s(m_v+1, j_v) ≠ δ/3. Then T(v) = {(m_v−a, j_v+b) : a,b ≥ 0,
    3^a 2^b |δ| < 1/5} consists of gray points with value exactly 3^a2^bδ, since there is no
    wrap-around below ½.
  - **(3a)** Every gray point lies in such a T(v). Proof: by the doubling and tripling relations, a
    gray value x has lower neighbour x/2 or x/2 + ½, and right neighbour x/3 or (x±1)/3. The values
    other than x/2 and x/3 have norm ≥ 1/4 > 1/5. So the parent chain is unique, and values shrink
    by a factor ≥ 2 at each step, bounded below by 3^{−level}.
  - **(3b)** Below the bottom edge, at columns with 3^a|δ| < 1/5:
    - s(m_v−a, j_v−1) = ½ + 3^aδ/2 ∈ (0.4, 0.6);
    - s(m_v−a, j_v−2) ∈ ¼ + 3^aδ/4 + {0, ½}, which is white.
    - Proof: at the vertex the lower value is δ/2 + ½. Multiply by 3^a (odd).
  - **(3c)** A Q-step (−2, −b) from a point of T(v) that stays in the cone rows (b ≤ height) but
    leaves T(v) lands on value 9·2^{−b}·x with |x| < 1/5. That has norm in [1/5, 9/20], so it is
    white.
  - Brute-force check in `lemma_check.py`: 1124 vertices, 518 hypotenuse exits, 0 failures.
- **Lemma 4: separation** [P].
  - **Statement.** If v₁ ≠ v₂ are vertices with m₁ ≤ m₂ and j₁ ≤ j₂, then
    2^{j₂−j₁}|δ₁| ≥ ½ or 3^{m₂−m₁}|δ₂| ≥ ½.
  - **Proof.** Otherwise the point w = (m₁, j₂) has value 2^{j₂−j₁}δ₁ = 3^{m₂−m₁}δ₂ (as reals). Its
    lower neighbour is 2^{j₂−j₁−1}δ₁ from v₁'s column, and ½ + 3^{m₂−m₁}δ₂/2 from v₂'s vertex row.
    These differ by ½, a contradiction. The case j₁ = j₂ contradicts v₁ being a vertex.
  - Checked: 346,283 ordered pairs, 0 failures.
- **Lemma 5: exit** [P*].
  - **Statement.** Suppose the Q-walk is in T(v) at height h ≥ 0 above the bottom row, and leaves
    T(v) at this step. Hypotenuse exits are white with probability 1 (3c). For a bottom exit, if the
    landing column is still in T's bottom range, P(land on row j_v−1 or j_v−2 | exit) ≥ 1/4.
    Proof: P(b=c | b≥c) = (c−1)/(2c), and for h = 0, P(b=2) = 1/4. Those rows are white (3b).
  - **Gap:** the "tip" case. The landing column is past the left tip of T. This forces h ≤ 3 within
    2 columns of the tip. The landing colour is then not controlled. Not resolved.
- **Lemma 6: capacity** [P]. Any triangle restricted to levels ≤ m has height ≤ m·log₂3 − log₂5,
  because |s| ≥ 3^{−m}. One triangle can hold the walk from level m to the end only on a
  large-deviation event: Σb ≤ 1.585m over (1−ε)m/2 pairs, whose mean is 2(1−ε)m.

## 5. The obstacle, made precise [P for the scenario, C for its absence]

- Combine Lemmas 5 and 6.
  - **What one triangle takes.** A walk at level m that falls into a triangle of relative height
    γ·m stays about γm/2 levels and exits at level (1−γ/2)m. That is typical behaviour, not a large
    deviation. Example: the main triangle, γ = log₂3, exit level ≈ 0.21n.
  - **The chain.** Suppose the arithmetic puts, right under each exit, another triangle of relative
    height ≥ γ₀. Then the gray stretches have lengths that form a geometric series summing to all
    n levels. There are only O(log n / log(1/(1−γ₀/2))) exits, so V = O(log n).
- Lemma 2 then gives only a polynomial bound n^{−c}. **This is the polynomial barrier** (compare
  Tao's n^{−A}).
- Exponential decay needs, in the large-deviation sense (★), that triangles of height ≥ H take up a
  fraction ≤ e^{−cH} of the walk's time.
- **Diophantine form of a two-link chain.**
  - Triangles T₁ and T₂ with vertex residues r₁, r₂ (2^{j_i} ≡ r_i mod 3^{m_i}), with T₂ under T₁'s
    exit (m₂ ≈ 0.21m₁, d = j₁ − j₂ in a window of width O(√m)), require 2^d r₂ ≡ r₁ (mod 3^{m₂}).
  - For bounded r_i, Yu's 3-adic Baker bound gives v₃(2^d r₂ − r₁) ≪ log d, which excludes it.
  - For near-maximal triangles (|r_i| ≤ 2^{ηm}), the Baker bound is about (ηm)²·log d. That is
    useless. Heuristically such chains are absent once 2^{3ηn} < 3^{0.21n}, but nothing here proves
    that.

## 6. Numerics [N]

- **`majorant.py` (exact, all k, n ≤ 15).**
  - max F_n is 0.649, 0.522, 0.439, 0.384, 0.339, 0.301, 0.269 for even n = 2…14, against
    max|Φ| ≈ 0.019 at n = 14.
- **`band_majorant.py N J0 floor`.**
  - Band J ∈ [−J0, 1.6N]. Floor 1 below the band gives a rigorous majorant for J in the band.
  - For k ≥ 0, floor 1 and floor 0 give identical values to all printed digits, for J0 = 300 and
    for J0 = 600.
  - max_{k≥0} F_n: 8.16e-3 (n=100), 2.23e-4 (200), 6.01e-6 (300), 1.39e-7 (400), 3.35e-9 (500),
    8.72e-11 (600). This is at k ≈ n·log₂3 − 4.
  - Ratio per 50 levels: 0.153 → 0.162, rising slowly. That fits e^{−In}·n^{−c} with I ≈ 0.034 and
    c ≈ 1.3.
  - **Float caveat.** Accumulated rounding is ≲ 10⁻¹⁴ per level, since the recursion has total
    weight ≤ 1. Not interval-certified.
- **`qwalk_mc.py` (exact big-int phases).**

  | start | V/P | 1% quantile |
  |---|---|---|
  | random k | 0.60 (gray density 2/5) | ≥ 0.44 (n = 100–400) |
  | top k | 0.13–0.145 | 0.02–0.07 |

  For top k, white visits come only after the exit from the main triangle at level ≈ 0.21n.
- **`block_gray.py 300 400 L`.**
  - Worst-case P(the next L Q-positions are all gray) over 2-adic starts (m ≤ 300, J ∈ [−J0+8L, 0)):
    1.0 (L=2), 0.74 (L=4), 0.080 (L=8), 1.1e-4 (L=16), 6.6e-8 (L=24).
  - Over main-triangle starts it is 1.0 for every L.
- **`triangle_census.py 600 1200`** (m ≤ 600, −1200 ≤ J < 0).
  - 96,000 vertices. #(H ≥ h) matches the generic law 96000·2^{−h} for every h from 0 to 15.
  - The largest height is 15.7, about log₂(#vertices). Nothing anomalous below the main triangle.

## 7. Lemmas that fail

- **(F1)** sup_J of the L-step majorant operator applied to 1 is 1 − o(1) for every L. Witness: J
  deep in the main triangle, where all phases are within 2^{−depth} of 1. So no uniform
  sup-norm block contraction exists. This is the stated obstacle, confirmed.
- **(F2)** H(L,q): "after the first white visit, every L-block contains a white Q with conditional
  probability ≥ 1−q". It would give (★), hence Theorem A, with ρ = μ_B^{0.105/L} for the top start,
  where μ_B = 1 − (1−q)(1−λ)/4.
  - With the window values (L = 16, q = 1.2e-4) this is ρ ≈ 0.99968.
  - But H(L,q) is **false globally** for each fixed L. For large m there are r-triangles (r small) of
    height ≈ m·log₂3 − log₂r anywhere in the period. Points deep inside them can be reached from a
    white point just above them.
  - So only an averaged statement like (★) can be true.

## 8. What would close the gap

Prove (★). Equivalently, by Lemmas 5–6, exclude chains of triangles of linear relative height under
the walk with exponential probability, i.e. show that the walk's gray time is mostly spent in small
triangles. The census in §6 says this is true in every window computed. A proof needs
equidistribution of 2^j mod 3^m over linear-length windows of j, or a Baker-type input for
near-maximal residues. Neither is available.

## Files

- `majorant.py`: exact |Φ| against the majorants F and U, n ≤ 15.
- `band_majorant.py`: banded F, with a rigorous floor (1) or the heuristic floor (0).
- `lemma_check.py`: exact brute-force checks of Lemmas 3–4.
- `qwalk_mc.py`: Monte Carlo of the Q-walk.
- `block_gray.py`: worst-case block gray persistence.
- `triangle_census.py`: census of triangle heights.
