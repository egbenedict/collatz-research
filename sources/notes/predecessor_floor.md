# Predecessor-basin floor d_z: z-dependence in Shaik (2610.00523), and an exponential-decay version

Context: ideas/log.md I-008 next step 1 (floor D(Z) ≥ Z^(−B)). Line numbers refer to
sources/pdfs/shaik2026/shaik2026.txt. Tags: **[RS]** = read in Shaik; **[MD]** = my derivation
(not in the paper, not checked line by line against every hypothesis of Sections 7–10).

## 0. Inputs (RS)

- Reference law (T.law, l.610): C_n = Σ_{i≤n} 3^(i−1)·2^(−S_i) mod 3^n, K_i iid Geom(1/2).
  This is Tao's Syrac(Z/3^n) in law (reverse the iid word), so I-007's numerics apply to F_n.
- Primitive decay (E.primitive, l.7103; M.primitive, l.1943): |F_n(ξ)| ≤ C·n^(−6409) for 3∤ξ.
- Mixing envelope (Thm 6.11, l.2300): Osc(m,n) := ‖μ_n − lift(μ_m)‖₁ ≤ ε_m ≤ C6·m^(−A),
  A = 2314/25 ≈ 92.56, C6 = 2·C·20^6409 + 2.
- Radial counting (D.radial, proof ends l.3845; I.conversion l.390): d ≥ D(u) = (32/225)(e^(2u) − 1) ≥ 64u/225.
- Moment ceiling (F.fracsimple, l.3417): ⟨H_m^(3/2)⟩ ≤ B_m, B_m² ≤ 5·(57/50)^m.
  Sharper certified per-level rates: 2.1136^(1/6) ≈ 1.1328 (12-step), 1.1065 (16-step) (l.170-176).

## 1. How W_z and the density depend on z

**W_z ≍ 1/z (RS + one line of arithmetic).** Thm 10.10 (l.3978-4002):
W_y = min{Wmax, K_q·P/M_max}, M_max = X_(t0+2·3^q−1)(y), X_t(y) = (2^(2t+e(y))·y − 1)/3
(W.constructor, l.572), e(y) ∈ {1,2}, K_q = 2·3^(q−1), q = q_N ≥ 223 fixed (R.pool, l.2850),
P = Π − V(N) > 2^(−1575N−6305) (R.denominator, l.2844), Wmax = 27/227.
t0 is only needed so that X_t0(y) ≥ 16^b, b = 225, so t0 = 0 once y ≳ 16^225. Then
  W_y > 3·K_q·P·2^(−4·3^q)/y   (two periods),
  W_y^one > (3/2)·K_q·P·2^(−2·3^q)/y   (Cor 10.11, l.4004; one period).
Every root is ≥ y/3, so W_y ≤ 3·K_q/y as well. So W_z = A_0/z up to z-independent factors, with
A_0 = 2^(−Θ(3^q)): astronomically small, but a constant. Nothing else in the target step depends on z.

**Density as a function of W (RS).** Cor R.2 = Cor 10.14 (l.421-432, l.4103-4135):
  m° = max{2, ⌈(5·C6/(2W))^(25/2314)⌉},   d ≥ c_FM(W) > (4/1125)·W³·(57/50)^(−m°).
Proof (l.4119): pick m° with κ·ε_m° ≤ 2W/5 (κ = 16 ln2/9), so p_m = W − κε_m > W/2. Then
u = p³/(L²B_m²) > W³/(80·(57/50)^m) by (I.generalp) at s = 3/2. Elementary version (Prop 10.5,
l.3850): d ≥ 4W/(25·2^m), with κε_m < W/2.

**Where the loss comes from (RS):** there are two competing scalings.
1. The coarse level m must make the mixing error smaller than the score. With polynomial mixing,
   ε_m ~ m^(−92.56), this forces m ≍ W^(−25/2314) = W^(−0.0108).
2. Every level costs a fixed factor (57/50, or 2 in the elementary bound), because B_m grows
   exponentially in m (the multifractal growth of the reference fan H_m).
Together, d_z ≈ exp(−c·z^(0.0108)), a stretched exponential. Shaik states this as intrinsic to the
recipe (Prop 11.3, Q.logscale, l.4383-4402):
  a/95 − O(1) ≤ log₂log₂(1/c_FM(W)) ≤ 25a/2314 + O(1),   a = log₂(1/W),
using the envelope floor ε_n ≥ n^(−95)/300 (Q.error, l.4154). So for this recipe,
1/d_z = 2^(z^θ) with θ ∈ [1/95, 25/2314] = [0.01053, 0.01080].

**Important (RS):** the polynomial rate is not only the n^(−6409) input. The local estimate
(M.high, l.2182) is Osc ≤ 2D·n^(−a(v)) + 2·Pr(F_{n,v}). The typicality windows have width
v(√(rL) + L) with L = log n (l.2034), so Pr(F_{n,v}) ≤ (16/(q²L))·n^(1−q) + n^(−q) + n^(−q−1) with
q = v·ln2 (M.failure, l.2231), which is itself polynomial. If you substitute exponential Fourier
decay into Section 6 as written, the Fourier term becomes super-polynomially small but this
failure term stays: ε_m remains polynomial and the floor stays stretched-exponential. The
typicality event has to be redesigned too (§2).

## 2. Assume |F_n(ξ)| ≤ C·ρ^n for all n and all 3∤ξ, with ρ = e^(−I) (I-007: I ≈ 0.05498)

**Local mixing, one level at a time [MD].** Follow Lemmas 6.2–6.5 with n = m+1, using Shaik's
first-crossing head/tail split (l.2026-2198). Replace typicality by a linear-deviation lower
envelope with slope α ∈ (λ, 2), λ = log₂3, and offset a = βn (in bits):
  (T) S_i ≥ α·i − a for all i ≤ n, and K_i ≤ a for all i ≤ n.
- Failure: the walk S_i − αi has drift 2 − α > 0. The Lundberg martingale e^(−θ(S_i−αi)) gives
  Pr(inf_i (S_i − αi) < −a) ≤ e^(−θa), where θ = θ(α) > 0 solves e^(θα) = 2e^θ − 1.
  Also Pr(max K_i > a) ≤ n·2^(−a).
- Injectivity (replaces Lemma 6.3): B(w) = Σ 3^(i−1)·2^(l−S_i) ≤ 2^(l+a)/(2^α − 3). Take the gate
  Q = λn − 2a − log₂(2/(2^α−3)). Then l ≤ Q + a gives B(w) < 3^n, so Lemma 6.4 holds:
  Σ p_{k,l}² ≤ 2^(−l).
- Head length: S_k ≤ Q together with (T) gives k ≤ (Q + a)/α ≤ (λn − a)/α. The tail conductor is
  therefore r = n − k − 1 ≥ (1 − λ/α)·n + a/α − 1. (Shaik guarantees only r ≥ n/20, M.tailconductor l.2173.)
- Collision bound (M.collision): per slice, O ≤ C·ρ^r·√(3^n·2^(−l)) ≤ C·ρ^r·2^a·√(2/(2^α−3)).
  There are at most n² slices.
- Total: Osc(n−1, n) ≤ poly(n)·[e^(−n((1−λ/α)I + β(I/α − ln2))) + e^(−θβn) + 2^(−βn)].
  Balancing gives β* = (1−λ/α)·I/(θ + ln2 − I/α) and local rate a_loc = θ·β*.
- Chaining (M.triangle): ε_m ≤ Σ_{j≥m} Osc(j, j+1) ≤ C_mix·e^(−a_mix·m) for any a_mix < a_loc.

Numerically (scratch script, optimizing over α): the optimum is α* ≈ 1.79, θ ≈ 0.264, and
  a_mix ≈ 0.0316·I.  At I = 0.05498 this is a_mix ≈ 1.80e−3.
This method can never do better than a_mix ≤ (1 − λ/2)·I = 0.2075·I ≈ 0.0114, the conductor cap
with zero deviation.

**Coarse level and floor [MD], on top of the RS formulas of §1.**
m° = ⌈ln(2κ·C_mix/W)/a_mix⌉ gives κε_m° ≤ W/2. Then Cor 10.14's chain gives
  d ≥ (4/1125)·W³·(57/50)^(−m°) ≥ c·W^(3 + ln(57/50)/a_mix).
With W_z ≥ A_0/z:
  **d_z ≥ c_0 · z^(−B),   B = 3 + ln(57/50)/a_mix ≈ 3 + 4.15/I,**
  c_0 ≈ (4/1125)·(57/50)^(−1)·A_0^B·(2κ·C_mix)^(−(B−3)).
- At I = 0.05498: **B ≈ 76.** Using the 16-step rate 1.1065 instead of 1.14 (c_FM,16; the
  exponent of W is still 3 if that branch has the same p³ shape, not checked) gives B ≈ 59.
- The elementary branch (4W/(25·2^m)) gives B = 1 + ln2/a_mix ≈ 387, which is worse.
- Ideal limit of this route (a_mix = 0.2075·I, ρ = e^(−I)): B ≈ 3 + 0.131/(0.2075·I) ≈ 14.5.
- General s-moment [MD, from (I.generalp)]: if ⟨H_m^s⟩ ≲ β_s^m, then
  B(s) = s/(s−1) + ln β_s/((s−1)·a_mix). Only s = 3/2, s = 2 and s = ∞ are certified in the paper.
- Sanity: I = 0.1 gives B ≈ 42, and I = 0.2 gives B ≈ 21. B is always ≥ 3 here, far above the
  heuristic B ≈ 1 of I-008. That gap is the cost of the exponential growth of B_m.

Caveats [MD]:
- The hypothesis must hold for every level n ≥ 1 and every unit frequency. I-007 shows only the
  asymptotic rate; a constant C must absorb small n. Proving the hypothesis at all is open: it is
  Tao's question.
- I did not re-derive Sections 7–8. With better mixing the startup q_N, P and K_q would change,
  but they are z-independent, so they affect only c_0 and A_0, not B.

## 3. Other z-dependence

- W_z: linear (§1). The cycle removal (Lemma 1.5, two-period or one-period block) costs a factor
  ≤ 2 and is uniform in z. t0 matters only for z ≲ 16^225, where W saturates.
- The root pool, startup N and q_N, the moment tables, κ, L and ε_m are all fixed before the
  target. None depends on z (Thm 10.8, l.3953-3965: "none depends on which fixed target is eventually
  visited").
- The only super-polynomial mechanism is §1 (polynomial ε_m against exponential per-level
  moment cost). I found no other.
- Outside the d_z floor, the constants are horrible but z-free: A_0 = 2^(−Θ(3^q)) with q ≥ 223, so
  c_0 is like 2^(−B·Θ(3^223)). This is harmless for the I-008 inequality f(y) < D(1.5y+1) as an
  asymptotic statement, but useless numerically.
- The cutoff X past which #(G_c(z) ∩ [1,X]) ≥ d·X holds is ineffective in z (via graft
  tolerances, Lemma 10.3). This is irrelevant for lower density, but matters for any finite-X use.

## Bottom line
- Shaik as written: d_z ≥ exp(−c·z^(0.0108)), and the exponent is intrinsic to his envelope.
- With exponential primitive decay ρ = e^(−I), plus a redesigned linear-deviation typicality:
  d_z ≥ c_0·z^(−B), B ≈ 3 + ln(57/50)/(0.0316·I), which is ≈ 76 at I = 0.055.
- Confidence: §1 high. §2 medium: the argument is short and follows Shaik's lemmas, but it is
  unrefereed. The constants 0.0316 and 76 are crude and improvable.
