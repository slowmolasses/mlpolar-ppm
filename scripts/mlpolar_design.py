"""
mlpolar_design.py
=================
Multi-Level Polar Code Design Tool for L-ary PPM + Photon-Detection Channels.

Theory background
-----------------
We consider a classical channel induced by L-pulse-position modulation (PPM)
with a photon-counting receiver.  The PPM alphabet has size L = 2^K, and the
channel transition matrix is:

    W(y | x)  =  p            if y == x          (correct detection)
              =  q            if y == L           (erasure / no photon)
              =  (1-p-q)/(L-1)  otherwise         (wrong-slot / dark count)

where x ∈ {0,...,L-1} is the transmitted PPM slot and y ∈ {0,...,L} is the
received symbol (L = erasure).

This is the discrete classical channel that arises from:
  - Coherent state L-PPM at the transmitter (block length L Hadamard-encoded BPSK)
  - A Green Machine (Hadamard stage) + photon detection at the receiver

K-fold Binary Sequential Partition (K-SBP)
--------------------------------------------
A K-SBP (also called a K-level set-partition labeling) is a sequence of nested
binary partitions of the input alphabet X = {0,...,L-1}:

    P_0 ⊃ P_1 ⊃ ... ⊃ P_K

where P_0 = {X} (trivial) and P_K = {{x} : x ∈ X} (singletons), and each
refinement P_i → P_{i+1} splits every cell of P_i into exactly two equal halves.

Each refinement step defines a **bit channel** B^(i)_λ  (Seidl et al. [1] Eq. 4):
- Input:   b_i ∈ {0, 1}   (which half within the cell of P_i)
- Output:  (Y, b_0, ..., b_{i-1})  (channel output + all previously decoded bits)

The **capacity** of B^(i)_λ is the conditional mutual information (Seidl Eq. 4):

    I(B^(i)_λ) = I(b_i ; Y | b_0, ..., b_{i-1})

These K capacities satisfy the chain rule:

    Σ_{i=0}^{K-1} I(B^(i)_λ) = I(X; Y)                              (1)

because b_0,...,b_{K-1} = λ(X) is a BIJECTION of X (any bijective labeling
preserves I(X;Y) by data-processing equality applied to λ).  This is the
fundamental reason a K-SBP decomposes channel capacity without loss.

SP (Set-Partitioning / Natural Binary) Labeling
------------------------------------------------
The SP labeling maps slot index x to its K-bit binary representation:

    λ_SP(x) = (b_0, b_1, ..., b_{K-1})   where  b_j = (x >> (K-1-j)) & 1

This means:
  - Level 0 (b_0 = MSB): splits X into {0,...,L/2-1} vs {L/2,...,L-1}
  - Level 1 (b_1):       splits each half into two quarters
  - ...
  - Level K-1 (b_{K-1} = LSB): splits each pair into singletons

SP labeling is chosen (over Gray code, etc.) because it maximises the variance
of the individual bit-channel capacities [1, Sec. IV-D], which makes the polar
codes at each level operate closer to their respective channel capacities.

Multi-Level Polar Coding
-------------------------
For each bit level i, a BINARY polar code of length N = 2^n is constructed for
the bit channel B^(i)_λ with rate R_i = I(B^(i)_λ)  (capacity rule [1] Sec.
IV-A).  The frozen set for level i is chosen by applying Arikan's channel
polarization to N copies of B^(i)_λ, selecting the N - K_info^(i) least
reliable synthetic channels as frozen (fixed to zero).

Density Evolution via Bhattacharyya Parameters
-----------------------------------------------
For a binary-input channel W with transition probabilities P(y|0), P(y|1),
the Bhattacharyya parameter is:

    Z(W) = Σ_y √(P(y|0) · P(y|1))

It satisfies Z(W) ∈ [0, 1] and  I(W) ≥ 1 - h(Z(W))  (h = binary entropy).

Under Arikan's polar transform F_2^{⊗n}, the Bhattacharyya parameters of the
two child synthetic channels satisfy the upper bounds [2]:

    Z(W^−) ≤ 2·Z(W) - Z(W)^2     (XOR / "check" channel  — worse)
    Z(W^+)  = Z(W)^2              (copy / "variable" channel — better, exact)

These bounds are EXACT for BEC channels and are valid upper bounds for all
binary-input symmetric channels.  We compute Z(B^(i)_λ) exactly from the PPM
channel and then apply these recursive bounds to obtain reliability estimates
for all N synthetic channels.  Sorting by Z ascending gives the frozen set.

References
----------
[1] M. Seidl, A. Schenk, C. Stierstorfer, J. B. Huber, "Polar-Coded
    Multi-Level Modulation," IEEE Trans. Commun., 2013.
[2] E. Arikan, "Channel Polarization: A Method for Constructing
    Capacity-Achieving Codes for Symmetric Binary-Input Memoryless Channels,"
    IEEE Trans. Inf. Theory, 2009.
"""

import numpy as np
from typing import Optional


# ---------------------------------------------------------------------------
# 1.  PPM+PD channel transition matrix
# ---------------------------------------------------------------------------

def ppm_channel(K: int, p: float, q: float) -> np.ndarray:
    """Build the transition matrix of the L-ary PPM + photon-detection channel.

    Parameters
    ----------
    K : int
        Number of bit levels;  L = 2^K  is the PPM alphabet size.
    p : float
        Probability of correct detection: W(x|x) = p.
    q : float
        Probability of erasure (photon lost): W(L|x) = q for all x.

    Returns
    -------
    T : np.ndarray, shape (L, L+1)
        T[x, y] = W(y | x).
        Rows x ∈ {0,...,L-1}: transmitted slot.
        Cols y ∈ {0,...,L-1}: received slot;  y=L: erasure.
    """
    L = 1 << K           # 2^K
    r = (1.0 - p - q) / (L - 1)   # wrong-slot probability
    if r < 0:
        raise ValueError(f"p+q={p+q} > 1 is impossible; need p+q < 1.")
    T = np.full((L, L + 1), r, dtype=np.float64)
    for x in range(L):
        T[x, x] = p       # correct detection
        T[x, L] = q       # erasure
    return T


# ---------------------------------------------------------------------------
# 2.  K-SBP under SP labeling: transition matrices of bit channels
# ---------------------------------------------------------------------------

def ksbp_transition_matrices(K: int, p: float, q: float) -> list[np.ndarray]:
    """Compute the K transition matrices of the K-SBP bit channels B^(i)_λ.

    Uses the SP (natural binary / set-partitioning) labeling, which maps slot
    index x to its K-bit binary representation, MSB first.

    The bit channel B^(i)_λ has:
      - Binary input:  b_i ∈ {0, 1}
      - Output space:  Y × {0,1}^i   (PPM output Y  and  decoded prefix)

    Its transition matrix has shape (2, (L+1)·2^i):

        M_i[b_i,  prefix * (L+1) + y]
          = P(Y=y, b_{0:i-1}=prefix | b_i)
          = (1 / 2^i) · mean_{x: prefix(x)=prefix, bit_i(x)=b_i} W(y | x)

    where "prefix(x)" = the top i bits of x under SP labeling, and
    "bit_i(x)" = the (i+1)-th most-significant bit of x.

    The factor 1/2^i comes from P(prefix | b_i) = 1/2^i, which holds because
    under a uniform input and bijective SP labeling the prefix is independent
    of b_i.

    Parameters
    ----------
    K, p, q : as in ppm_channel().

    Returns
    -------
    matrices : list of K np.ndarray, each of shape (2, (L+1)*2^i).
        matrices[i] is the stochastic transition matrix of B^(i)_λ.
        Row sums equal 1 for each row.

    Why is this a valid K-SBP?
    ~~~~~~~~~~~~~~~~~~~~~~~~~~
    The SP labeling λ_SP is a bijection X → {0,1}^K, so the K bits b_0,...,b_{K-1}
    are a sufficient statistic for X.  At each level i the split is:
      {x : b_j(x)=s_j for j<i, b_i(x)=0}  vs
      {x : b_j(x)=s_j for j<i, b_i(x)=1}
    Each cell of the partition P_i has exactly 2^(K-i) elements, split equally
    by bit i.  This produces K balanced binary partitions, each further refining
    the previous one — the defining property of a K-SBP.  The chain-rule
    identity Σ I(B^(i)_λ) = I(X;Y) holds because the bijection preserves all
    information: I(λ(X);Y) = I(X;Y).
    """
    L = 1 << K
    T = ppm_channel(K, p, q)          # shape (L, L+1)
    n_out = L + 1

    matrices = []
    for i in range(K):
        group_size = L >> i            # = 2^(K-i): total inputs per prefix group
        half       = group_size >> 1   # = 2^(K-i-1): inputs per (prefix, bit) cell
        n_prefix   = 1 << i            # = 2^i: number of distinct prefixes

        # Output width: (L+1) channel outputs × 2^i prefixes
        M = np.zeros((2, n_out * n_prefix), dtype=np.float64)

        for prefix in range(n_prefix):
            start   = prefix * group_size
            # SP labeling: within a prefix group, the first half has b_i=0
            #              and the second half has b_i=1 (MSB-first binary order)
            inputs_0 = np.arange(start,        start + half)        # b_i = 0
            inputs_1 = np.arange(start + half, start + group_size)  # b_i = 1

            # Average W(y|x) over each half — this gives P(Y=y | b_i=v, prefix=prefix)
            p0_y = T[inputs_0, :].mean(axis=0)   # shape (L+1,)
            p1_y = T[inputs_1, :].mean(axis=0)   # shape (L+1,)

            # Weight by P(prefix | b_i) = 1/n_prefix  (prefix ⊥ b_i under uniform X)
            col_off = prefix * n_out
            M[0, col_off : col_off + n_out] = p0_y / n_prefix
            M[1, col_off : col_off + n_out] = p1_y / n_prefix

        # Sanity check: each row must sum to 1
        assert np.allclose(M.sum(axis=1), 1.0, atol=1e-10), \
            f"Row sum error at level {i}: {M.sum(axis=1)}"

        matrices.append(M)

    return matrices


# ---------------------------------------------------------------------------
# 3.  Bit-channel capacities  I(B^(i)_λ)  =  I(b_i ; Y | b_0,...,b_{i-1})
# ---------------------------------------------------------------------------

def _mi_binary(p0_y: np.ndarray, p1_y: np.ndarray) -> float:
    """MI between a uniform binary input and an output with P(·|0)=p0_y, P(·|1)=p1_y.

    I(b; Y) = H(Y) - H(Y|b)
            = H(mix) - 0.5*(H(p0_y) + H(p1_y))

    where mix = 0.5*(p0_y + p1_y).
    """
    mix = 0.5 * (p0_y + p1_y)
    def H(v): return -np.sum(v * np.log2(v, where=v > 0, out=np.zeros_like(v)))
    return H(mix) - 0.5 * (H(p0_y) + H(p1_y))


def bit_channel_capacities(K: int, p: float, q: float) -> list[float]:
    """Compute the K symmetric bit-channel capacities I(B^(i)_λ) under SP labeling.

    Implements the conditional mutual information definition (Seidl Eq. 4):

        I(B^(i)_λ) = I(b_i ; Y | b_0,...,b_{i-1})
                   = E_{b_{0:i-1}} [ I(b_i ; Y | b_{0:i-1} = s) ]

    where the expectation is uniform over the 2^i possible prefixes s.

    Parameters
    ----------
    K, p, q : as in ppm_channel().

    Returns
    -------
    capacities : list of K floats.
        capacities[i] = I(B^(i)_λ) in bits.

    Note: Σ capacities = I(X;Y) by the chain rule (verified internally).
    """
    L = 1 << K
    T = ppm_channel(K, p, q)

    capacities = []
    for i in range(K):
        group_size = L >> i
        half       = group_size >> 1
        n_prefix   = 1 << i

        # Average MI over all prefixes (each prefix has equal probability 1/n_prefix)
        cap_i = 0.0
        for prefix in range(n_prefix):
            start    = prefix * group_size
            p0_y = T[start        : start + half,        :].mean(axis=0)
            p1_y = T[start + half : start + group_size,  :].mean(axis=0)
            cap_i += _mi_binary(p0_y, p1_y) / n_prefix

        capacities.append(cap_i)

    return capacities


# ---------------------------------------------------------------------------
# 4.  Bhattacharyya parameters of the bit channels
# ---------------------------------------------------------------------------

def bit_channel_bhattacharyya(K: int, p: float, q: float) -> list[float]:
    """Compute Bhattacharyya parameters Z(B^(i)_λ) for each bit level.

    Z(W) = Σ_y √(P(y|0) · P(y|1))

    For the bit channel B^(i)_λ with joint output (y, prefix), using the fact
    that prefix ⊥ b_i under uniform input:

        Z(B^(i)_λ) = E_{prefix} [ Σ_y √(P(y|b_i=0,prefix) · P(y|b_i=1,prefix)) ]

    Parameters
    ----------
    K, p, q : as in ppm_channel().

    Returns
    -------
    z_params : list of K floats in [0, 1].
    """
    L = 1 << K
    T = ppm_channel(K, p, q)

    z_params = []
    for i in range(K):
        group_size = L >> i
        half       = group_size >> 1
        n_prefix   = 1 << i

        Z_i = 0.0
        for prefix in range(n_prefix):
            start    = prefix * group_size
            p0_y = T[start        : start + half,       :].mean(axis=0)
            p1_y = T[start + half : start + group_size, :].mean(axis=0)
            Z_i += np.sum(np.sqrt(p0_y * p1_y)) / n_prefix

        z_params.append(float(Z_i))

    return z_params


# ---------------------------------------------------------------------------
# 5.  Bhattacharyya density evolution for polar codes
# ---------------------------------------------------------------------------

def _polar_bhatt_evolution(Z0: float, n: int) -> list[float]:
    """Evolve the Bhattacharyya parameter Z0 of a binary-input channel W through
    n stages of Arikan's polar transform, yielding 2^n synthetic channel parameters.

    Recursive upper bounds (exact for BEC; upper bounds for general channels):
        Z(W^−) ≤ 2·Z − Z^2       (XOR/check branch, indices 2k)
        Z(W^+)  = Z^2             (direct/copy branch, indices 2k+1)

    The resulting array Z_out[j] approximates the Bhattacharyya parameter of
    the j-th synthetic channel.  Smaller Z → more reliable → likely info channel.

    Returns list of 2^n values, ordered by their natural tree index.
    """
    channels = [Z0]
    for _ in range(n):
        new = []
        for Z in channels:
            new.append(min(2.0 * Z - Z * Z, 1.0))   # worse  (XOR child)
            new.append(Z * Z)                         # better (copy child)
        channels = new
    return channels


# ---------------------------------------------------------------------------
# 6.  Frozen set computation for each bit-channel polar code
# ---------------------------------------------------------------------------

def compute_frozen_sets(
    K: int,
    p: float,
    q: float,
    N: int = 256,
) -> tuple[list[list[int]], list[list[int]], list[int]]:
    """Compute the frozen and information sets for each of the K component polar codes.

    For each bit level i:
      1. Compute the exact Bhattacharyya parameter Z(B^(i)_λ).
      2. Evolve Z through log2(N) stages of Arikan's polar transform to obtain
         N synthetic-channel Bhattacharyya parameters (reliability estimates).
      3. Set K_info^(i) = round(I(B^(i)_λ) · N)  (capacity rule).
      4. Select the K_info^(i) most reliable (smallest Z) synthetic channels
         as the information set; the rest are frozen.

    Parameters
    ----------
    K : int
        Number of bit levels (L = 2^K PPM alphabet size).
    p : float
        Correct-detection probability.
    q : float
        Erasure probability.
    N : int
        Component polar code block length (must be a power of 2).

    Returns
    -------
    info_sets   : list of K lists of ints (sorted synthetic-channel indices)
    frozen_sets : list of K lists of ints
    k_infos     : list of K ints, K_info^(i) = len(info_sets[i])
    """
    n_stages = int(round(np.log2(N)))
    assert (1 << n_stages) == N, "N must be a power of 2."

    caps   = bit_channel_capacities(K, p, q)
    z_bits = bit_channel_bhattacharyya(K, p, q)

    info_sets   = []
    frozen_sets = []
    k_infos     = []

    for i in range(K):
        k_info = int(round(caps[i] * N))
        # Evolve Z(B^i) through n_stages polar transform stages
        Z_synthetic = _polar_bhatt_evolution(z_bits[i], n_stages)
        # Rank synthetic channels by reliability (ascending Z = more reliable)
        ranked = sorted(range(N), key=lambda j: Z_synthetic[j])
        info_idx   = sorted(ranked[:k_info])
        frozen_idx = sorted(ranked[k_info:])

        info_sets.append(info_idx)
        frozen_sets.append(frozen_idx)
        k_infos.append(k_info)

    return info_sets, frozen_sets, k_infos


# ---------------------------------------------------------------------------
# 7.  Generate frozen masks as 256-bit integers (for SystemVerilog pkg)
# ---------------------------------------------------------------------------

def frozen_masks_as_int(frozen_sets: list[list[int]]) -> list[int]:
    """Convert frozen sets to 256-bit integer masks.

    mask[level][i] == 1  iff synthetic channel i is frozen (carries 0).
    Bit 0 = least-significant bit = synthetic channel 0 (always frozen;
    channel 0 is the worst after polarization).

    Parameters
    ----------
    frozen_sets : as returned by compute_frozen_sets().

    Returns
    -------
    masks : list of K ints, each fitting in 256 bits.
    """
    masks = []
    for fset in frozen_sets:
        mask = 0
        for idx in fset:
            mask |= (1 << idx)
        masks.append(mask)
    return masks


def print_sv_frozen_masks(masks: list[int], K: int) -> None:
    """Print the FROZEN_MASK localparam block for pkg_mlpolar.sv."""
    print("// -----------------------------------------------------------------------")
    print("// Frozen bit masks (256-bit per level) for SystemVerilog pkg_mlpolar.sv")
    print("// Generated by mlpolar_design.py")
    print("// -----------------------------------------------------------------------")
    print("localparam logic [255:0] FROZEN_MASK [0:K_LEVELS-1] = '{")
    for i, mask in enumerate(masks):
        comma = "," if i < K - 1 else ""
        print(f"  256'h{mask:064x}{comma}  // level {i}")
    print("};")


# ---------------------------------------------------------------------------
# 8.  Full design report
# ---------------------------------------------------------------------------

def design_report(K: int, p: float, q: float, N: int = 256) -> None:
    """Print a complete ML polar code design summary.

    Reproduces:
      - The SBP bit-channel capacities (matching the table in the system design)
      - The frozen masks in SystemVerilog format (matching pkg_mlpolar.sv)

    Parameters
    ----------
    K, p, q, N : as above.
    """
    L = 1 << K
    r = (1.0 - p - q) / (L - 1)

    print("=" * 72)
    print(f"ML Polar Code Design  |  K={K}, L={L}, N={N}")
    print(f"PPM channel: p={p}, q={q}, r={r:.4e}")
    print("=" * 72)

    # --- Channel capacity ---
    T   = ppm_channel(K, p, q)
    py  = T.mean(axis=0)
    def H(v): return float(-np.sum(v * np.log2(v, where=v > 0, out=np.zeros_like(v))))
    hy  = H(py)
    hyx = H(T[0])     # H(Y|X=0) = H(Y|X=x) for all x by symmetry
    IXY = hy - hyx
    print(f"\nChannel capacity  I(X;Y) = {IXY:.6f} bits/channel use")
    print(f"  H(Y)   = {hy:.6f} bits")
    print(f"  H(Y|X) = {hyx:.6f} bits")

    # --- SBP bit-channel capacities ---
    caps   = bit_channel_capacities(K, p, q)
    z_bits = bit_channel_bhattacharyya(K, p, q)
    k_infos_raw = [int(round(c * N)) for c in caps]
    total_info = sum(k_infos_raw)

    print(f"\nK-SBP bit-channel capacities (SP / natural-binary labeling):")
    print(f"  {'Level':>5}  {'I(B^i) [bits]':>14}  {'Z(B^i)':>10}  "
          f"{'K_info':>7}  {'Frozen':>7}")
    print("  " + "-" * 50)
    for i in range(K):
        print(f"  {i:5d}  {caps[i]:14.6f}  {z_bits[i]:10.6f}  "
              f"{k_infos_raw[i]:7d}  {N - k_infos_raw[i]:7d}")
    print(f"  {'SUM':>5}  {sum(caps):14.6f}  {'—':>10}  "
          f"{total_info:7d}  {K*N - total_info:7d}")
    print(f"\n  Chain-rule check: Σ I(B^i) = {sum(caps):.6f}, "
          f"I(X;Y) = {IXY:.6f},  match={np.isclose(sum(caps), IXY, atol=1e-6)}")
    print(f"\n  Total info bits / block: {total_info}")
    print(f"  Total coded bits / block: {K * N}")
    print(f"  Effective rate: {total_info / (K * N):.6f}  "
          f"(target = {IXY / K:.6f} bits/bit)")

    # --- Frozen sets ---
    info_sets, frozen_sets, k_infos = compute_frozen_sets(K, p, q, N)
    masks = frozen_masks_as_int(frozen_sets)

    print(f"\nFrozen sets (Bhattacharyya density evolution, N={N}):")
    print(f"  {'Level':>5}  {'K_info':>7}  {'first 5 info indices':>30}  "
          f"{'first 5 frozen indices':>25}")
    print("  " + "-" * 75)
    for i in range(K):
        fi = frozen_sets[i][:5]
        ii = info_sets[i][:5]
        print(f"  {i:5d}  {k_infos[i]:7d}  {str(ii):>30}  {str(fi):>25}")

    print()
    print_sv_frozen_masks(masks, K)

    # --- Transition matrices summary ---
    mats = ksbp_transition_matrices(K, p, q)
    print(f"\nK-SBP transition matrix shapes:")
    for i, M in enumerate(mats):
        print(f"  B^({i})_λ:  {M.shape}  (2 inputs, {M.shape[1]} = (L+1)·2^{i} outputs)")

    print("\n" + "=" * 72)


# ---------------------------------------------------------------------------
# 9.  Main: run for K=10, p=0.7, q=0.05
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    # -----------------------------------------------------------------------
    # Default parameters matching the ASIC design target
    # -----------------------------------------------------------------------
    K   = 10          # log2(L) = 10 bit levels
    p   = 0.70        # correct-detection probability
    q   = 0.05        # erasure probability
    N   = 256         # component polar code block length

    design_report(K, p, q, N)

    # -----------------------------------------------------------------------
    # Additionally: write out the transition matrices for inspection
    # -----------------------------------------------------------------------
    print("\nTransition matrix sample: B^(0)_λ (first 4 output columns):")
    mats = ksbp_transition_matrices(K, p, q)
    M0 = mats[0]
    print("  P(y=0 | b_0=0) =", M0[0, 0])
    print("  P(y=0 | b_0=1) =", M0[1, 0])
    print("  P(y=L | b_0=0) =", M0[0, (1<<K)])   # erasure column
    print("  P(y=L | b_0=1) =", M0[1, (1<<K)])

    # -----------------------------------------------------------------------
    # Verify chain rule explicitly
    # -----------------------------------------------------------------------
    caps = bit_channel_capacities(K, p, q)
    L    = 1 << K
    T    = ppm_channel(K, p, q)
    py   = T.mean(axis=0)
    def H(v): return float(-np.sum(v * np.log2(v, where=v>0, out=np.zeros_like(v))))
    IXY  = H(py) - H(T[0])
    assert np.isclose(sum(caps), IXY, atol=1e-6), \
        f"Chain rule FAILED: Σ caps = {sum(caps):.8f}, I(X;Y) = {IXY:.8f}"
    print(f"\n✓ Chain rule verified: Σ I(B^i) = {sum(caps):.6f} = I(X;Y) = {IXY:.6f}")
