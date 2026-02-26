// =============================================================================
// ppm_llr_compute.sv
// Compute bit-level LLRs for the 1024-PPM channel under SP labeling
//
// INPUT: Photon detection result y ∈ {0,...,1023, 1024}
//   y ∈ {0..1023}: photon detected in slot y  (correct or wrong detection)
//   y = 1024:      no photon detected           (erasure event)
//
// The SP (natural binary) labeling maps PPM slot k to its K=10-bit binary
// representation [b_0,...,b_9] (MSB to LSB).
//
// MULTI-STAGE DECODING (MSD) LLRs (Seidl et al. Sec. II-A, eq. 3-4):
// At level i with side information b_0,...,b_{i-1} already decoded:
//
//   L^(i)(y, b_0,...,b_{i-1}) = log[ P(b_i=0 | y, b_0,...,b_{i-1}) ]
//                                   / [ P(b_i=1 | y, b_0,...,b_{i-1}) ]
//
// For the 1024-PPM channel with SP labeling, these are computed as:
//
//   P(b_i=v | y, prefix) ∝ Σ_{k: b_0(k)=prev_0,...,b_{i-1}(k)=prev_{i-1},
//                                   b_i(k)=v} P(y|k)
//
// where the transition probabilities are:
//   P(y|k) = p    if y==k (correct slot)
//           = q    if y==1024 (erasure)
//           = r    otherwise (wrong slot)
//
// KEY INSIGHT: For SP labeling of a symmetric channel, each bit level LLR
// depends only on whether the OBSERVED slot y falls in the "same half" as
// the input, given the decoded prefix. This allows efficient computation.
//
// HARDWARE STRUCTURE:
//   For each PPM symbol (10 bit levels):
//   1. Detect which SP half each possible k falls into at level i
//   2. Sum P(y|k) over both halves conditioned on decoded prefix
//   3. Compute log ratio → quantized LLR
//
// For the FPGA implementation we exploit the structure of P(y|k):
//   - If y=1024 (erasure): LLR=0 for all levels (symmetric erasure)
//   - If y=k: contributes p to the correct half at every level
//   - Otherwise: contributes r to both halves equally (cancels in LLR)
//
// This simplifies massively: the LLR for bit i given observation y (non-erasure)
// is determined entirely by whether bit i of y equals the decoded half.
// =============================================================================

`include "pkg_mlpolar.sv"

module ppm_llr_compute
  import pkg_mlpolar::*;
(
  input  logic              clk,
  input  logic              rst_n,
  input  logic              valid_in,

  // Raw channel observation: which PPM slot was detected
  // Encoded as 10-bit slot index (0..1023) or 1024 for erasure
  input  logic [K_LEVELS:0] obs_slot,     // 11 bits: 0..1024

  // Previously decoded bits from MSD (bit levels 0..i-1 are known)
  // Updated by msd_controller as each level decodes
  input  logic [K_LEVELS-1:0] decoded_prefix, // bits decoded so far
  input  logic [3:0]           current_level,  // which level is requesting LLRs

  // Output: LLR for bit 'current_level' of each of the N=256 component code
  // positions. Since all 256 positions share the SAME PPM observation
  // (one PPM symbol per time slot), this outputs a single LLR value that
  // is broadcast to all N=256 decoder positions for this level.
  //
  // NOTE: In MLC, the N=256 component code bits at each level are
  // transmitted across N=256 SEPARATE PPM symbols. So we output one LLR
  // per PPM symbol → the decoder buffers all N=256 LLRs per level.
  // This module computes ONE LLR per cycle (one PPM symbol per clock).
  output llr_t             llr_out,
  output logic             valid_out
);

  // -------------------------------------------------------------------------
  // LLR computation for 1024-PPM with SP labeling under MSD
  //
  // Given observation y and decoded prefix b_{0:i-1}, the LLR for bit i is:
  //
  //   Numerator:   sum_{k in C0} P(y|k)  where C0 = {k : b_0(k)..b_{i-1}(k)
  //                                                   = prefix, b_i(k)=0}
  //   Denominator: sum_{k in C1} P(y|k)  where C1 = same but b_i(k)=1
  //
  // |C0| = |C1| = 2^(K-i-1) slots
  //
  // Transition probabilities:
  //   P(y|k) = p if y==k, q if y==1024, r otherwise  (r << p, r << q)
  //
  // ERASURE (y=1024):
  //   P(y=1024|k) = q for ALL k → ratio = 1 → LLR = 0
  //
  // CORRECT SLOT (y ∈ {0..1023}):
  //   One k* = y has P(y|k*) = p, all others have P(y|k) = r
  //   The key contribution of p comes from k* = y.
  //   k* ∈ C0 iff bits 0..i-1 of y == prefix AND bit i of y == 0
  //   k* ∈ C1 iff bits 0..i-1 of y == prefix AND bit i of y == 1
  //   k* ∉ C0 or C1 iff prefix bits of y ≠ decoded_prefix
  //
  // Case A: prefix bits of y match decoded_prefix AND bit i of y == 0:
  //   Num = p + (|C0|-1)*r,  Den = |C1|*r
  //   LLR ≈ log((p + (|C0|-1)*r) / (|C1|*r))
  //        ≈ log(p / (|C1|*r))  since p >> r
  //        = log(p) - log(|C1|) - log(r)
  //
  // Case B: prefix bits match AND bit i of y == 1: (symmetric) LLR → negative
  //
  // Case C: prefix bits of y DON'T match decoded_prefix:
  //   k* ∉ C0 ∪ C1 → both halves contribute r * |C0| = r * |C1|
  //   Num = |C0|*r,  Den = |C1|*r → LLR = 0  (no information from this obs)
  //
  // IMPLEMENTATION: We precompute the LLR magnitudes for each level using
  // the formula log(p / (2^(K-i-1) * r)), quantized to 8-bit signed.
  // The sign is determined by bit i of the observed slot (after prefix check).
  // -------------------------------------------------------------------------

  // Precomputed LLR magnitudes per level (offline computed):
  // LLR_mag[i] = round(log2(p / (2^(K-i-1) * r)) * SCALE)
  // With p=0.7, r=2.444e-4, scale to 8-bit (factor 4 → 0.25 bit resolution):
  //   i=0: log2(0.7 / (512  * 2.444e-4)) = log2(5.59)  = 2.48 → *4 = 9.9  ≈ 10
  //   i=1: log2(0.7 / (256  * 2.444e-4)) = log2(11.18) = 3.48 → *4 = 13.9 ≈ 14
  //   i=2: log2(0.7 / (128  * 2.444e-4)) = log2(22.35) = 4.48 → *4 = 17.9 ≈ 18
  //   i=3: log2(0.7 / (64   * 2.444e-4)) = log2(44.70) = 5.48 → *4 = 21.9 ≈ 22
  //   i=4: log2(0.7 / (32   * 2.444e-4)) = log2(89.40) = 6.48 → *4 = 25.9 ≈ 26
  //   i=5: log2(0.7 / (16   * 2.444e-4)) = log2(178.8) = 7.48 → *4 = 29.9 ≈ 30
  //   i=6: log2(0.7 / (8    * 2.444e-4)) = log2(357.6) = 8.48 → *4 = 33.9 ≈ 34
  //   i=7: log2(0.7 / (4    * 2.444e-4)) = log2(715.2) = 9.48 → *4 = 37.9 ≈ 38
  //   i=8: log2(0.7 / (2    * 2.444e-4)) = log2(1430)  = 10.48→ *4 = 41.9 ≈ 42
  //   i=9: log2(0.7 / (1    * 2.444e-4)) = log2(2864)  = 11.48→ *4 = 45.9 ≈ 46
  //   (clipped to LLR_MAX=64 for saturation)

  localparam logic signed [LLR_BITS-1:0] LLR_MAG [0:K_LEVELS-1] = '{
    8'sd10, 8'sd14, 8'sd18, 8'sd22, 8'sd26,
    8'sd30, 8'sd34, 8'sd38, 8'sd42, 8'sd46
  };

  // -------------------------------------------------------------------------
  // Pipeline registers
  // -------------------------------------------------------------------------
  logic [K_LEVELS:0]    obs_r;
  logic [K_LEVELS-1:0]  prefix_r;
  logic [3:0]           level_r;
  logic                 valid_r;

  always_ff @(posedge clk or negedge rst_n) begin
    if (!rst_n) begin
      obs_r    <= '0; prefix_r <= '0;
      level_r  <= '0; valid_r  <= 1'b0;
    end else begin
      obs_r    <= obs_slot;
      prefix_r <= decoded_prefix;
      level_r  <= current_level;
      valid_r  <= valid_in;
    end
  end

  // -------------------------------------------------------------------------
  // Combinational LLR computation (1-cycle latency after input register)
  // -------------------------------------------------------------------------
  logic                  is_erasure;
  logic                  prefix_match;
  logic                  obs_bit_i;        // bit 'level' of observed slot
  llr_t                  llr_comb;

  always_comb begin
    // Erasure detection
    is_erasure = obs_r[K_LEVELS];  // obs_r == 1024

    // Check if prefix bits of observed slot match decoded prefix
    // prefix bits = bits [K_LEVELS-1 : current_level+1] of obs_r
    // For level i, we check bits K-1 downto i+1
    // In SP labeling, bit j of slot index k = (k >> (K-1-j)) & 1
    // i.e., bit 0 = MSB of 10-bit slot index
    prefix_match = 1'b1;
    for (int j = 0; j < K_LEVELS; j++) begin
      if (j < int'(level_r)) begin
        // Check bit j of obs_r against decoded_prefix[j]
        // Bit j of slot: obs_r[K_LEVELS-1-j]
        if (obs_r[K_LEVELS-1-j] !== prefix_r[K_LEVELS-1-j])
          prefix_match = 1'b0;
      end
    end

    // Bit at current level from observation
    obs_bit_i = obs_r[K_LEVELS-1-level_r];

    // Compute LLR
    if (is_erasure) begin
      // Erasure: symmetric → LLR = 0
      llr_comb = '0;
    end else if (!prefix_match) begin
      // Observed slot not in conditioning set → no information → LLR = 0
      llr_comb = '0;
    end else begin
      // Slot in conditioning set: high-reliability LLR
      // Sign: positive if obs_bit_i==0 (b_i=0 hypothesis favored), negative if ==1
      if (!obs_bit_i)
        llr_comb = LLR_MAG[level_r];   // positive: b_i=0 likely
      else
        llr_comb = -LLR_MAG[level_r];  // negative: b_i=1 likely
    end
  end

  always_ff @(posedge clk or negedge rst_n) begin
    if (!rst_n) begin
      llr_out   <= '0;
      valid_out <= 1'b0;
    end else begin
      llr_out   <= llr_comb;
      valid_out <= valid_r;
    end
  end

endmodule
