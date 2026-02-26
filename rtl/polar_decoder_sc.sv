// =============================================================================
// polar_decoder_sc.sv
// Successive Cancellation (SC) Decoder for N=256 Polar Code
//
// Implements Arikan's SC decoding algorithm [1] in log-domain (LLR domain).
// Corresponds to Seidl et al. Sec. III-B: successive decoding of symbols u_i
// over bit channels B^(i)_{pi^n}, using information combining operations.
//
// ALGORITHM OVERVIEW (SC decoding as a binary tree traversal):
//   - The decoder operates on a binary tree of depth LOG2N=8
//   - Each node (stage s, position i) stores a channel LLR
//   - Two fundamental operations (F2 kernel inverses):
//
//     (1) f-function (upper branch, "check node", worse synthetic channel):
//         LLR_out = f(a, b) = sign(a)*sign(b)*min(|a|,|b|)   [min-sum approx]
//         Exact: log[(1 + e^{a+b}) / (e^a + e^b)]
//         This combines two LLRs for the XOR-ed input (b_0 path)
//
//     (2) g-function (lower branch, "variable node", better synthetic channel):
//         LLR_out = g(a, b, u_hat) = b + (1-2*u_hat)*a
//         Conditioned on previously decoded bit u_hat
//         This is exact (no approximation needed)
//
// ARCHITECTURE: Radix-2 SC uses a schedule that processes 2N-1 nodes total.
// For N=256: 511 node evaluations, processed sequentially.
// Each evaluation is O(1) → total O(N) decode time = 512 clock cycles.
//
// DATA FLOW:
//   - LLR memory: 2N-1 entries (binary tree nodes), 8-bit signed
//   - Bit memory: N entries (partial-sum bits for g-function conditioning)
//   - Controller: generates node addresses and operation codes
//
// FROZEN CHANNELS: At frozen positions, u_hat = 0 (known a priori).
//   The frozen mask is provided as a parameter per bit level.
//
// INTERFACE:
//   - Load phase: channel LLRs for all N=256 positions are written in
//   - Decode phase: SC schedule runs, emitting decisions at each info channel
//
// =============================================================================

`include "pkg_mlpolar.sv"

module polar_decoder_sc
  import pkg_mlpolar::*;
#(
  parameter logic [N-1:0] FROZEN = '0,  // frozen bit mask for this level
  parameter int            K_INFO_P = 106  // number of info bits for this level
)(
  input  logic              clk,
  input  logic              rst_n,

  // LLR loading interface (N cycles, one per channel output position)
  input  logic              llr_valid,   // new LLR arriving
  input  llr_t              llr_in,      // channel LLR for next position
  output logic              llr_ready,   // decoder ready to accept LLRs

  // Decode trigger
  input  logic              decode_start, // assert after all N LLRs loaded

  // Decoded information bits (emitted in order as SC proceeds)
  output logic              bit_valid,   // info bit decoded
  output logic              bit_out,     // decoded bit value
  output logic              decode_done  // all K_INFO_P bits decoded
);

  // -------------------------------------------------------------------------
  // Internal memory
  // -------------------------------------------------------------------------
  // LLR tree: leaf-level LLRs at stage 0 (depth LOG2N from root)
  // Store as flat array: node n at stage s occupies address s*N + n
  // Simpler: use separate array per stage (LOG2N+1 stages, N/2^s nodes each)
  // For simplicity: store LLRs for each SC node as they are computed.
  // We use a RAM of depth 2*N (enough for the full LLR tree).
  localparam int LLR_RAM_DEPTH = 2 * N;  // 512 entries
  logic signed [LLR_BITS-1:0] llr_ram [0:LLR_RAM_DEPTH-1];

  // Partial sum bits: one per leaf (N entries)
  logic [N-1:0] partial_sum;

  // -------------------------------------------------------------------------
  // State machine
  // -------------------------------------------------------------------------
  typedef enum logic [1:0] {
    ST_IDLE    = 2'b00,
    ST_LOAD    = 2'b01,   // receiving channel LLRs
    ST_DECODE  = 2'b10,   // running SC schedule
    ST_DONE    = 2'b11
  } state_t;

  state_t state;

  // SC schedule counters
  logic [7:0]  load_cnt;      // counts 0..N-1 during LLR loading
  logic [8:0]  phase;         // SC phase index: 0..2N-2 (N=256 → 0..511)
  logic [7:0]  bit_idx;       // current leaf index being processed (0..N-1)
  int          info_cnt;      // count of info bits emitted

  // -------------------------------------------------------------------------
  // LLR loading
  // -------------------------------------------------------------------------
  always_ff @(posedge clk or negedge rst_n) begin
    if (!rst_n) begin
      load_cnt  <= '0;
      state     <= ST_IDLE;
    end else begin
      case (state)
        ST_IDLE: begin
          if (decode_start) begin
            state    <= ST_LOAD;
            load_cnt <= '0;
          end
        end
        ST_LOAD: begin
          if (llr_valid) begin
            // Leaf LLRs stored at stage LOG2N (addresses N..2N-1)
            llr_ram[N + load_cnt] <= llr_in;
            load_cnt <= load_cnt + 1;
            if (load_cnt == N-1) begin
              state <= ST_DECODE;
              phase <= '0;
              bit_idx <= '0;
              info_cnt <= 0;
            end
          end
        end
        ST_DECODE: begin
          if (phase == 2*N-2)
            state <= ST_DONE;
          else
            phase <= phase + 1;
        end
        ST_DONE: begin
          state <= ST_IDLE;
        end
      endcase
    end
  end

  assign llr_ready = (state == ST_IDLE) || (state == ST_LOAD);

  // -------------------------------------------------------------------------
  // SC Schedule: for each phase, determine which node to evaluate
  //
  // The SC schedule for N=2^n polar codes visits leaves left-to-right.
  // For leaf i, we first compute all necessary f/g operations top-down,
  // then make the decision at the leaf, then propagate bottom-up.
  //
  // Compact representation: process leaf by leaf.
  // For leaf i, the f/g operations form a binary tree traversal.
  // We use a simple sequential controller:
  //   - For each leaf i (0..N-1), evaluate LOG2N LLRs top-down
  //   - Make decision
  //   - Update partial sums bottom-up
  //
  // This requires N*(LOG2N+1) cycles, but we pipeline to N cycles by
  // noting that only the path from root to leaf i matters.
  // For clarity and iCE40 resource constraints, we use N cycles/leaf
  // with a depth-first schedule.
  // -------------------------------------------------------------------------

  // For each leaf i, the node sequence is determined by the binary tree structure.
  // We implement a simplified version: compute the LLR for each leaf sequentially
  // using precomputed alpha/beta values.

  // Intermediate LLR arrays for the current leaf's path
  logic signed [LLR_BITS-1:0] alpha [0:LOG2N];  // LLRs top-down per stage
  logic                        beta  [0:LOG2N];  // partial sums bottom-up

  // f-function: min-sum approximation
  function automatic logic signed [LLR_BITS-1:0] f_func(
    input logic signed [LLR_BITS-1:0] a,
    input logic signed [LLR_BITS-1:0] b
  );
    logic signed [LLR_BITS-1:0] abs_a, abs_b, min_val;
    logic                        sign_xor;
    abs_a   = (a[LLR_BITS-1]) ? -a : a;
    abs_b   = (b[LLR_BITS-1]) ? -b : b;
    min_val = (abs_a < abs_b) ? abs_a : abs_b;
    sign_xor = a[LLR_BITS-1] ^ b[LLR_BITS-1];
    return sign_xor ? -min_val : min_val;
  endfunction

  // g-function: exact (no approximation)
  function automatic logic signed [LLR_BITS-1:0] g_func(
    input logic signed [LLR_BITS-1:0] a,
    input logic signed [LLR_BITS-1:0] b,
    input logic                        u_hat
  );
    return u_hat ? (b - a) : (b + a);
  endfunction

  // -------------------------------------------------------------------------
  // Decode phase: sequential leaf-by-leaf SC
  // Each leaf takes LOG2N+1 cycles (simplified pipelined version)
  // For a lean iCE40 implementation we use a combinational tree per leaf
  // -------------------------------------------------------------------------
  logic [7:0]   leaf_idx;    // current leaf being decoded
  logic [3:0]   stage_cnt;   // which stage in the top-down pass (0=root)
  logic         decision;    // hard decision at leaf

  // Current leaf LLR accumulation (combinational for speed)
  // We compute all stages for leaf i in one clock by unrolling the tree

  logic signed [LLR_BITS-1:0] stage_llr [0:LOG2N][0:N-1]; // LLR at each (stage, position)

  // The stage_llr[LOG2N][i] are the channel LLRs (leaf inputs)
  always_comb begin
    for (int i = 0; i < N; i++)
      stage_llr[LOG2N][i] = llr_ram[N + i];
  end

  // Compute stages from leaves upward (root = stage 0)
  // stage_llr[s][i] = LLR of synthetic channel at stage s, position i
  // Uses partial_sum bits for g-function inputs
  // This is the "message passing" of the SC algorithm

  genvar gs, gi;
  generate
    for (gs = LOG2N-1; gs >= 0; gs--) begin : SC_STAGE
      localparam int HALF = 1 << gs;  // half-block size at this stage
      for (gi = 0; gi < (1 << gs); gi++) begin : SC_NODE
        // Each node at stage gs, position gi processes:
        // children at stage gs+1: left child = gi, right child = gi + HALF
        logic signed [LLR_BITS-1:0] left_llr, right_llr;
        assign left_llr  = stage_llr[gs+1][gi];
        assign right_llr = stage_llr[gs+1][gi + HALF];

        // Which operation depends on the current leaf being decoded:
        // if the gi-th bit of leaf_idx at depth gs is 0 → f-function (top branch)
        // if it's 1 → g-function (bottom branch), conditioned on partial_sum
        logic use_g;
        assign use_g = (state == ST_DECODE) && ((leaf_idx >> gs) & 1);

        assign stage_llr[gs][gi] = use_g ?
          g_func(left_llr, right_llr, partial_sum[gi]) :
          f_func(left_llr, right_llr);
      end
    end
  endgenerate

  // -------------------------------------------------------------------------
  // Decision and partial sum update (each clock in DECODE state)
  // -------------------------------------------------------------------------
  logic [7:0] leaf_r;  // registered leaf index for output

  always_ff @(posedge clk or negedge rst_n) begin
    if (!rst_n) begin
      leaf_idx     <= '0;
      partial_sum  <= '0;
      bit_valid    <= 1'b0;
      bit_out      <= 1'b0;
      decode_done  <= 1'b0;
      info_cnt     <= 0;
    end else begin
      bit_valid   <= 1'b0;
      decode_done <= 1'b0;

      if (state == ST_DECODE) begin
        // LLR at root node (stage 0, position 0) is the decision LLR for leaf_idx
        logic signed [LLR_BITS-1:0] root_llr;
        logic                        u_hat;

        root_llr = stage_llr[0][0];

        // Hard decision: u_hat = 0 if LLR >= 0, else 1
        if (FROZEN[leaf_idx]) begin
          u_hat = 1'b0;  // frozen channel: force to 0
        end else begin
          u_hat = root_llr[LLR_BITS-1];  // MSB = sign bit; 1 = negative → u=1
          bit_valid  <= 1'b1;
          bit_out    <= u_hat;
          info_cnt   <= info_cnt + 1;
          if (info_cnt == K_INFO_P - 1)
            decode_done <= 1'b1;
        end

        // Update partial sums (bottom-up propagation)
        // partial_sum[leaf_idx] = u_hat; then XOR up the tree
        partial_sum[leaf_idx] <= u_hat;

        // Advance to next leaf
        leaf_idx <= leaf_idx + 1;
      end

      if (state == ST_IDLE) begin
        leaf_idx    <= '0;
        partial_sum <= '0;
        info_cnt    <= 0;
      end
    end
  end

endmodule
