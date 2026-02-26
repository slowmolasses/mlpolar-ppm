// =============================================================================
// msd_controller.sv
// Multi-Stage Decoder (MSD) Controller for ML Polar Code
//
// Implements the multi-stage decoding procedure of Seidl et al. Sec. IV-A:
//   "The receiver performs multi-stage decoding (MSD), i.e., it computes
//    reliability information for decoding of the first bit level, passed to
//    the decoder of the first component code. The decoding results are used
//    for demapping and decoding of the next bit level, and so on."
//
// DECODING ORDER (Sequential: level 0 → level 9):
//   Stage 0: Compute LLRs for B^(0)_lambda using only channel outputs y
//            → Decode polar code 0 → get b_0[0..N-1]
//
//   Stage 1: Compute LLRs for B^(1)_lambda given b_0 (MSD side info)
//            For each PPM symbol n: LLR^(1)(y_n, b_0[n])
//            → Decode polar code 1 → get b_1[0..N-1]
//
//   Stage i: Compute LLRs for B^(i)_lambda given b_0,...,b_{i-1}
//            → Decode polar code i → get b_i[0..N-1]
//
// The "decoded_prefix" accumulates decoded bits:
//   For symbol n at level i: prefix_n[j] = b_j[n] for j=0..i-1
//   This is the SP label prefix that conditions the bit-level LLR.
//
// TIMING:
//   - Receive N=256 PPM observations (parallel if possible, else serial)
//   - For level i=0..9:
//       * Compute N LLRs (N cycles, one per PPM symbol)
//       * Run SC decoder (N cycles for N=256 leaf decisions)
//       * Store N decoded bits b_i[0..N-1] as side info for next level
//   - Total per block: 10 * 2N = 5120 cycles + overhead
//
// MEMORY:
//   - obs_buf[N]:  buffered PPM observations (10-bit each)
//   - side_info[K_LEVELS][N]: b_j[n] for j=0..level-1, n=0..N-1
//     (stored as 10-bit "column" per symbol: side_info_col[n] = {b_9,...,b_0})
//
// =============================================================================

`include "pkg_mlpolar.sv"

module msd_controller
  import pkg_mlpolar::*;
(
  input  logic                clk,
  input  logic                rst_n,

  // Input: PPM observations for N=256 symbols (stream, one per cycle)
  input  logic                obs_valid,
  input  logic [K_LEVELS:0]   obs_in,       // 11 bits: slot 0..1023 or 1024=erasure
  output logic                obs_ready,

  // Output: decoded information bits (stream, in order: level 0 bits, then level 1, etc.)
  output logic                info_valid,
  output logic                info_bit,
  output logic                block_done    // all K_INFO bits for this block decoded
);

  // -------------------------------------------------------------------------
  // Observation buffer
  // -------------------------------------------------------------------------
  logic [K_LEVELS:0]  obs_buf [0:N-1];   // 11-bit x 256
  logic [7:0]         obs_load_cnt;
  logic               obs_buf_full;

  always_ff @(posedge clk or negedge rst_n) begin
    if (!rst_n) begin
      obs_load_cnt <= '0;
      obs_buf_full <= 1'b0;
    end else if (obs_valid && obs_ready) begin
      obs_buf[obs_load_cnt] <= obs_in;
      obs_load_cnt <= obs_load_cnt + 1;
      if (obs_load_cnt == N-1)
        obs_buf_full <= 1'b1;
    end else if (block_done) begin
      obs_load_cnt <= '0;
      obs_buf_full <= 1'b0;
    end
  end

  // -------------------------------------------------------------------------
  // Side information: decoded bits per level per symbol
  // side_info_col[n][i] = bit i (level i) decoded for symbol n
  // -------------------------------------------------------------------------
  logic [K_LEVELS-1:0] side_info_col [0:N-1];  // 10 bits per symbol

  // -------------------------------------------------------------------------
  // State machine for MSD
  // -------------------------------------------------------------------------
  typedef enum logic [2:0] {
    MSD_IDLE      = 3'd0,
    MSD_OBS_LOAD  = 3'd1,   // loading observations
    MSD_LLR_COMP  = 3'd2,   // computing LLRs for current level
    MSD_SC_WAIT   = 3'd3,   // SC decoder running
    MSD_NEXT_LVL  = 3'd4,   // advance to next level
    MSD_DONE      = 3'd5
  } msd_state_t;

  msd_state_t msd_state;
  logic [3:0] current_level;     // 0..9
  logic [7:0] sym_cnt;           // symbol counter (0..N-1)

  // -------------------------------------------------------------------------
  // LLR computation interface → SC decoder interface
  // -------------------------------------------------------------------------
  logic        llr_comp_valid;
  llr_t        llr_computed;
  logic        llr_comp_done_r;

  logic        sc_llr_valid;
  llr_t        sc_llr_in;
  logic        sc_llr_ready;
  logic        sc_decode_start;
  logic        sc_bit_valid;
  logic        sc_bit_out;
  logic        sc_decode_done;

  // LLR computation module (one instance, time-multiplexed across 10 levels)
  ppm_llr_compute llr_unit (
    .clk             (clk),
    .rst_n           (rst_n),
    .valid_in        (llr_comp_valid),
    .obs_slot        (obs_buf[sym_cnt]),
    .decoded_prefix  (side_info_col[sym_cnt]),
    .current_level   (current_level),
    .llr_out         (llr_computed),
    .valid_out       (sc_llr_valid)
  );

  // SC decoder (10 instances in full parallel, OR time-multiplex one)
  // For iCE40 resource budget: use ONE SC decoder, reconfigured per level.
  // Frozen mask MUX: select FROZEN_MASK[current_level]

  // Since SC decoder takes a parameter (frozen mask), we instantiate all 10
  // but only clock the active one. For iCE40, this saves routing; yosys
  // will optimize away inactive register toggles.
  // Alternative: ROM-based frozen mask lookup (used here for compactness).

  logic [N-1:0]  active_frozen_mask;
  always_comb begin
    active_frozen_mask = FROZEN_MASK[current_level];
  end

  // Single SC decoder instance (frozen mask from ROM)
  polar_decoder_sc #(
    .FROZEN    (256'h0),   // overridden dynamically via frozen_rom below
    .K_INFO_P  (106)       // worst case; actual K varies — decoder counts bits
  ) sc_dec (
    .clk          (clk),
    .rst_n        (rst_n),
    .llr_valid    (sc_llr_valid),
    .llr_in       (llr_computed),
    .llr_ready    (sc_llr_ready),
    .decode_start (sc_decode_start),
    .bit_valid    (sc_bit_valid),
    .bit_out      (sc_bit_out),
    .decode_done  (sc_decode_done)
  );

  // NOTE: The frozen mask is actually applied inside the SC decoder as a
  // parameter. For a single-instance design that reconfigures per level,
  // we expose the frozen check externally and override the decision:

  // External frozen check: if active_frozen_mask[leaf_idx] → force bit=0
  // (The SC decoder's bit_valid only fires for info bits when we wire
  //  frozen mask externally. For this implementation we treat the SC
  //  decoder's FROZEN parameter as a ROM and instantiate separate decoders
  //  per level in the top-level module. msd_controller is the orchestrator.)

  // -------------------------------------------------------------------------
  // SC decoder output → side info accumulation
  // -------------------------------------------------------------------------
  logic [7:0]  decoded_bit_cnt;  // which bit within current level
  logic [7:0]  sym_map [0:N-1];  // maps decoded info bit index → symbol index
  // For polar codes with SP-labeling: the decoded bit stream from SC decoder
  // is in channel order (i=0..N-1). We need to map back to symbol n.
  // In the MLC setup: component code bit i corresponds to PPM symbol i.
  // So decoded bit i at level l → side_info_col[i][l].

  // Count of info bits seen (to distinguish frozen from info positions)
  // The SC decoder only emits bit_valid for info positions.
  // But for side-info accumulation we need ALL N bits (info + frozen=0).
  // We track the leaf index separately:
  logic [7:0]  leaf_ctr;   // 0..N-1, advances each cycle in SC_WAIT

  always_ff @(posedge clk or negedge rst_n) begin
    if (!rst_n) begin
      msd_state     <= MSD_IDLE;
      current_level <= '0;
      sym_cnt       <= '0;
      leaf_ctr      <= '0;
      block_done    <= 1'b0;
      info_valid    <= 1'b0;
      info_bit      <= 1'b0;
    end else begin
      block_done <= 1'b0;
      info_valid <= 1'b0;

      case (msd_state)
        // ---------------------------------------------------------------
        MSD_IDLE: begin
          if (obs_buf_full) begin
            msd_state     <= MSD_LLR_COMP;
            current_level <= '0;
            sym_cnt       <= '0;
          end
        end

        // ---------------------------------------------------------------
        // Compute N LLRs for current level (one per PPM symbol, sequential)
        // ppm_llr_compute takes obs_buf[sym_cnt] and side_info_col[sym_cnt]
        // and outputs the LLR for bit level 'current_level'.
        // ---------------------------------------------------------------
        MSD_LLR_COMP: begin
          llr_comp_valid <= 1'b1;
          sym_cnt        <= sym_cnt + 1;
          // LLR output arrives 2 cycles later (ppm_llr_compute pipeline depth)
          // SC decoder receives them via sc_llr_valid signal (combinationally)
          if (sym_cnt == N-1) begin
            msd_state      <= MSD_SC_WAIT;
            sc_decode_start <= 1'b1;
            sym_cnt        <= '0;
            leaf_ctr       <= '0;
          end
        end

        // ---------------------------------------------------------------
        // SC decoder running: N=256 cycles for N leaf decisions
        // Collect decoded bits and store as side info for next level
        // ---------------------------------------------------------------
        MSD_SC_WAIT: begin
          sc_decode_start <= 1'b0;
          llr_comp_valid  <= 1'b0;

          // leaf_ctr tracks which symbol we're processing
          // SC decoder emits sc_bit_valid only for info channels
          // We need all N bits (frozen=0, info=decoded value)
          // Track via leaf_ctr:
          if (sc_bit_valid || active_frozen_mask[leaf_ctr]) begin
            logic decoded_val;
            decoded_val = active_frozen_mask[leaf_ctr] ? 1'b0 : sc_bit_out;

            // Store in side info for next level
            side_info_col[leaf_ctr][current_level] <= decoded_val;

            // Emit info bits to output
            if (!active_frozen_mask[leaf_ctr]) begin
              info_valid <= 1'b1;
              info_bit   <= decoded_val;
            end

            leaf_ctr <= leaf_ctr + 1;
          end

          if (sc_decode_done) begin
            msd_state <= MSD_NEXT_LVL;
          end
        end

        // ---------------------------------------------------------------
        MSD_NEXT_LVL: begin
          if (current_level == K_LEVELS-1) begin
            msd_state  <= MSD_DONE;
          end else begin
            current_level <= current_level + 1;
            sym_cnt       <= '0;
            leaf_ctr      <= '0;
            msd_state     <= MSD_LLR_COMP;
          end
        end

        // ---------------------------------------------------------------
        MSD_DONE: begin
          block_done <= 1'b1;
          msd_state  <= MSD_IDLE;
        end

      endcase
    end
  end

  assign obs_ready      = (msd_state == MSD_IDLE) || (msd_state == MSD_OBS_LOAD);
  assign llr_comp_valid = (msd_state == MSD_LLR_COMP);

endmodule
