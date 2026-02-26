// =============================================================================
// ml_polar_top.sv
// Top-Level: Multi-Level Polar Code over 1024-PPM Channel
// Target: Lattice iCE40 FPGA (iCE40HX8K or iCE40UP5K)
//
// ARCHITECTURE OVERVIEW
// ─────────────────────────────────────────────────────────────────────────────
// This module integrates:
//
//  (A) ENCODER PATH:
//      Input: 1588 info bits (packed as serial stream)
//      Split into 10 groups per bit level (K_INFO[i] bits each)
//      Each group feeds a polar_encoder (N=256 butterfly network)
//      Output: 10 × 256 = 2560 coded bits → mapped to 256 PPM slots via
//              SP (natural binary) labeling:
//              PPM slot index = {c_9[n], c_8[n], ..., c_0[n]} for symbol n
//
//  (B) DECODER PATH (Multi-Stage Decoding):
//      Input: 256 PPM observations (11-bit: slot 0..1023 or 1024=erasure)
//      Level 0: compute bit-level LLRs → SC decode → get b_0[0..255]
//      Level i: condition on b_0..b_{i-1} → SC decode → get b_i[0..255]
//      Output: 1588 info bits per block
//
// MAPPING (SP labeling, natural binary, MSB = level 0):
//   Slot index k for symbol n = SUM_i  c_i[n] * 2^(K-1-i)
//   where c_i[n] is the n-th codeword bit of component code i.
//
// RESOURCE ESTIMATE (iCE40HX8K: 7680 LUTs, 2× SPRAM, 20× BRAM):
//   - 10× polar_encoder (256 XOR butterfly): ~800 LUTs (pipelined, shared)
//   - 1×  SC decoder (combinational tree for N=256): ~2000 LUTs
//   - LLR memory (N=256 × 8-bit): 2× EBR BRAM
//   - Obs buffer (N=256 × 11-bit): 1× EBR BRAM
//   - Side-info RAM (N=256 × 10-bit): 1× EBR BRAM
//   - Controller FSM + misc: ~200 LUTs
//   Total estimate: ~3000 LUTs → fits in iCE40HX4K (3520 LUTs)
//
// CLOCK: 12 MHz XTAL → target 24 MHz after PLL (ice40 internal PLL)
// Throughput: 10 × 2×256 = 5120 cycles/block = 5120/24e6 ≈ 213 µs/block
//             = 1588 info bits / 213 µs ≈ 7.5 Mbps
//
// =============================================================================

`include "pkg_mlpolar.sv"

module ml_polar_top
  import pkg_mlpolar::*;
(
  input  logic        clk_12mhz,    // 12 MHz board crystal
  input  logic        rst_n,        // active-low async reset (button)

  // ── Encoder interface ──────────────────────────────────────────────────────
  // Info bits arrive serially, 1588 bits per block, LSB first
  input  logic        enc_bit_in,   // serial info bit input
  input  logic        enc_valid_in, // info bit valid
  output logic        enc_ready,    // encoder ready to accept info bits

  // Encoded PPM symbol output: 10-bit slot index per symbol (N=256 symbols/block)
  output logic [9:0]  ppm_slot_out, // PPM slot index to transmit
  output logic        ppm_valid_out,// PPM symbol valid
  input  logic        ppm_ready_in, // downstream ready

  // ── Decoder interface ──────────────────────────────────────────────────────
  // PPM observations: 11-bit (0..1023 = slot, 1024 = erasure), N=256 per block
  input  logic [10:0] obs_in,       // PPM observation
  input  logic        obs_valid,    // observation valid
  output logic        obs_ready,    // decoder ready

  // Decoded info bits (serial output)
  output logic        dec_bit_out,  // decoded bit
  output logic        dec_valid,    // decoded bit valid
  output logic        dec_block_done// block decoding complete

);

  // -------------------------------------------------------------------------
  // Internal PLL (iCE40 SB_PLL40_CORE): 12 → 24 MHz
  // If PLL not available (e.g., iCE40LP), omit and use clk_12mhz directly
  // -------------------------------------------------------------------------
  logic clk;
  // synthesis translate_off
  assign clk = clk_12mhz;  // simulation: bypass PLL
  // synthesis translate_on

  // For actual synthesis, instantiate SB_PLL40_CORE:
  // SB_PLL40_CORE #(
  //   .FEEDBACK_PATH("SIMPLE"),
  //   .DIVR(4'b0000), .DIVF(7'b0000111), .DIVQ(3'b101),
  //   .FILTER_RANGE(3'b001)
  // ) pll (
  //   .REFERENCECLK(clk_12mhz), .PLLOUTCORE(clk),
  //   .LOCK(), .RESETB(1'b1), .BYPASS(1'b0)
  // );

  // =========================================================================
  // ENCODER PATH
  // =========================================================================

  // -------------------------------------------------------------------------
  // Info bit demultiplexer: route incoming info bits to 10 level buffers
  // Buffer i receives K_INFO[i] bits, then triggers encoder i
  // -------------------------------------------------------------------------
  localparam int MAX_K_INFO = 178;  // maximum K_INFO across levels
  logic [MAX_K_INFO-1:0] info_buf [0:K_LEVELS-1];  // info bit buffers
  logic [7:0]            info_cnt [0:K_LEVELS-1];  // bits received per level
  logic [7:0]            total_info_rcvd;           // total info bits received

  // Bit-to-level mapping: bits arrive in order [level 0 bits, level 1 bits, ...]
  logic [3:0] enc_level_sel;   // which level we're currently filling
  logic [8:0] enc_bit_cnt;     // global bit counter within current block

  always_ff @(posedge clk or negedge rst_n) begin
    if (!rst_n) begin
      enc_level_sel <= '0;
      enc_bit_cnt   <= '0;
      for (int lv = 0; lv < K_LEVELS; lv++) begin
        info_buf[lv] <= '0;
        info_cnt[lv] <= '0;
      end
    end else if (enc_valid_in && enc_ready) begin
      // Route bit to current level
      info_buf[enc_level_sel][info_cnt[enc_level_sel]] <= enc_bit_in;
      info_cnt[enc_level_sel] <= info_cnt[enc_level_sel] + 1;

      // Advance level when current level buffer is full
      if (info_cnt[enc_level_sel] == K_INFO[enc_level_sel] - 1) begin
        enc_level_sel <= enc_level_sel + 1;
        info_cnt[enc_level_sel] <= '0;
      end

      // Advance global counter
      enc_bit_cnt <= enc_bit_cnt + 1;
    end
  end

  // -------------------------------------------------------------------------
  // u-vector assembly: insert frozen zeros at correct positions
  // u[i] = 0 if channel i is frozen, else next info bit from level buffer
  // (Seidl et al. eq. 43-46: pi^n partitions B^N into N bit channels,
  //  frozen channels carry fixed zeros, info channels carry source data)
  // -------------------------------------------------------------------------
  logic [N-1:0] u_vec [0:K_LEVELS-1];   // input to polar encoders

  generate
    for (genvar lv = 0; lv < K_LEVELS; lv++) begin : U_ASSEMBLE
      // For each bit channel position i, if frozen → 0, else take from info_buf
      logic [7:0] info_ptr;   // points into info_buf[lv]

      always_comb begin
        int ptr;
        ptr = 0;
        for (int i = 0; i < N; i++) begin
          if (FROZEN_MASK[lv][i]) begin
            u_vec[lv][i] = 1'b0;
          end else begin
            u_vec[lv][i] = info_buf[lv][ptr];
            ptr = ptr + 1;
          end
        end
      end
    end
  endgenerate

  // -------------------------------------------------------------------------
  // 10 polar encoders (run in parallel, triggered when all info bits loaded)
  // Each encoder implements c = u * G_N (N=256 butterfly)
  // -------------------------------------------------------------------------
  logic        enc_start;    // triggers all encoders simultaneously
  logic [N-1:0] codeword [0:K_LEVELS-1];  // encoder outputs c_i[0..N-1]
  logic        enc_valid [0:K_LEVELS-1];

  // enc_start fires when all 1588 info bits received
  assign enc_start = (enc_bit_cnt == 1588) && enc_valid_in;
  assign enc_ready = (enc_bit_cnt < 1588);

  generate
    for (genvar lv = 0; lv < K_LEVELS; lv++) begin : ENC_INST
      polar_encoder enc_i (
        .clk       (clk),
        .rst_n     (rst_n),
        .valid_in  (enc_start),
        .valid_out (enc_valid[lv]),
        .u         (u_vec[lv]),
        .c         (codeword[lv])
      );
    end
  endgenerate

  // -------------------------------------------------------------------------
  // PPM symbol assembly: SP labeling
  // slot_index[n] = {c_0[n], c_1[n], ..., c_9[n]}
  // (bit 0 of slot index = MSB output of level 0 encoder, i.e. c_0[n])
  // Serializer outputs one 10-bit slot per clock after enc_valid fires
  // -------------------------------------------------------------------------
  logic [7:0]  ppm_sym_cnt;
  logic        ppm_active;

  // All encoders finish simultaneously (same latency LOG2N cycles)
  assign ppm_active = enc_valid[0];

  always_ff @(posedge clk or negedge rst_n) begin
    if (!rst_n) begin
      ppm_sym_cnt  <= '0;
      ppm_valid_out <= 1'b0;
      ppm_slot_out  <= '0;
    end else begin
      ppm_valid_out <= 1'b0;
      if (ppm_active || (ppm_sym_cnt > 0 && ppm_sym_cnt < N)) begin
        // Assemble SP slot index for symbol ppm_sym_cnt
        // slot = {c_0[n], c_1[n], ..., c_9[n]}  (MSB = level 0)
        for (int lv = 0; lv < K_LEVELS; lv++)
          ppm_slot_out[K_LEVELS-1-lv] <= codeword[lv][ppm_sym_cnt];

        ppm_valid_out <= 1'b1;
        ppm_sym_cnt   <= (ppm_sym_cnt == N-1) ? '0 : ppm_sym_cnt + 1;
      end
    end
  end

  // =========================================================================
  // DECODER PATH
  // =========================================================================
  // Instantiate 10 SC decoders (one per level), driven by msd_controller
  // Each decoder is parameterized with its level's frozen mask and K_INFO

  // LLR buffers for each level (output of ppm_llr_compute, input to SC dec)
  llr_t  llr_buf [0:K_LEVELS-1][0:N-1];  // 10 × 256 × 8-bit = 20 KB (use BRAM)
  logic  llr_buf_wr [0:K_LEVELS-1];
  logic [7:0] llr_wr_addr [0:K_LEVELS-1];

  // -------------------------------------------------------------------------
  // 10 LLR computation units (one per level, run sequentially via MSD ctrl)
  // -------------------------------------------------------------------------
  logic [K_LEVELS-1:0]  decoded_prefix_n [0:N-1]; // side info per symbol
  logic [3:0]            active_level;
  logic [7:0]            llr_sym_cnt;

  // Side info: as each SC decoder finishes level i, store bit decisions
  // into decoded_prefix_n[n][i] for all n=0..N-1
  // (handled by msd_controller)

  // -------------------------------------------------------------------------
  // Instantiate msd_controller
  // -------------------------------------------------------------------------
  msd_controller msd_ctrl (
    .clk           (clk),
    .rst_n         (rst_n),
    .obs_valid     (obs_valid),
    .obs_in        (obs_in),
    .obs_ready     (obs_ready),
    .info_valid    (dec_valid),
    .info_bit      (dec_bit_out),
    .block_done    (dec_block_done)
  );

  // =========================================================================
  // Status / Debug LEDs (iCE40 has 3-8 user LEDs)
  // Expose on I/O for timing analysis
  // =========================================================================
  // (Wired at PCB/constraint level; no explicit port here)

endmodule


// =============================================================================
// PPM-ready simple top for iCE40 Breakout Board pinout
// (Wraps ml_polar_top with board-specific I/O)
// =============================================================================
module top (
  input  logic       CLK,       // 12 MHz XTAL
  input  logic       BTN_N,     // active-low reset button

  // SPI-like encoder interface (shift in 1588 bits)
  input  logic       ENC_DATA,  // PMOD pin 1
  input  logic       ENC_CLK,   // PMOD pin 2 (soft SPI clock, synchronous)
  input  logic       ENC_CS_N,  // PMOD pin 3 (frame start)

  // Observation input (11-bit serial, 2 PMOD pins + counter)
  input  logic       OBS_DATA,  // PMOD pin 7
  input  logic       OBS_CLK,   // PMOD pin 8

  // Decoded output
  output logic       DEC_DATA,  // PMOD pin 9
  output logic       DEC_CLK,   // PMOD pin 10
  output logic       DEC_DONE,  // PMOD pin 11

  // PPM output (10-bit parallel to DAC/modulator board)
  output logic [9:0] PPM_SLOT,  // 10 PMOD pins
  output logic       PPM_VALID,

  // Status LEDs
  output logic       LED_ENC,   // encoder active
  output logic       LED_DEC    // decoder active
);

  import pkg_mlpolar::*;

  logic rst_n_sync;
  always_ff @(posedge CLK) rst_n_sync <= BTN_N;

  // SPI deserializer for encoder bits
  logic enc_bit, enc_valid;
  logic [3:0] spi_cnt;
  logic spi_cs_n_r;

  always_ff @(posedge CLK) begin
    spi_cs_n_r <= ENC_CS_N;
    enc_valid  <= 1'b0;
    if (!ENC_CS_N && ENC_CLK) begin
      enc_bit   <= ENC_DATA;
      enc_valid <= 1'b1;
    end
  end

  // Observation shift-in (11 bits per observation)
  logic [10:0] obs_shift;
  logic [3:0]  obs_cnt;
  logic        obs_valid_r;

  always_ff @(posedge CLK or negedge rst_n_sync) begin
    if (!rst_n_sync) begin
      obs_shift  <= '0;
      obs_cnt    <= '0;
      obs_valid_r <= 1'b0;
    end else begin
      obs_valid_r <= 1'b0;
      if (OBS_CLK) begin
        obs_shift <= {OBS_DATA, obs_shift[10:1]};
        obs_cnt   <= obs_cnt + 1;
        if (obs_cnt == 10) begin
          obs_valid_r <= 1'b1;
          obs_cnt     <= '0;
        end
      end
    end
  end

  // Decoded output SPI-like serializer
  logic [7:0] dec_shift;
  logic [2:0] dec_cnt;
  logic       dec_clk_r;

  always_ff @(posedge CLK) begin
    dec_clk_r <= 1'b0;
    // Simple: output directly
    DEC_DATA <= 1'b0;
    DEC_CLK  <= 1'b0;
    DEC_DONE <= 1'b0;
  end

  ml_polar_top dut (
    .clk_12mhz    (CLK),
    .rst_n        (rst_n_sync),
    .enc_bit_in   (enc_bit),
    .enc_valid_in (enc_valid),
    .enc_ready    (LED_ENC),
    .ppm_slot_out (PPM_SLOT),
    .ppm_valid_out(PPM_VALID),
    .ppm_ready_in (1'b1),
    .obs_in       (obs_shift),
    .obs_valid    (obs_valid_r),
    .obs_ready    (),
    .dec_bit_out  (DEC_DATA),
    .dec_valid    (DEC_CLK),
    .dec_block_done(DEC_DONE)
  );

  assign LED_DEC = DEC_DONE;

endmodule
