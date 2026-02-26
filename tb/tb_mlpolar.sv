// =============================================================================
// tb_mlpolar.sv
// Testbench: Encode a random block, transmit over 1024-PPM channel,
// decode with ML polar SC decoder, check bit error rate.
// =============================================================================

`timescale 1ns/1ps
`include "pkg_mlpolar.sv"

module tb_mlpolar;
  import pkg_mlpolar::*;

  // Clock and reset
  logic clk = 0;
  logic rst_n;
  always #20.833 clk = ~clk;   // 24 MHz

  // DUT I/O
  logic        enc_bit_in, enc_valid_in, enc_ready;
  logic [9:0]  ppm_slot_out;
  logic        ppm_valid_out;
  logic [10:0] obs_in;
  logic        obs_valid, obs_ready;
  logic        dec_bit_out, dec_valid, dec_block_done;

  ml_polar_top dut (
    .clk_12mhz    (clk),
    .rst_n        (rst_n),
    .enc_bit_in   (enc_bit_in),
    .enc_valid_in (enc_valid_in),
    .enc_ready    (enc_ready),
    .ppm_slot_out (ppm_slot_out),
    .ppm_valid_out(ppm_valid_out),
    .ppm_ready_in (1'b1),
    .obs_in       (obs_in),
    .obs_valid    (obs_valid),
    .obs_ready    (obs_ready),
    .dec_bit_out  (dec_bit_out),
    .dec_valid    (dec_valid),
    .dec_block_done(dec_block_done)
  );

  // -------------------------------------------------------------------------
  // Random info bits
  // -------------------------------------------------------------------------
  localparam int TOTAL_INFO = 1588;
  logic [TOTAL_INFO-1:0] tx_bits;
  logic [TOTAL_INFO-1:0] rx_bits;
  int rx_cnt;

  // -------------------------------------------------------------------------
  // Channel simulation: 1024-PPM
  // -------------------------------------------------------------------------
  real p_det = 0.70;   // correct detection
  real q_era = 0.05;   // erasure
  // r = (1-p-q)/(L-1)

  function automatic logic [10:0] ppm_channel(input logic [9:0] tx_slot);
    real rnd;
    rnd = $urandom_range(0, 65535) / 65535.0;
    if (rnd < p_det)
      return {1'b0, tx_slot};         // correct detection
    else if (rnd < p_det + q_era)
      return 11'd1024;                 // erasure
    else begin
      // wrong slot: uniform over L-1 others
      logic [9:0] wrong;
      do wrong = $urandom_range(0, 1023); while (wrong == tx_slot);
      return {1'b0, wrong};
    end
  endfunction

  // -------------------------------------------------------------------------
  // Test procedure
  // -------------------------------------------------------------------------
  int bit_errors = 0;
  int frame_cnt  = 0;

  initial begin
    $dumpfile("tb_mlpolar.vcd");
    $dumpvars(0, tb_mlpolar);

    rst_n        = 0;
    enc_valid_in = 0;
    obs_valid    = 0;
    repeat(4) @(posedge clk);
    rst_n = 1;
    repeat(2) @(posedge clk);

    // ── Encode and transmit 3 blocks ──────────────────────────────────────
    for (int blk = 0; blk < 3; blk++) begin
      $display("=== BLOCK %0d ===", blk);

      // Generate random info bits
      for (int i = 0; i < TOTAL_INFO; i++)
        tx_bits[i] = $urandom_range(0,1);

      // Feed info bits to encoder
      enc_valid_in = 1;
      for (int i = 0; i < TOTAL_INFO; i++) begin
        @(posedge clk);
        while (!enc_ready) @(posedge clk);
        enc_bit_in = tx_bits[i];
      end
      enc_valid_in = 0;

      // Wait for PPM symbols to appear, channel-corrupt, feed to decoder
      obs_valid = 0;
      rx_cnt = 0;

      fork
        // Thread 1: capture PPM output, apply channel, feed to decoder
        begin
          int sym_cnt = 0;
          while (sym_cnt < N) begin
            @(posedge clk);
            if (ppm_valid_out) begin
              logic [10:0] obs;
              obs = ppm_channel(ppm_slot_out);
              // Feed observation to decoder
              @(posedge clk);
              obs_in    = obs;
              obs_valid = 1;
              @(posedge clk);
              obs_valid = 0;
              sym_cnt++;
            end
          end
        end

        // Thread 2: collect decoded bits
        begin
          while (!dec_block_done) begin
            @(posedge clk);
            if (dec_valid) begin
              rx_bits[rx_cnt] = dec_bit_out;
              rx_cnt++;
            end
          end
        end
      join

      // Count bit errors
      for (int i = 0; i < TOTAL_INFO; i++) begin
        if (tx_bits[i] !== rx_bits[i]) begin
          bit_errors++;
          $display("  BIT ERROR at position %0d: tx=%b rx=%b", i, tx_bits[i], rx_bits[i]);
        end
      end

      frame_cnt++;
      $display("  Block %0d: rx_cnt=%0d, errors so far=%0d", blk, rx_cnt, bit_errors);
      $display("  BER = %0.4f", real'(bit_errors) / real'(frame_cnt * TOTAL_INFO));
    end

    $display("\n=== FINAL RESULTS ===");
    $display("Frames: %0d, Total info bits: %0d", frame_cnt, frame_cnt * TOTAL_INFO);
    $display("Bit errors: %0d", bit_errors);
    $display("BER: %0.6f", real'(bit_errors) / real'(frame_cnt * TOTAL_INFO));

    #1000;
    $finish;
  end

  // Timeout
  initial begin
    #10_000_000;
    $display("TIMEOUT");
    $finish;
  end

endmodule
