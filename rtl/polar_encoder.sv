// =============================================================================
// polar_encoder.sv
// N=256 polar encoder: computes c = u * G_N where G_N = B_N * F_N
//
// Architecture: 8-stage combinational butterfly network (no registers).
// The generator matrix G_N = B_N * F_2^{⊗8} implements the polar transform.
//
// Each "butterfly" stage l (0..7) performs the kernel F_2 = [[1,0],[1,1]]:
//   x_out[2i]   = x_in[2i] XOR x_in[2i+1]   (check node: worse channel)
//   x_out[2i+1] = x_in[2i+1]                 (copy:  better channel)
// applied to pairs separated by stride 2^l.
//
// After all 8 stages, a bit-reversal permutation B_N is applied.
//
// In the MLC framework (Seidl et al. Sec. IV-B, eq. 54-55):
//   This encoder implements the N-SBP pi^n for one component code.
//   The labeling L_lambda maps bit-level outputs to PPM slot indices.
//   All K=10 component encoders run in parallel (one per bit level).
//
// For iCE40: N=256 XOR butterfly = 8 stages × 128 XOR gates = 1024 LUTs
//   Each stage is a combinational pipeline stage behind a register layer.
// =============================================================================

`include "pkg_mlpolar.sv"

module polar_encoder
  import pkg_mlpolar::*;
(
  input  logic              clk,
  input  logic              rst_n,
  // Control
  input  logic              valid_in,    // u[] is valid
  output logic              valid_out,   // c[] is valid (8 cycles later)
  // Data: u = info+frozen bits in natural order (frozen=0)
  input  logic [N-1:0]      u,
  // Codeword output (bit-reversed generator output)
  output logic [N-1:0]      c
);

  // -------------------------------------------------------------------------
  // Butterfly pipeline: 8 stages, each registered for timing closure on iCE40
  // Stage s operates on stride = 2^s
  // -------------------------------------------------------------------------
  logic [N-1:0] stage [0:LOG2N];  // stage[0]=u, stage[8]=pre-bitrev output
  logic         valid_pipe [0:LOG2N];

  // Stage 0: input
  always_ff @(posedge clk or negedge rst_n) begin
    if (!rst_n) begin
      stage[0]      <= '0;
      valid_pipe[0] <= 1'b0;
    end else begin
      stage[0]      <= u;
      valid_pipe[0] <= valid_in;
    end
  end

  // Generate butterfly stages 1..8
  // At stage s (producing stage[s] from stage[s-1]):
  //   stride = 2^(s-1)  (stride in the butterfly graph)
  //   Within each block of size 2*stride, XOR pairs at distance stride
  genvar s, i;
  generate
    for (s = 1; s <= LOG2N; s++) begin : STAGE
      logic [N-1:0] stg_comb;
      localparam int STRIDE = 1 << (s-1);

      // Combinational butterfly for this stage
      for (i = 0; i < N; i++) begin : BUTTERFLY
        // Within each block of 2*STRIDE elements, bit i:
        // If (i % (2*STRIDE)) < STRIDE: it is the "XOR output" position
        //   → receives XOR of itself and the element STRIDE positions higher
        // Otherwise: copy
        if ((i % (2*STRIDE)) < STRIDE) begin
          assign stg_comb[i] = stage[s-1][i] ^ stage[s-1][i + STRIDE];
        end else begin
          assign stg_comb[i] = stage[s-1][i];
        end
      end

      always_ff @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
          stage[s]      <= '0;
          valid_pipe[s] <= 1'b0;
        end else begin
          stage[s]      <= stg_comb;
          valid_pipe[s] <= valid_pipe[s-1];
        end
      end
    end
  endgenerate

  // -------------------------------------------------------------------------
  // Bit-reversal permutation B_N (applied after all butterfly stages)
  // B_N reverses the LOG2N-bit index of each element.
  // B_N is its own inverse: B_N^2 = I
  // Example for N=8: index 6 (110) -> 3 (011)
  // -------------------------------------------------------------------------
  function automatic int bit_reverse(input int idx, input int nbits);
    int result;
    result = 0;
    for (int b = 0; b < nbits; b++)
      result |= ((idx >> b) & 1) << (nbits - 1 - b);
    return result;
  endfunction

  generate
    for (i = 0; i < N; i++) begin : BITREV
      assign c[i] = stage[LOG2N][bit_reverse(i, LOG2N)];
    end
  endgenerate

  assign valid_out = valid_pipe[LOG2N];

endmodule
