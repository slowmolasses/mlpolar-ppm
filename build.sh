#!/usr/bin/env bash
# =============================================================================
# build.sh
# Full FPGA Build Flow: SystemVerilog → Yosys → nextpnr-ice40 → iceprog
#
# Requirements:
#   sudo apt install yosys nextpnr-ice40 icestorm iverilog
#   (or build from source: https://github.com/YosysHQ/oss-cad-suite-build)
#
# Usage:
#   ./build.sh          # synthesize + P&R + bitstream + program
#   ./build.sh synth    # synthesis only
#   ./build.sh sim      # run testbench simulation (iverilog)
#   ./build.sh clean    # remove build artifacts
#
# Target: iCEBreaker board (iCE40UP5K-SG48)
# =============================================================================

set -e  # exit on first error

PROJ=mlpolar_ppm
DEVICE=up5k
PACKAGE=sg48
TOP=top
PCF=mlpolar_ppm.pcf

# Source files (order matters for yosys: package first)
SV_FILES=(
  pkg_mlpolar.sv
  polar_encoder.sv
  ppm_llr_compute.sv
  polar_decoder_sc.sv
  msd_controller.sv
  ml_polar_top.sv
)

# ─────────────────────────────────────────────────────────────────────────────
# STEP 0: Tool version check
# ─────────────────────────────────────────────────────────────────────────────
check_tools() {
  echo "=== Checking tools ==="
  for tool in yosys nextpnr-ice40 icepack iceprog iverilog; do
    if command -v $tool &>/dev/null; then
      echo "  ✓ $tool: $(${tool} --version 2>&1 | head -1)"
    else
      echo "  ✗ $tool NOT FOUND"
      echo "    Install: sudo apt install yosys nextpnr-ice40 icestorm iverilog"
      echo "    Or use OSS CAD Suite: https://github.com/YosysHQ/oss-cad-suite-build"
      exit 1
    fi
  done
}

# ─────────────────────────────────────────────────────────────────────────────
# STEP 1: Simulation with iverilog (SystemVerilog via -g2012)
# ─────────────────────────────────────────────────────────────────────────────
sim() {
  echo ""
  echo "=== STEP 1: Simulation (iverilog + vvp) ==="
  mkdir -p build

  # Compile testbench + DUT
  iverilog -g2012 -Wall \
    -I . \
    "${SV_FILES[@]}" \
    tb_mlpolar.sv \
    -o build/${PROJ}_sim.vvp

  # Run simulation
  vvp build/${PROJ}_sim.vvp | tee build/sim.log
  echo "  → Simulation log: build/sim.log"
  echo "  → Waveform: tb_mlpolar.vcd (open with gtkwave tb_mlpolar.vcd)"
}

# ─────────────────────────────────────────────────────────────────────────────
# STEP 2: Synthesis with Yosys
#
# Key yosys commands explained:
#   read_verilog -sv   : Parse SystemVerilog (SV2012 subset)
#   synth_ice40        : iCE40-specific synthesis pass:
#     - proc            : Convert processes (always blocks) to RTL
#     - flatten         : Inline all module hierarchies
#     - tribuf          : Handle tri-state buffers
#     - alumacc         : Map adder/multiplier to SB_MAC16
#     - share           : Resource sharing optimization
#     - opt             : General logic optimization
#     - fsm             : FSM extraction and re-encoding
#     - techmap         : Map to iCE40 technology cells (SB_LUT4, SB_FF, etc.)
#     - abc             : Berkeley ABC for logic optimization and LUT mapping
#     - ice40_opt       : iCE40-specific post-ABC optimization
#     - autoname/stat   : Naming and statistics
#   write_json         : Output netlist for nextpnr
# ─────────────────────────────────────────────────────────────────────────────
synth() {
  echo ""
  echo "=== STEP 2: Synthesis (yosys) ==="
  mkdir -p build

  yosys -p "
    # Read all source files as SystemVerilog
    read_verilog -sv -I. $(echo "${SV_FILES[@]}")

    # Elaborate design from top module
    hierarchy -check -top ${TOP}

    # Synthesize for iCE40 UP5K
    # -top: specify top module
    # -json: output netlist format for nextpnr
    # -abc2: use ABC optimization (better for data paths)
    # -dsp: attempt to use iCE40 DSP blocks (SB_MAC16) for multiplications
    synth_ice40 -top ${TOP} -json build/${PROJ}.json -abc2

    # Print resource utilization report
    stat
  " 2>&1 | tee build/yosys.log

  echo "  → Netlist: build/${PROJ}.json"
  echo "  → Yosys log: build/yosys.log"

  # Extract LUT count from log
  LUT_COUNT=$(grep "SB_LUT4" build/yosys.log | grep -oP '\d+' | tail -1)
  echo "  → SB_LUT4 count: ${LUT_COUNT:-unknown}"
}

# ─────────────────────────────────────────────────────────────────────────────
# STEP 3: Place and Route with nextpnr-ice40
#
# nextpnr is a timing-driven FPGA P&R tool.
# Key options:
#   --up5k          : iCE40 UP5K device (5280 LUTs, 30 EBR BRAMs, 8 DSPs)
#     alternatives:   --hx8k (iCE40HX8K), --lp8k, --u4k
#   --package sg48  : SG48 package (48-pin QFN, used on iCEBreaker)
#   --json          : Synthesized netlist from yosys
#   --pcf           : Physical constraints (pin assignments)
#   --asc           : Output bitstream in ASCII format
#   --freq 24       : Target clock frequency 24 MHz (for timing analysis)
#   --seed          : Random seed for P&R (change if timing fails)
#   --placer heap   : Use HeAP analytic placer (good for dense designs)
#   --router router2: nextpnr's second-generation router (better congestion)
#   --timing-allow-fail: Don't abort if timing constraints not met (for dev)
# ─────────────────────────────────────────────────────────────────────────────
pnr() {
  echo ""
  echo "=== STEP 3: Place and Route (nextpnr-ice40) ==="

  nextpnr-ice40 \
    --${DEVICE} \
    --package ${PACKAGE} \
    --json build/${PROJ}.json \
    --pcf ${PCF} \
    --asc build/${PROJ}.asc \
    --freq 24 \
    --placer heap \
    --router router2 \
    --seed 42 \
    2>&1 | tee build/nextpnr.log

  echo "  → ASCII bitstream: build/${PROJ}.asc"
  echo "  → nextpnr log: build/nextpnr.log"

  # Extract timing summary
  grep -A3 "Max frequency" build/nextpnr.log || true
}

# ─────────────────────────────────────────────────────────────────────────────
# STEP 4: Bitstream packing with icepack (part of Project IceStorm)
#
# icepack converts nextpnr's ASCII bitstream to binary .bin format
# that iceprog and the iCE40 boot ROM understand.
# ─────────────────────────────────────────────────────────────────────────────
pack() {
  echo ""
  echo "=== STEP 4: Bitstream Packing (icepack) ==="
  icepack build/${PROJ}.asc build/${PROJ}.bin
  echo "  → Binary bitstream: build/${PROJ}.bin"
  ls -lh build/${PROJ}.bin
}

# ─────────────────────────────────────────────────────────────────────────────
# STEP 5: Program the FPGA with iceprog
#
# iceprog uses the FTDI FT2232H USB bridge on iCEBreaker to:
#   1. Send bitstream over SPI to the iCE40's configuration flash (AT25SF081)
#   2. Assert CRESET_B to force the iCE40 to reconfigure from flash
#
# Requires: USB rules - add user to 'plugdev' group, or run with sudo
#   sudo sh -c 'echo "ATTRS{idVendor}==\"0403\", ATTRS{idProduct}==\"6010\",
#     MODE=\"0660\", GROUP=\"plugdev\"" > /etc/udev/rules.d/53-ftdi.rules'
#   sudo udevadm control --reload-rules && sudo udevadm trigger
#
# iceprog options:
#   -d i:0x0403:0x6010  : Use FTDI device (iCEBreaker VID:PID)
#   (default)           : Write to flash address 0x000000
#   -o 0x100000         : Write to flash offset (for multi-image)
# ─────────────────────────────────────────────────────────────────────────────
program() {
  echo ""
  echo "=== STEP 5: Programming FPGA (iceprog) ==="
  echo "  Ensure iCEBreaker is connected via USB..."

  # Check device is present
  if ! lsusb | grep -q "0403:6010"; then
    echo "  ERROR: FTDI USB device not found. Check USB connection."
    exit 1
  fi

  iceprog build/${PROJ}.bin
  echo "  → FPGA programmed successfully!"
}

# ─────────────────────────────────────────────────────────────────────────────
# STEP 6: Timing analysis (optional, post-P&R)
# ─────────────────────────────────────────────────────────────────────────────
timing() {
  echo ""
  echo "=== Timing Analysis ==="
  icetime -d ${DEVICE} -mtr build/${PROJ}_timing.rpt build/${PROJ}.asc
  cat build/${PROJ}_timing.rpt
}

# ─────────────────────────────────────────────────────────────────────────────
# Clean
# ─────────────────────────────────────────────────────────────────────────────
clean() {
  echo "=== Cleaning build artifacts ==="
  rm -rf build/ tb_mlpolar.vcd
}

# ─────────────────────────────────────────────────────────────────────────────
# Main dispatch
# ─────────────────────────────────────────────────────────────────────────────
case "${1:-all}" in
  sim)     check_tools; sim ;;
  synth)   check_tools; synth ;;
  pnr)     check_tools; pnr ;;
  pack)    check_tools; pack ;;
  program) check_tools; pack; program ;;
  timing)  check_tools; timing ;;
  clean)   clean ;;
  all)
    check_tools
    sim
    synth
    pnr
    pack
    program
    ;;
  *)
    echo "Usage: $0 [all|sim|synth|pnr|pack|program|timing|clean]"
    exit 1
    ;;
esac

echo ""
echo "=== Build complete ==="
