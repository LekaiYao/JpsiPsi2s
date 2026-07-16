#!/bin/bash
# run_MC_all.sh — produce UnweightMC/Unweight_*.root for all SPS/DPS channels.
# Usage: after cmssw-el7 && cmsenv, run from NtupleAnalyzer directory:
#   bash run_MC_all.sh

CHANNELS=(
    DPS
    ggpsi1psi2
    ggpsi1psi2g
    ggpsi1psi2g_gpt0p4
    ggpsi1psi2g_gpt0p8
    ggpsi1psi2g_gpt1p2
    ggpsi1psi2g_gpt1p6
    ggpsi1psi2g_gpt2p0
    ggpsi1psi2g_gpt2p4
)

echo "=== run_MC_all.sh ==="
echo "Channels: ${#CHANNELS[@]}"
echo "Output dir: UnweightMC/"
echo ""

# Enable MC_MODE (uncomment the line)
sed -i 's|^// #define MC_MODE|#define MC_MODE|' rephrase.cpp

for ch in "${CHANNELS[@]}"; do
    echo "--- $ch ---"
    sed -i "s|string MC_CHANNEL = .*;|string MC_CHANNEL = \"$ch\";|" rephrase.cpp
    root -l -q -b rephrase.cpp
    echo ""
done

# Restore Data mode (re-comment MC_MODE)
sed -i 's|^#define MC_MODE|// #define MC_MODE|' rephrase.cpp

echo "=== All done ==="
echo "Output files:"
ls -la UnweightMC/
