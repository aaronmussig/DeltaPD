#!/bin/zsh

set -e

cd /Users/aaron/phd/DeltaPDNew

# mamba activate deltapd

maturin develop --release
