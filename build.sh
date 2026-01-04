#!/bin/bash
# Build script to avoid conda compiler conflicts

# Unset conda environment variables that interfere with compilation
unset CC
unset CXX
unset AR
unset CFLAGS
unset CXXFLAGS
unset LDFLAGS

# Use system compiler explicitly
export CC=gcc
export CXX=g++

# Clean and build
cargo clean
cargo build "$@"
