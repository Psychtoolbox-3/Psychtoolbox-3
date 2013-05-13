#!/bin/bash
# Convert | to || and & to && - driver. Just iterates over all
# M-Files and calls boolandtologand.sh for each file.
# Usage: ./convertlogicals.sh <Path to Psychtoolbox folder.
find $1 -name "*.m" -exec ./boolandtologand.sh {} \;
echo Done.
