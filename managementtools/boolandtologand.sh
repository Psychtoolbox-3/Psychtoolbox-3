#!/bin/bash
#
# Convert boolean operators from & to && and
# | to ||
# Converted files need manual review!
echo Converting $1 ...

# Solution used for conversion:
#sed 's/ \& / \&\& /g' $1 > /tmp/unixfile.m
#sed 's/ | / || /g' < /tmp/unixfile.m > $1

# Tobias solution, probably more robust, but untested on OSX:
sed -i -e 's/\(\b\|\s\)\&\(\b\|\s\)/ \&\& /g;s/\(\b\|\s\)|\(\b\|\s\)/ || /g' $1
