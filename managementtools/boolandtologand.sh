#!/bin/bash
#
# Convert boolean operators from & to && and
# | to ||
# Converted files need manual review!
echo Converting $1 ...
sed 's/ \& / \&\& /g' $1 > /tmp/unixfile.m
sed 's/ | / || /g' < /tmp/unixfile.m > $1
