#!/bin/bash

#
# 
#

# Run the simulation
echo "Running simulation..."
yank script --yaml=hladr1-tpi-explicit.yaml

# Analyze the data
#echo "Analyzing data..."
#yank analyze --store=hladr1-tpi-explicit-output
