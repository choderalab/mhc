#!/bin/bash
# set directory
echo "working directory set"
cd ~/boehm/hladr1-tpi/run-yank-scripts

# Run the simulation
echo "Running simulation..."
yank script --yaml=1klu.yaml

# Analyze the data
#echo "Analyzing data..."
#yank analyze --store=../output-files/1klu
