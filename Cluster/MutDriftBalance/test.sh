#!/bin/bash

tar -xzf slimbuild.tar.gz

./build/slim -d fmin=0 -d fmax=1 -d npops=10 -d nloci=1000 -d popsize=1000 -d scaleT0=0 -d scales=0 -d seed=$1 -d fitnessFunction=$2 -d a=$3 -d s=$4 -d r=$5 -d b=$6 -d mu=$7 -d std=$8 MutDrift.slim | tail -n +14 > ff=$2_seed=$1_a=$3_s=$4_r=$5_b=$6_mu=$7_std=$8.csv







