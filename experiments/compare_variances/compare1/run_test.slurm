#!/bin/bash
#SBATCH -n1
#SBATCH -c40
#SBATCH --time=01:00:00
#SBATCH --mail-user=wleeb@math.princeton.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --mem=100000
#
matlab -nodisplay -r "run12; quit"
