#!/bin/bash
#SBATCH -n1
#SBATCH -c20
#SBATCH --time=00:45:00
#SBATCH --mail-user=wleeb@math.princeton.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --mem=100000
#
matlab -nodisplay -r "run11; quit"
