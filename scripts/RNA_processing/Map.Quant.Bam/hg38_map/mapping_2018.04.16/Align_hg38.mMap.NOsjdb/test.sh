#!/bin/bash
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -A lsens2018-3-3
#SBATCH -t 03:00:00
#SBATCH -J 256.m.noS.STAR
#SBATCH -o 256.m.noS.STAR.out
#SBATCH -e 256.m.noS.STAR.err

touch test.worked.txt
