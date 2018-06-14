#!/bin/bash

#BATCH --job-name 2-7
#SBATCH --partition=titanx-short
#SBATCH --sockets-per-node=1
#SBATCH --cores-per-socket=1
#SBATCH --time 0-4:00:00
#SBATCH --output=res_%j.txt  # output file
#SBATCH -e res_%j.err        # File to which STDERR will be written
#SBATCH --gres=gpu:4
#SBATCH --mem=0

srun python models_tgcn_tanh_attn.py > ./finaloutput-422.txt
