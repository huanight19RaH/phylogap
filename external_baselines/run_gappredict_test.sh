#!/bin/bash
set -e
export TF_USE_LEGACY_KERAS=1
export PYTHONPATH=/mnt/d/THStudy/UniversityStudy/Research/phylogap/external_baselines/GapPredict/lib
mkdir -p /mnt/d/THStudy/UniversityStudy/Research/phylogap/external_baselines/GapPredict/test_output
cd /mnt/d/THStudy/UniversityStudy/Research/phylogap/external_baselines/GapPredict/lib
~/.local/share/mamba/envs/gappredict_env/bin/python GapPredict.py \
  -o /mnt/d/THStudy/UniversityStudy/Research/phylogap/external_baselines/GapPredict/test_output/ \
  -fa /mnt/d/THStudy/UniversityStudy/Research/phylogap/external_baselines/GapPredict/lib/data/real_gaps/sealer_filled/7391826_358-1408.fasta \
  -fq /mnt/d/THStudy/UniversityStudy/Research/phylogap/external_baselines/GapPredict/lib/data/real_gaps/sealer_filled/7391826_358-1408.fastq \
  -e 2 \
  -gpu -1
