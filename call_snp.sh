#!/bin/bash
#SBATCH --job-name=snakemake_job
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=56  # 使用 56 线程
#SBATCH --nodes=1
#SBATCH --time=48:00:00  # 根据实际需要调整
#SBATCH --mem=0  # 让 SLURM 自动分配所有内存
#SBATCH --output=slurm-%j.out

# 激活 conda
source ~/miniconda3/etc/profile.d/conda.sh

# 运行 snakemake
snakemake --cores 56 --use-conda
