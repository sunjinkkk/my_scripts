# my_scripts
我用于生物信息学分析的一些脚本。
## fasta_stat.py
用于统计fasta格式文件
```bash
python fasta_stat.py -a genome.fa -o result.txt
```
## call_snp.smk
使用snakemake写的call snp pipeline
```bash
snakemake --cores 56 --use-conda
```
## call_snp.sh
用于提交到hpc的call_snp shell
```bash
bash call_snp.sh
```
## gtf_extract.pl
根据给定的 gtf 文件和 id.txt 文件，筛选出特定的基因和转录本信息，并输出到新的 GTF 文件中。
```bash
perl gtf_extract.pl  input.gtf id.txt > output.gtf
```
**id.txt**需要一个包含两列的文本文件。第一列是 transcript_id（转录本 ID），第二列是 gene_id（基因 ID），两列之间用制表符（tab）分隔。
