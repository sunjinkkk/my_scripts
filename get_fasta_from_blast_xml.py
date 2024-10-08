import os
import xml.etree.ElementTree as ET
import argparse

def extract_hsp_hseq_to_fasta(input_file, output_file):
    # 从输入文件名中提取样本名称，假设样本名称在前缀部分，如W0123-1
    sample_name = os.path.basename(input_file).split('.')[0]

    # 解析 XML 文件
    tree = ET.parse(input_file)
    root = tree.getroot()

    # 列表存储所有的 Hsp_hseq 序列
    hseq_list = []

    # 遍历所有的 Hsp 标签
    for hsp in root.iter('Hsp'):
        # 提取 Hsp_hseq 标签的内容
        hseq = hsp.find('Hsp_hseq')
        if hseq is not None:
            hseq_list.append(hseq.text)

    # 倒序排列
    hseq_list.reverse()

    # 打开输出文件，保存为 FASTA 格式
    with open(output_file, 'w') as out_f:
        for i, hseq in enumerate(hseq_list):
            # 使用样本名称作为序列ID，序号从1开始
            out_f.write(f">{sample_name}_sequence_{i+1}\n")
            out_f.write(hseq + "\n")

    print(f"从 {input_file} 提取的序列已保存为 {output_file}")

def main():
    # 使用 argparse 处理命令行参数
    parser = argparse.ArgumentParser(description="提取 BLAST XML 文件中的 Hsp_hseq 序列并保存为 FASTA 格式")
    parser.add_argument("-i", "--input", required=True, help="输入的 XML 文件路径")
    parser.add_argument("-o", "--output", required=True, help="输出的 FASTA 文件路径")

    args = parser.parse_args()

    # 调用函数进行处理
    extract_hsp_hseq_to_fasta(args.input, args.output)

if __name__ == "__main__":
    main()
