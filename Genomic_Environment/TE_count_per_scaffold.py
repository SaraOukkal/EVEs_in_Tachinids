import argparse
import os

def count_TE_per_scaffold(input_file, output_file):
    scaffold_TE_count = {}

    with open(input_file, "r") as infile:
        for line in infile:
            scaffold = line.split("\t")[0]
            if scaffold in scaffold_TE_count:
                scaffold_TE_count[scaffold] += 1
            else:
                scaffold_TE_count[scaffold] = 1

    with open(output_file, "w") as outfile:
        outfile.write("Scaffold\tTE_count\n")
        for scaffold, count in scaffold_TE_count.items():
            outfile.write(f"{scaffold}\t{count}\n")

def main():
    parser = argparse.ArgumentParser(description='Count number of TE per scaffold.')
    parser.add_argument('-i', '--input', type=str, help='Input BED file path')
    parser.add_argument('-o', '--output', type=str, help='Output file path')
    args = parser.parse_args()

    count_TE_per_scaffold(args.input, args.output)

if __name__ == "__main__":
    main()

