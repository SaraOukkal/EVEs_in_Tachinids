import argparse

def process_file(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            columns = line.strip().split('\t')
            scaffold = columns[0]
            qstart = int(columns[2])
            qend = int(columns[3])
            tcov = float(columns[4])
            
            if tcov >= 0.8:
                outfile.write(f'{scaffold}\t{qstart-1}\t{qend}\n')

def main():
    parser = argparse.ArgumentParser(description='Process MMseqs2 output file and generate BED file')
    parser.add_argument('-i', '--input', type=str, help='Input file path')
    parser.add_argument('-o', '--output', type=str, help='Output file path')
    args = parser.parse_args()

    process_file(args.input, args.output)

if __name__ == "__main__":
    main()
