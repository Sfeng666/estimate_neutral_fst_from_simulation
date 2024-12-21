import optparse

def count_alleles(seq):
    count_0 = seq.count('0')
    count_1 = seq.count('1')
    return count_0, count_1

def process_ms_file(sample_size, ms_result_file, count_table_pop1, count_table_pop2, win_info):
    # Read the ms result file and process each line
    with open(ms_result_file, 'r') as infile, open(count_table_pop1, 'w') as out_ct1, open(count_table_pop2, 'w') as out_ct2, open(win_info, 'w') as out_win:
        lines = infile.readlines()

        pop1_sequences = []
        pop2_sequences = []
        positions = False

        for line in lines:
            line = line.strip()
            # first add genotype of all samples to the corresponding list
            if line.startswith("//"):
                pop1_sequences = []
                pop2_sequences = []
                positions = False
            elif line.startswith("segsites:"):
                num_segsites = int(line.split(":")[1].strip())  # number of segregating sites within each simulated window
                out_win.write(f"{num_segsites}\n")
            elif line.startswith("positions:"):
                positions = True
            elif positions and line:
                if len(pop1_sequences) < sample_size:   # assuming the first sample_size number of samples in the ms result are from population 1
                    pop1_sequences.append(line)
                else:   # the rest of the samples are from population 2
                    pop2_sequences.append(line)

                # then count the number of ancestral and derived alleles
                if len(pop1_sequences) == sample_size and len(pop2_sequences) == sample_size:
                    for i in range(num_segsites):
                        pop1_column = ''.join([seq[i] for seq in pop1_sequences])
                        pop2_column = ''.join([seq[i] for seq in pop2_sequences])
                        count_0_pop1, count_1_pop1 = count_alleles(pop1_column)
                        count_0_pop2, count_1_pop2 = count_alleles(pop2_column)
                        out_ct1.write(f"{count_0_pop1}\t{count_1_pop1}\n")
                        out_ct2.write(f"{count_0_pop2}\t{count_1_pop2}\n")
                    pop1_sequences = []
                    pop2_sequences = []
                    positions = False

        print(f"Count tables generated: {count_table_pop1}, {count_table_pop2}")

def main():
    parser = optparse.OptionParser()
    parser.add_option('--sample_size', dest='sample_size', type='int', help='Input parameter of sample size for each population', metavar = "INT")
    parser.add_option('--in_ms_result', dest='ms_result_file', help='Input ms result file', metavar = "PATH")
    parser.add_option('--out_ct_pop1', dest='count_table_pop1', help='Output count table file for population 1 (.ct)', metavar = "PATH")
    parser.add_option('--out_ct_pop2', dest='count_table_pop2', help='Output count table file for population 2 (.ct)', metavar = "PATH")
    parser.add_option('--out_win', dest='win_info', help='Output table of the number of segregating sites across populations within each window (.win)', metavar = "PATH")

    (options, args) = parser.parse_args()

    if not all([options.sample_size, options.ms_result_file, options.count_table_pop1, options.count_table_pop2, options.win_info]):
        parser.error('Missing required arguments')

    process_ms_file(options.sample_size, options.ms_result_file, options.count_table_pop1, options.count_table_pop2, options.win_info)

if __name__ == '__main__':
    main()
