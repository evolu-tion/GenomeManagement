#!/usr/bin/python

"""
Copyright (c) 2016 King Mongkut's University technology Thonburi
Author: Nattawet Sriwichai
Contact: nattawet.sri@mail.kmutt.ac.th
Version: 2.0.0
License: MIT License
"""

import gzip
import codecs
from operator import itemgetter

VERSION = "2.0.0"
UTF8_READER = codecs.getreader('UTF-8')


class FastaManager:
    def __init__(self, fasta_file, show_genome_stat=False):
        self.chromosome_length = {}
        self.chromosome_seq = {}
        self.chromosome_statistics = {}  # Length, GC, AT, N
        count_bases = {
            'A': 0, 'T': 0, 'C': 0, 'G': 0,
            'R': 0, 'S': 0, 'W': 0, 'K': 0,
            'M': 0, 'B': 0, 'D': 0, 'H': 0,
            'V': 0, 'Y': 0, 'N': 0
        }

        if fasta_file.endswith('.gz'):
            file_gz = gzip.open(fasta_file, 'rb')
            self.file = UTF8_READER(file_gz)
        else:
            self.file = open(fasta_file, 'r')
        
        fasta_content = self.file.read().split('>')
        fasta_content = fasta_content[1:]
        
        for chromosome in fasta_content:
            if chromosome[:50].find(' ') < 0:
                header = chromosome[:chromosome[:50].find('\n')]
            else:
                header = chromosome[:chromosome[:50].find(' ')]
            
            sequence = chromosome[chromosome.find('\n'):].replace('\n', '')
            length = len(sequence)
            
            self.chromosome_seq[header] = sequence
            self.chromosome_length[header] = length

            if show_genome_stat:
                for base in sequence:
                    count_bases[base.upper()] += 1

        if show_genome_stat:
            total_bases = sum(count_bases.values())
            print("Total sequence length:", "{:0,}".format(total_bases))
            print("Total ungapped length:", "{:0,}".format(total_bases - count_bases['N']))
            print("Total spanned gaps:", "{:0,}".format(count_bases['N']))
            print("Number of chromosomes/scaffolds/contigs: ", "{:0,}".format(len(fasta_content)))
            
            sum_gc = (count_bases['G'] + count_bases['C'] + count_bases['S'] + 
                      count_bases['Y']/2 + count_bases['K']/2 + count_bases['M']/2 + 
                      count_bases['B']*2/3 + count_bases['D']/3 + count_bases['H']/3 + 
                      count_bases['V']*2/3 + count_bases['N']/2)
            
            if total_bases > 0:
                print("GC content (%):", "{:0,.2f}".format(sum_gc * 100 / total_bases))
                print("N content (%):", "{:0,.2f}".format(count_bases['N'] * 100 / total_bases))
            
            scaffold_len = sorted(self.chromosome_length.values(), reverse=True)
            half_sum_len = sum(scaffold_len) / 2

            current_sum_len = 0
            i = 0
            while i < len(scaffold_len) and current_sum_len < half_sum_len:
                current_sum_len += scaffold_len[i]
                i += 1

            if i > 0:
                print("N50:", "{:0,}".format(scaffold_len[i-1]))
                print("L50:", "{:0,}".format(i))

    def check_chromosome(self, chromosome, start=0, end=1):
        if start > end:
            print(f"Error: check_chromosome({chromosome}, {start}, {end}), start > end")
            return False
        
        if chromosome in self.chromosome_length:
            if end <= self.chromosome_length[chromosome]:
                return True
            else:
                print(f"Not found {chromosome} at {end}, please try again.")
                return False
        else:
            print(f"Not found {chromosome}, please check chromosome name.")
            return False

    def get_gc_content(self, sequence):
        gc = sequence.lower().count('g') + sequence.lower().count('c')
        at = sequence.lower().count('a') + sequence.lower().count('t')
        total = at + gc
        return float(gc) * 100 / total if total > 0 else 0.0

    def get_gc(self, sequence):
        return sequence.lower().count('g') + sequence.lower().count('c')

    def get_statistic_sequence(self, sequence):
        gc = sequence.lower().count('g') + sequence.lower().count('c')
        at = sequence.lower().count('a') + sequence.lower().count('t')
        n = sequence.lower().count('n')
        total = at + gc
        gc_percent = float(gc) * 100 / total if total > 0 else 0.0
        return [len(sequence), gc, at, n, gc_percent]

    def get_statistic_seq_from_genome(self, chromosome, start, end, strand):
        seq_length = self.get_chromosome_length(chromosome)
        if 0 < start < seq_length + 1 and end < seq_length + 1:
            if strand == '+':
                return self.get_statistic_sequence(self.chromosome_seq[chromosome][start - 1:end])
            else:
                reverse = self.chromosome_seq[chromosome][start - 1:end]
                reverse = self.complementary(reverse[::-1])
                return self.get_statistic_sequence(reverse)
        else:
            print("Out of length in seq please check again")
            print(f"chromosome {chromosome} length: {seq_length}")
            print(f"gene position: {start} to {end} on {strand} strand")
            exit(1)

    def get_chromosome_length(self, chromosome_name):
        return self.chromosome_length.get(chromosome_name, 0)

    def get_sequence(self, chromosome, start, end, strand):
        if self.check_chromosome(chromosome, start, end):
            if strand == '+':
                return self.chromosome_seq[chromosome][start - 1:end]
            else:
                reverse = self.chromosome_seq[chromosome][start - 1:end]
                return self.complementary(reverse[::-1])
        return ""

    def get_chr_sequence(self, chromosome):
        return self.chromosome_seq.get(chromosome, "")

    def complementary(self, seq):
        mapping = {
            'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
            'a': 't', 't': 'a', 'g': 'c', 'c': 'g'
        }
        return "".join([mapping.get(base, base) for base in seq])

    def search_seq_in_chromosome(self, chromosome_name, pattern):
        pattern = pattern.upper()
        len_pattern = len(pattern)
        index_found = []
        
        # Search plus strand
        seq = self.chromosome_seq[chromosome_name]
        index = seq.find(pattern)
        while index > -1:
            index_found.append([index + 1, index + len_pattern, '+'])
            index = seq.find(pattern, index + 1)
            
        # Search minus strand
        pattern_rev = self.complementary(pattern)[::-1]
        index = seq.find(pattern_rev)
        while index > -1:
            index_found.append([index + 1, index + len_pattern, '-'])
            index = seq.find(pattern_rev, index + 1)
            
        return index_found

    def search_seq_in_genome(self, pattern):
        pattern = pattern.upper()
        len_pattern = len(pattern)
        index_found = []
        
        pattern_rev = self.complementary(pattern)[::-1]
        
        for chromosome_name, seq in sorted(self.chromosome_seq.items()):
            # Plus strand
            index = seq.find(pattern)
            while index > -1:
                index_found.append([chromosome_name, index + 1, index + len_pattern, '+'])
                index = seq.find(pattern, index + 1)
            
            # Minus strand
            index = seq.find(pattern_rev)
            while index > -1:
                index_found.append([chromosome_name, index + 1, index + len_pattern, '-'])
                index = seq.find(pattern_rev, index + 1)
        return index_found


class GffManager:
    def __init__(self, file_name):
        self.data = []  # all line
        self.gene_struc = {}  # {'gene1': [line1 ,line2, line 3]}
        self.chromosome_contain_gene = {}  # {[(chr1,'+')]: [..........]}

        gene_name = ""
        if file_name.endswith('.gz'):
            file_gz = gzip.open(file_name, 'rb')
            gff_file = UTF8_READER(file_gz)
        else:
            gff_file = open(file_name, 'r')
            
        gene_annotation = []
        
        for line in gff_file:
            if not line.startswith('#') and line.strip():
                parts = line.split()
                # Assuming standard GFF3/GTF columns
                # col 3 and 4 are start/end (1-based)
                # But in original code it accessed index 3 and 4 after split.
                # GFF: seqid source type start end score strand phase attributes
                # indices: 0     1      2    3     4   5     6      7     8
                
                parts[3] = int(parts[3])
                parts[4] = int(parts[4])
                parts[8] = parts[8].split(';')
                self.data.append(parts)
                
                if parts[2] != 'gene':
                    gene_annotation.append(parts)
                else:
                    if gene_name:
                        self.gene_struc[gene_name] = gene_annotation
                    gene_annotation = [parts]
                    # Attempt to extract gene name from attributes
                    # Original code: line[8][1][5:] -> assuming specific format "ID=...;Name=..."
                    try:
                        gene_name = parts[8][1][5:] 
                    except IndexError:
                        gene_name = "Unknown"

        if gene_name:
            self.gene_struc[gene_name] = gene_annotation

        # Build chromosome contain gene index
        table = self.get_table_specific_type("gene")
        table = sorted(table, key=itemgetter(0, 6, 3, 4))
        
        for line in table:
            key = (line[0], line[6])
            if key not in self.chromosome_contain_gene:
                self.chromosome_contain_gene[key] = []
            self.chromosome_contain_gene[key].append(line)
            
        for key, value in self.chromosome_contain_gene.items():
            if key[1] == '-':
                self.chromosome_contain_gene[key] = sorted(value, key=itemgetter(4, 3), reverse=True)

    def get_number_of_gff_lines(self):
        return len(self.data)

    def get_table(self):
        return self.data

    def get_table_specific_type(self, gene_struc_type):
        return [line for line in self.data if line[2] == gene_struc_type]

    def get_table_specific_type_and_strand(self, gene_struc_type, strand):
        return [line for line in self.data if line[2] == gene_struc_type and line[6] == strand]

    def print_data(self, feature_type="five_prime_UTR"):
        count_line = 0
        for line in self.data:
            if line[2] == feature_type:
                print(f"{line[0]}\t{line[2]}\t{line[3]}\t{line[4]}\t{line[6]}\t{line[8][0]}")
                count_line += 1

    def get_table_data_of_gene_and_type(self, gene_name, feature_type):
        if gene_name not in self.gene_struc:
            return []
        return [i for i in self.gene_struc[gene_name] if i[2] == feature_type]

    def get_table_data_of_gene(self, gene_name):
        return self.gene_struc.get(gene_name, [])

    def get_transcript_have_5utr(self):
        print("gene", "transcript", "label5UTR", "lengthOf5UTR", "strand", "start", "stop", sep='\t')
        gene_name = ""
        for line in self.data:
            if line[2] == 'gene':
                # Parsing logic from original code
                # gene_name = line[8][0][3:] -> ID=gene:...
                gene_name = line[8][0][3:]
            elif line[2] in ('five_prime_UTR', '5-UTR'):
                transcript_name = line[8][0][3:26]
                label_5utr = line[8][0][-1:]
                start_5utr = int(line[3])
                stop_5utr = int(line[4])
                len_5utr = stop_5utr - start_5utr + 1
                strand = line[6]
                print(gene_name, transcript_name, label_5utr, len_5utr, strand, start_5utr, stop_5utr, sep='\t')

    def get_gene_list(self):
        return sorted(list(self.gene_struc.keys()))

    def get_data_specific_type(self, gene_component):
        return [line for line in self.data if line[2] == gene_component]

    def get_transcript(self):
        for line in self.data:
            if line[2] == 'mRNA':
                print(line[8][0][3:])

    def check_gene(self, gene_name):
        return gene_name in self.gene_struc

    def get_gene_forward(self, gene_name):
        # return end position of forward gene, if don't have forward gene return False
        if gene_name not in self.gene_struc:
            return False
            
        x = self.gene_struc[gene_name]
        chromosome = x[0][0]
        strand = x[0][6]
        start = x[0][3]
        end = x[0][4]
        
        table_gene = self.chromosome_contain_gene.get((chromosome, strand), [])
        
        if strand == "+":
            i = 0
            while i < len(table_gene) and table_gene[i][3] < start:
                i += 1
            i -= 1
            if i == -1:
                return False
            return table_gene[i][4]
        else:
            i = 0
            while i < len(table_gene) and end < table_gene[i][4]:
                i += 1
            i -= 1
            if i == -1:
                return False
            return table_gene[i][3]


class GenomeManager(FastaManager, GffManager):
    def __init__(self, fasta_file, gff_file):
        self.fasta_file = fasta_file
        FastaManager.__init__(self, fasta_file)
        GffManager.__init__(self, gff_file)
        self.list_of_gene_no_promoter = []

    def get_list_of_gene_no_promoter(self):
        return self.list_of_gene_no_promoter

    def get_gc_content_in_transcript(self, feature_type):
        sum_gc = 0
        sum_at = 0
        for line in self.data:
            if line[2] == feature_type:
                statistic = self.get_statistic_seq_from_genome(line[0], line[3], line[4], line[6])
                sum_gc += statistic[1]
                sum_at += statistic[2]
        
        total = sum_gc + sum_at
        if total > 0:
            print(f"Summary GC content in {feature_type} : {float(sum_gc) * 100 / total}")

    def selected_tss_protein(self, upstream, downstream):
        import re
        file_write = open(f"{self.fasta_file[:-6]}_upstream_-{upstream}to+{downstream}.fa", 'w')
        gene_list_selected = []
        gene_count = 0
        five_prime_utr = []
        count_five_prime_utr_selected = 0
        count_five_prime_utr_total = 0
        count_upstream_out_of_criteria = 0
        count_seq = 0
        gene_name = ""

        # Processing logic preserved from original but standardized
        for line in self.data:
            if line[2] == 'gene':
                gene_name = line[8][0][3:]
                gene_count += 1
            elif line[2] == 'mRNA':
                count_five_prime = len(five_prime_utr)
                if count_five_prime > 0:
                    count_five_prime_utr_selected += 1
                    count_five_prime_utr_total += count_five_prime
                    if gene_name not in gene_list_selected:
                        gene_list_selected.append(gene_name)
                    
                    if five_prime_utr[0][6] == '+':
                        five_prime_utr.sort(key=itemgetter(3, 4))
                        selected_five_prime = five_prime_utr[count_five_prime - 1]
                    else:
                        five_prime_utr.sort(key=itemgetter(4, 3))
                        selected_five_prime = five_prime_utr[0]
                    
                    text = self.get_promoter_of_gene(upstream, downstream, selected_five_prime)
                    if text is False:
                        count_upstream_out_of_criteria += 1
                    else:
                        file_write.writelines(text)
                        count_seq += 1
                
                five_prime_utr = []
            elif line[2] in ('five_prime_UTR', '5-UTR'):
                five_prime_utr.append(line)
        
        # Last group
        if five_prime_utr:
            count_five_prime = len(five_prime_utr)
            count_five_prime_utr_selected += 1
            count_five_prime_utr_total += count_five_prime
            if gene_name not in gene_list_selected:
                gene_list_selected.append(gene_name)
            
            if five_prime_utr[0][6] == '+':
                five_prime_utr.sort(key=itemgetter(3, 4))
                selected_five_prime = five_prime_utr[count_five_prime - 1]
            else:
                five_prime_utr.sort(key=itemgetter(4, 3))
                selected_five_prime = five_prime_utr[0]
                
            text = self.get_promoter_of_gene(upstream, downstream, selected_five_prime)
            if text is False:
                count_upstream_out_of_criteria += 1
            else:
                file_write.writelines(text)
                count_seq += 1

        file_write.close()
        print(f"Statistic of genome {self.fasta_file[:-6]}_upstream_-{upstream}to+{downstream}.fa")
        print("Number of annotated gene:", gene_count)
        print("Number of 5'UTR of known gene:", len(gene_list_selected))
        print("Number of alternative 5'UTR transcript:", count_five_prime_utr_total)
        print("Number of selected 5'UTR transcript (unique):", count_five_prime_utr_selected)
        print("Upstream correct:", count_seq)
        print("Upstream out of criteria:", count_upstream_out_of_criteria)

    def get_promoter_of_gene(self, upstream, downstream, five_prime_utr):
        chromosome = five_prime_utr[0]
        start = five_prime_utr[3]
        end = five_prime_utr[4]
        strand = five_prime_utr[6]
        
        if strand == '+':
            seq = self.get_sequence(chromosome, start - upstream, start + downstream, strand)
        else:
            seq = self.get_sequence(chromosome, end - downstream, end + upstream, strand)
            
        if seq is False:
            return False
        
        if seq.count('N') == 0:
            name = five_prime_utr[8][0][3:]
            if strand == '+':
                header = f">{name}|{start - upstream}|{start + downstream}|+\n"
            else:
                header = f">{name}|{end - downstream}|{end + upstream}|-\n"
            
            if len(seq) != upstream + downstream + 1:
                print("\nLength of sequence not correct please check code it again.")
                exit(1)
            return header + str(seq) + "\n"
        else:
            return False

    def get_all_promoter_known_tss(self, upstream, downstream):
        not_selected = 0
        not_selected_poly_n = 0
        count_seq = 0
        for line in self.data:
            if line[2] == 'five_prime_UTR':
                chromosome = line[0]
                start = line[3]
                end = line[4]
                strand = line[6]
                name = line[8][0][3:]
                
                if strand == '+':
                    if start > upstream:
                        seq = self.get_sequence(chromosome, start - upstream, start + downstream, strand)
                    else:
                        seq = self.get_sequence(chromosome, 1, start + downstream, strand)
                else:
                    chr_len = self.get_chromosome_length(chromosome)
                    if end + upstream <= chr_len:
                        seq = self.get_sequence(chromosome, end - downstream, end + upstream, strand)
                    else:
                        seq = self.get_sequence(chromosome, end - downstream, chr_len, strand)
                
                if seq is False:
                    not_selected += 1
                else:
                    if seq.count('N') == 0:
                        if len(seq) == upstream + downstream + 1:
                            if strand == '+':
                                print(f">{name}|{chromosome}|{start-upstream}|{start+downstream}|+")
                            else:
                                print(f">{name}|{chromosome}|{end-downstream}|{end+upstream}|-")
                            print(seq)
                            count_seq += 1
                        else:
                            not_selected += 1
                    else:
                        not_selected_poly_n += 1
        print("not selected sequence:", not_selected)
        print("not selected sequence because N:", not_selected_poly_n)
        print("It including ", count_seq, "sequences for next step")

    def check_correct_position(self, chromosome_name, gene_name, prom_start, prom_end, strand, min_len, removed_n_gap):
        import re
        chromosome_len = self.get_chromosome_length(chromosome_name)
        if prom_start < 1:
            prom_start = 1
        elif prom_start > chromosome_len:
            prom_start = chromosome_len
            prom_end = chromosome_len
        elif prom_end > chromosome_len:
            prom_end = chromosome_len
        if prom_end < 1:
            prom_end = 1

        forward_end_pos = self.get_gene_forward(gene_name)
        if forward_end_pos is not False:
            if strand == '+':
                if prom_start < forward_end_pos + 1 and prom_end > forward_end_pos:
                    prom_start = forward_end_pos + 1
                elif prom_end == forward_end_pos:
                    prom_start = prom_end
                elif prom_end < forward_end_pos:
                    return False
            elif strand == '-':
                if prom_end > forward_end_pos - 1 and prom_start < forward_end_pos:
                    prom_end = forward_end_pos - 1
                elif prom_start == forward_end_pos:
                    prom_end = prom_start
                elif prom_start > forward_end_pos:
                    return False

        if removed_n_gap:
            sequence = self.get_sequence(chromosome_name, prom_start, prom_end, strand)
            if not sequence:
                return False
                
            if sequence[0] == 'N':
                pos = re.search('[ATGCatgc]+', sequence)
                if pos:
                    seq = sequence[pos.start():]
                    if strand == '+' and len(seq) >= min_len:
                        prom_start += pos.start()
                        return {'promoter_start': prom_start, 'promoter_end': prom_end}
                    if strand == '-' and len(seq) >= min_len:
                        prom_end -= pos.start()
                        return {'promoter_start': prom_start, 'promoter_end': prom_end}
                    return False
                return False
            elif prom_end - prom_start + 1 >= min_len:
                return {'promoter_start': prom_start, 'promoter_end': prom_end}
            else:
                return False
        else:
            seq = self.get_sequence(chromosome_name, prom_start, prom_end, strand)
            if seq and len(seq) >= min_len:
                return {'promoter_start': prom_start, 'promoter_end': prom_end}
            return False

    def get_promoter_of_gene_from_tls(self, gene_name, upstream, downstream, promoter_min_len, removed_n_gap, output_format):
        gene_struc_table = self.get_table_data_of_gene_and_type(gene_name, "CDS")
        if not gene_struc_table:
            # print("Gene name is not correct or has no CDS") 
            # Original code exited here, but it might be better to return empty
            # But the original code was: exit()
            print("Gene name is not correct, please check it again")
            exit(1)
        
        strand = gene_struc_table[0][6]
        chromosome = gene_struc_table[0][0]
        
        if strand == '+':
            promoter_start = gene_struc_table[0][3] - upstream
            promoter_end = gene_struc_table[0][3] + downstream - 1
        else:
            gene_struc_table = sorted(gene_struc_table, key=itemgetter(4), reverse=True)
            promoter_start = gene_struc_table[0][4] - downstream + 1
            promoter_end = gene_struc_table[0][4] + upstream

        new_promoter_position = self.check_correct_position(chromosome, gene_name, promoter_start, promoter_end, strand, promoter_min_len, removed_n_gap)
        
        if new_promoter_position is False:
            self.list_of_gene_no_promoter.append(gene_name)
            return ''
        else:
            promoter_start = new_promoter_position['promoter_start']
            promoter_end = new_promoter_position['promoter_end']
            seq = self.get_sequence(chromosome, promoter_start, promoter_end, strand)
            
            text = ""
            if output_format.lower() in ('fasta', 'fa'):
                text = f">{gene_name}_promoter|{chromosome}|{promoter_start}|{promoter_end}|{strand}|length={promoter_end - promoter_start + 1}|Promoter from CDS|{VERSION}\n{seq}\n"
            elif output_format.lower() in ('gff', 'gff3'):
                text = f"{chromosome}\t{VERSION}\tpromoter\t{promoter_start}\t{promoter_end}\t.\t{strand}\t.\tID={gene_name}_promoter;Name={gene_name};length={promoter_end-promoter_start+1}\n"
            return text

    def get_all_promoter_of_gene_from_tls(self, upstream, downstream, promoter_min_len, removed_n_gap, output_format):
        for gene_name in self.get_gene_list():
            self.get_promoter_of_gene_from_tls(gene_name, upstream, downstream, promoter_min_len, removed_n_gap, output_format)
        print("\n-----------List of gene no promoter-----------")
        for gene_name in self.list_of_gene_no_promoter:
            print(gene_name)

    def get_promoter_of_gene_from_tss(self, gene_name, upstream, downstream, promoter_min_len, removed_n_gap, output_format):
        gene_struc_table = self.get_table_data_of_gene_and_type(gene_name, "five_prime_UTR")
        if not gene_struc_table:
            # Fallback to TLS
            return self.get_promoter_of_gene_from_tls(gene_name, upstream, downstream, promoter_min_len, removed_n_gap, output_format)
        
        strand = gene_struc_table[0][6]
        chromosome = gene_struc_table[0][0]
        
        if strand == '+':
            gene_struc_table = sorted(gene_struc_table, key=itemgetter(3), reverse=True)
            promoter_start = gene_struc_table[0][3] - upstream
            promoter_end = gene_struc_table[0][3] + downstream - 1
        else:
            gene_struc_table = sorted(gene_struc_table, key=itemgetter(4))
            promoter_start = gene_struc_table[0][4] - downstream + 1
            promoter_end = gene_struc_table[0][4] + upstream

        new_promoter_position = self.check_correct_position(chromosome, gene_name, promoter_start, promoter_end, strand, promoter_min_len, removed_n_gap)
        
        if new_promoter_position is False:
            self.list_of_gene_no_promoter.append(gene_name)
            return ''
        else:
            promoter_start = new_promoter_position['promoter_start']
            promoter_end = new_promoter_position['promoter_end']
            seq = self.get_sequence(chromosome, promoter_start, promoter_end, strand)

            text = ""
            if output_format.lower() in ('fasta', 'fa'):
                text = f">{gene_name}_promoter|{chromosome}|{promoter_start}|{promoter_end}|{strand}|length={promoter_end - promoter_start + 1}|Promoter from 5'UTR|{VERSION}\n{seq}\n"
            elif output_format.lower() in ('gff', 'gff3'):
                text = f"{chromosome}\t{VERSION}\tpromoter\t{promoter_start}\t{promoter_end}\t.\t{strand}\t.\tID={gene_name}_promoter;Name={gene_name};length={promoter_end-promoter_start+1}\n"
            return text

    def get_all_promoter_of_gene_from_tss(self, upstream, downstream, promoter_min_len, removed_n_gap, output_format):
        for gene_name in self.get_gene_list():
            self.get_promoter_of_gene_from_tss(gene_name, upstream, downstream, promoter_min_len, removed_n_gap, output_format)
            
        print("\n-----------List of gene no promoter-----------")
        for gene_name in self.list_of_gene_no_promoter:
            print(gene_name)
