from seq_manage import Fasta_manager
from seq_manage import Gff_manager
from seq_manage import Genome_manager

# ############################## MAIN ##################################
# # read genome sequence
# cassava = Fasta_manager('Mesculenta_305_v6.fa')

# # exmaple retrieving sequence use: cassava.getSequence(scaffold, start, end, strand)
# seq = cassava.getSequence('Chromosome01', 100, 200, '+')
# print(seq)


# gene_gff = Gff_manager('Mesculenta_305_v6.1.gene.gff3')
# cassava_genome = Genome_manager('Mesculenta_305_v6.fa','Mesculenta_305_v6.1.gene.gff3')
# cassava_genome.getAllPromoterKnownTSS(2000,0)


gff_file_name = 'Mesculenta_305_v6.1.gene.gff3'
genome_seq_file_name = 'Mesculenta_305_v6.fa'

# gff = Gff_manager(gff_file_name)
# # table = gff.getTableDataOfGene("Manes.01G001000")
# table = gff.getTableDataOfGeneAndType("Manes.01G001400","CDS")
# for i in table:
# 	print(i)

genome = Genome_manager(genome_seq_file_name,gff_file_name)
upstream = 2000
downstream = 0
promoter_minimum_length = 500
gene_exmple = 'Manes.01G000200'

# genome.getPromoterOfGeneFromTLS(gene_exmple, upstream,downstream, promoter_minimum_length)

# genome.getAllPromoterOfGeneFromTLS(upstream, downstream, promoter_minimum_length)
# genome.getPromoterOfGeneFromTSS(gene_exmple, upstream,downstream, promoter_minimum_length)
genome.getAllPromoterOfGeneFromTSS(upstream, downstream, promoter_minimum_length)
