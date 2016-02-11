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
cassava_genome = Genome_manager('Mesculenta_305_v6.fa','Mesculenta_305_v6.1.gene.gff3')
cassava_genome.getAllPromoterKnownTSS(2000,0)
