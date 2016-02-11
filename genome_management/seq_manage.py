import re
import gzip
import codecs
from operator import itemgetter, attrgetter

utf8 = codecs.getreader('UTF-8')


class Fasta_manager(object):
	def __init__(self, fastaFile, show_genome_stat=False):
		self.chromosomeLength = {}
		self.chromosomeSeq = {}
		self.chromosomeStatistics = {}  # Length, GC, AT, N
		sumGC = sumAT = sumN = sumLength = 0
		
		if(fastaFile.find('.gz') > 0):
			filegz = gzip.open(fastaFile, 'rb')
			self.file = utf8(filegz)
		else:
			self.file = open(fastaFile, 'r')
		fasta = self.file.read().split('>')
		fasta = fasta[1:]
		for chromosome in fasta:
			if (chromosome[:50].find(' ') < 0):
				header = chromosome[:chromosome[:50].find('\n')]
			else:
				header = chromosome[:chromosome[:50].find(' ')]
			sequence = chromosome[chromosome.find('\n'):-1].replace('\n', '')
			
			length = len(sequence)
			self.chromosomeSeq[header] = sequence
			self.chromosomeLength[header] = length

			if show_genome_stat: 
				GC = sequence.count('G')+sequence.count('C')+sequence.count('g')+sequence.count('c')
				AT = sequence.count('A')+sequence.count('T')+sequence.count('a')+sequence.count('t')
				N = sequence.count('N')+sequence.count('n')
				sumGC += GC
				sumAT += AT
				sumN += N
				sumLength += length
				self.chromosomeStatistics[header] = [length, GC, AT, N]
				# print(header ,length, GC, AT, N, sep='\t')
		if show_genome_stat:
			print("summary" ,sumLength, sumGC, sumAT, sumN, sep='\t')
			print("GC content = ", float(sumGC) / (sumAT+sumGC))
			print("ATGC count  = ", sumGC + sumAT)
			print("N content  = ", sumN/sumLength)

	def checkChromosome(self, chromosome):
		if(chromosome in self.chromosomeLength):
			return True 
		else:
			print("Not found "+chromosome+" ,please check chromosome again!!!")
			return False

	def checkChromosome(self, chromosome, start, end):
		if(start>end):
			print("Please the start nucleotide should be less")
			return False
		elif(chromosome in self.chromosomeLength):
			if(end<=self.chromosomeLength[chromosome]):
				return True
			else:
				print("Not found "+chromosome+" at "+ str(end) +", please try again. (the first nucleotide is position = 1)")
				return False
		else:
			print("Not found "+chromosome+" ,please check chromosome again!!!")
			return False

	def getGCcontent(self, sequence):
		GC = sequence.count('G') + sequence.count('C') + sequence.count('g') + sequence.count('c')
		AT = sequence.count('A') + sequence.count('T') + sequence.count('a') + sequence.count('t')
		return float(GC) * 100 / (AT + GC)
	def getGC(self, sequence):
		return sequence.count('G') + sequence.count('C') + sequence.count('g') + sequence.count('c')
	def getStatisticSequence(self, sequence):
		GC = sequence.count('G') + sequence.count('C') + sequence.count('g') + sequence.count('c')
		AT = sequence.count('A') + sequence.count('T') + sequence.count('a') + sequence.count('t')
		N = sequence.count('N') + sequence.count('n')
		return [len(sequence), GC, AT, N, float(GC) * 100 / (AT + GC)]
	def getStatisticSeqFromGenome(self, chromosome, start, end, stand):
		seqLength = self.getChromosomeLength(chromosome)
		if (start > 0 and start < seqLength + 1 and end < seqLength + 1):
			if(stand == '+'):
				return self.getStatisticSequence(self.chromosomeSeq[chromosome][start - 1:end])
			else:
				reverse = self.chromosomeSeq[chromosome][start - 1:end]
				reverse = self.complementary(reverse[::-1])
				return self.getStatisticSequence(reverse)
		else:
			print("Out of length in seq please check again")
			print("chromosome", chromosome, "length:", seqLength)
			print("gene position:", start, "to", end, "on", stand, "stand")
			exit()
	def getChromosomeLength(self, chromosome_name):
		return self.chromosomeLength[chromosome_name]
	def getSequence(self, chromosome, start, end, stand):
		if self.checkChromosome(chromosome, start, end):
			seqLength = self.getChromosomeLength(chromosome)
			if (start > 0 and start < seqLength + 1 and end < seqLength + 1):
				if(stand == '+'):
					return self.chromosomeSeq[chromosome][start - 1:end]
				else:
					reverse = self.chromosomeSeq[chromosome][start - 1:end]
					reverse = self.complementary(reverse[::-1])
					return reverse
			else:
				return False
				print("\nOut of chromosome length, please check again.")
				print("Chromosome length:", seqLength)
				print("Error command: getSequence(", chromosome, start, end, stand, ")", sep=', ')
		else:
			return ""

	def complementary(self, seq):
		new = ""
		for base in seq:
			if(base == 'A'):
				new = new + 'T'
			elif(base == 'T'):
				new = new + 'A'
			elif(base == 'G'):
				new = new + 'C'
			elif(base == 'C'):
				new = new + 'G'
			elif(base == 'a'):
				new = new + 't'
			elif(base == 't'):
				new = new + 'a'
			elif(base == 'g'):
				new = new + 'c'
			elif(base == 'c'):
				new = new + 'g'
			else:
				new = new + base
		return new
	def searchSeqInChromosome(self, chromosome_name, pattern):
		pattern = pattern.upper()
		len_pattern = len(pattern)
		index_found = []
		# Search pattern in plus stand
		index = self.chromosomeSeq[chromosome_name].find(pattern)
		while(index > -1):
			index_found.append([index + 1, index + len_pattern, '+'])
			index = self.chromosomeSeq[chromosome_name].find(pattern, index + 1)
		# Search pattern in minus stand
		pattern = self.complementary(pattern)[::-1]
		index = self.chromosomeSeq[chromosome_name].find(pattern)
		while(index > -1):
			index_found.append([index + 1, index + len_pattern, '-'])
			index = self.chromosomeSeq[chromosome_name].find(pattern, index + 1)
		# Return [fistMatch,endMatch,stand]
		return index_found
	def searchSeqInGenome(self, pattern):
		pattern = pattern.upper()
		len_pattern = len(pattern)
		index_found = []
		for chromosome_name, seq in sorted(self.chromosomeSeq.items()):
			# Search pattern in plus stand
			index = seq.find(pattern)
			while(index > -1):
				index_found.append([chromosome_name , index + 1, index + len_pattern, '+'])
				index = seq.find(pattern, index + 1)
			# Search pattern in minus stand
			pattern = self.complementary(pattern)[::-1]
			index = seq.find(pattern)
			while(index > -1):
				index_found.append([chromosome_name, index + 1, index + len_pattern, '-'])
				index = seq.find(pattern, index + 1)	
		return index_found

class Gff_manager(object):
	def __init__(self, file_name):
		self.data = []
		if(file_name.find('.gz') > 0):
			filegz = gzip.open(file_name, 'rb')
			gff_file = utf8(filegz)
		else:
			gff_file = open(file_name, 'r')
		for line in gff_file:
			if(line[0] != '#'):
				line = line.split()
				line[3] = int(line[3])
				line[4] = int(line[4])
				line[8] = line[8].split(';')
				if(line != ''):
					self.data.append(line)
	def getNumgerOfGffLine(self):
		return len(self.data)
	def getTable(self):
		return self.data
	def getTableSpecificType(self, gene_struc_type):
		table = []
		for line in self.data:
			if(line[2] == gene_struc_type):
				table.append(line)
		return table
	def printdata(self,type="five_prime_UTR"):
		countLine = 0
		for line in self.data:
			if(line[2] == 'five_prime_UTR'):
				print(line[0] + "\t" + line[2] + "\t" + str(line[3]) + "\t" + str(line[4]) + "\t" + line[6] + "\t" + line[8][0])
				countLine += 1
	def getTableDataOfGene(self, geneName):
		table = []
		found = False
		breaker = False
		for line in self.data:
			if (line[2] == 'gene'):
				if(line[8][0].find(geneName) > 0):
					found = True
					breaker = True
				elif(breaker == True): 
					found = False
					break
			if(found == True):
				table.append(line)
		return table
	def getTableDataOfGeneAndType(self, geneName, type):
		table = []
		found = False
		breaker = False
		for line in self.data:
			if (line[2] == 'mRNA'):
				if(line[8][4].find(geneName,7) > 0 and line[8][3][8].find('1')):
					found = True
					breaker = True
				elif(breaker == True): 
					found = False
					break
			if(found == True):
				if(line[2] == type):
					table.append(line)
		return sorted(table,key=itemgetter(4,5))
	def getTranscripthave5UTR(self):
		print("gene", "transcript", "label5UTR", "lengthOf5UTR", "strand", "start", "stop", sep='\t')
		for line in self.data:
			if(line[2] == 'gene'):
				geneName = line[8][0][3:]
			elif(line[2] == 'five_prime_UTR' or line[2] == '5-UTR'):
				transcriptName = line[8][0][3:26]
				label5UTR = line[8][0][-1:]
				start5UTR = int(line[3])
				stop5UTR = int(line[4])
				len5UTR = stop5UTR - start5UTR + 1
				stand = line[6]
				print(geneName, transcriptName, label5UTR, len5UTR, stand, start5UTR, stop5UTR, sep='\t')
	def getGeneList(self):
		for line in self.data:
			if(line[2] == 'gene'):
				print(line[8][0][3:])
	def getDataSpecificType(self,gene_component):
		table = []
		for line in self.data:
			if(line[2] == gene_component):
				table.append(line)
		return table
	def getTranscript(self):
		for line in self.data:
			if(line[2] == 'mRNA'):
				print(line[8][0][3:])

class Genome_manager(Fasta_manager, Gff_manager):
	def __init__(self, fastaFile, GffFile):
		self.fastaFile = fastaFile
		Fasta_manager.__init__(self, fastaFile)
		Gff_manager.__init__(self, GffFile)
	def getGCcontentInTranscript(self, type):
		sumGC = 0
		sumAT = 0
		for line in self.data:
			if(line[2] == type):
				# print(line[8][0][3:], line[0], line[3], line[4] , line[6], sep='\t',end = '\t')
				statistic = Fasta_manager.getStatisticSeqFromGenome(self, line[0], line[3], line[4] , line[6])
				# print(statistic[0], statistic[1], statistic[2], sep='\t')
				sumGC += statistic[1]
				sumAT += statistic[2]
		print("Summary GC content in", type, ":", float(sumGC) * 100 / (sumGC + sumAT))
	def selectedTSSProtein(self, upstream, downstream):
		file_write = open("%s_upstream_-%dto+%d.fa" % (self.fastaFile[:-6], upstream, downstream), 'w')
		statistic_of_5_prime_length = []
		geneListSelected = []
		geneCount = 0
		transcriptName = geneName = ''
		five_prime_UTR = []
		three_prime_UTR = []
		CDS = []
		count_five_prime_UTR_selected = 0
		count_five_prime_UTR_total = 0
		count_upstream_out_of_criteria = 0
		count_seq = 0

		for line in self.data:
			if(line[2] == 'gene'):
				geneName = line[8][0][3:]
				geneCount += 1
			elif(line[2] == 'mRNA'):
				count_five_prime = len(five_prime_UTR)
				if(count_five_prime > 0):
					# Gene have five_prime_UTR
					count_five_prime_UTR_selected += 1
					count_five_prime_UTR_total += count_five_prime
					if geneName not in geneListSelected:
						geneListSelected.append(geneName)
					if(five_prime_UTR[0][6] == '+'):
						five_prime_UTR.sort(key=itemgetter (3, 4))
						selected_five_prime = five_prime_UTR[count_five_prime - 1]
					else:
						five_prime_UTR.sort(key=itemgetter (4, 3))
						selected_five_prime = five_prime_UTR[0]
					sequence = Fasta_manager.getSequence(self, selected_five_prime[0], selected_five_prime[3], selected_five_prime[4], selected_five_prime[6])
					statistic_of_5_prime_length.append(len(sequence))
					# print(">", transcriptName, sep="")
					# print(sequence)
					text = self.getPromoterOfGene(upstream, downstream, selected_five_prime)
					if(text == False):
						count_upstream_out_of_criteria += 1
					else:
						file_write.writelines(text)
						count_seq += 1
				else:
					# Gene have not five_prime_UTR
					pass

				transcriptName = line[8][0][3:]
				five_prime_UTR = []
				three_prime_UTR = []
				CDS = []
			elif(line[2] == 'five_prime_UTR' or line[2] == '5-UTR'):
				five_prime_UTR.append(line)
			elif(line[2] == 'tree_prime_UTR' or line[2] == '3-UTR'):
				three_prime_UTR.append(line)
			elif(line[2] == 'CDS'):
				CDS.append(line)
		# lastLine imporve data
		count_five_prime = len(five_prime_UTR)
		if(count_five_prime > 0):
			count_five_prime_UTR_selected += 1
			count_five_prime_UTR_total += count_five_prime
			if geneName not in geneListSelected:
				geneListSelected.append(geneName)
			if(five_prime_UTR[0][6] == '+'):
				five_prime_UTR.sort(key=itemgetter (3, 4))
				selected_five_prime = five_prime_UTR[count_five_prime - 1]
			else:
				five_prime_UTR.sort(key=itemgetter (4, 3))
				selected_five_prime = five_prime_UTR[0]
			sequence = Fasta_manager.getSequence(self, selected_five_prime[0], selected_five_prime[3], selected_five_prime[4], selected_five_prime[6])
			statistic_of_5_prime_length.append(len(sequence))
			# print(">", transcriptName, sep="")
			# print(sequence)
			text = self.getPromoterOfGene(upstream, downstream, selected_five_prime)
			if(text == False):
				count_upstream_out_of_criteria += 1
			else:
				file_write.writelines(text)
				count_seq += 1

		# Get statistic
		print("Statistic of genome", "%s_upstream_-%dto+%d.fa" % (self.fastaFile[:-6], upstream, downstream))
		print("Number of annotated gene:", geneCount)
		print("Number of 5'UTR of known gene:", len(geneListSelected))
		print("Number of alternative 5'UTR transcript:", count_five_prime_UTR_total)
		print("Number of selected 5'UTR transcript (unique):", count_five_prime_UTR_selected)
		print("Upstream correct:", count_seq)
		print("Upstream out of criteria:", count_upstream_out_of_criteria)
		# Number of 5'UTR of selected transcript
	def getPromoterOfGene(self, upstream, downstream, five_prime_UTR):
		if(five_prime_UTR[6] == '+'):
			seq = Fasta_manager.getSequence(self, five_prime_UTR[0], five_prime_UTR[3] - upstream, five_prime_UTR[3] + downstream, five_prime_UTR[6])
		else:
			seq = Fasta_manager.getSequence(self, five_prime_UTR[0], five_prime_UTR[4] - downstream, five_prime_UTR[4] + upstream, five_prime_UTR[6])
			
		if(seq == False):
			return False
		else:
			if(seq.count('N') == 0):
				if(five_prime_UTR[6] == '+'):
					text = ">" + five_prime_UTR[8][0][3:] + "|" + str(five_prime_UTR[3] - upstream) + "|" + str(five_prime_UTR[3] + downstream) + "|+\n"
				else:
					text = ">" + five_prime_UTR[8][0][3:] + "|" + str(five_prime_UTR[4] - downstream) + "|" + str(five_prime_UTR[4] + upstream) + "|-\n"
				
				if(len(seq) != upstream + downstream + 1):
					print("\nLength of sequence not correct please check code it again.")
					exit()
				text = text + str(seq) + "\n"
				 # print(text)
				return text
			else:
				return False		
	def getAllPromoterKnownTSS(self, upstream, downstream):
		# Retrive upstream and downstream sequence from TSS
		not_selected = 0
		not_selected_polyN = 0
		count_seq = 0
		for line in self.data:
			if(line[2] == 'five_prime_UTR'):
				if(line[6] == '+'):
					if(line[3] > upstream):
						seq = Fasta_manager.getSequence(self, line[0], line[3] - upstream, line[3] + downstream, line[6])
					else:
						seq = Fasta_manager.getSequence(self, line[0], 1, line[3] + downstream, line[6])
				else:
					if(line[4]+upstream <= Fasta_manager.getChromosomeLength(self, line[0])):
						seq = Fasta_manager.getSequence(self, line[0], line[4] - downstream, line[4] + upstream, line[6])
					else:
						seq = Fasta_manager.getSequence(self, line[0], line[4] - downstream, Fasta_manager.getChromosomeLength(self, line[0]), line[6])
				if(seq == False):
					not_selected += 1
				else:
					if(seq.count('N') == 0):
						if(len(seq) == upstream + downstream + 1):
							if(line[6] == '+'):
								print(">", line[8][0][3:],"|",line[0],"|",line[3]-upstream,"|",line[3]+downstream,"|+" ,sep='')
							else:
								print(">", line[8][0][3:],"|",line[0],"|",line[4]-downstream,"|",line[4]+upstream, "|-" ,sep='')
							print(seq)
							count_seq += 1
						else:
							not_selected +=1
					else:
						not_selected_polyN += 1
		print("not selected sequence:", not_selected)
		print("not selected sequence because N:", not_selected_polyN)
		print("It including ", count_seq, "sequences for next step")
