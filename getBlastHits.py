#!/usr/bin/env python
import sys, os, argparse


def get_blast_hits(genes, samples, hits_files, fastq_input_filename, work_path, direction):
	id_dict = {}

	for hits_file, gene_name, sample_name in zip(hits_files, genes, samples):
		with open(hits_file) as hits_file_handle:
			for raw_fastq_id in hits_file_handle:
				fastq_id = raw_fastq_id.strip()
				if fastq_id not in id_dict:
					id_dict[fastq_id] = set()
				#need to remember naming convention
				hits_fastq_file = "{work_path}/fr_readHits/{gene_name}.{sample_name}.{direction}.fastq".format(**locals())
				if os.path.isfile(hits_fastq_file):
					os.remove(hits_fastq_file)
				
				id_dict[fastq_id].add(hits_fastq_file)

	inseq = None
	out_fastq_dict = {}
	with open(fastq_input_filename) as fastq_handle:
		for num, line in enumerate(fastq_handle):
			if num % 4 == 0:
				inseq = False
				if not line.startswith("@"):
					print("ERROR!! Malformed File")
				header = line[1:]
				fqid = line.split(" ")[0][1:] #split at space and exclude first @ symbol
				if fqid in id_dict:
					inseq = True
					export_file_names = id_dict[fqid]
					for file_name in export_file_names:
						if file_name not in out_fastq_dict:
							out_fastq_dict[file_name] = ""
						out_fastq_dict[file_name] += "@{header}".format(**locals())
						#
						#with open(file_name, "a") as out_file_handle: #REMOVE W+ LATER IT WILL BE A HUGE ISSUE
						#	out_file_handle.write("@{header}".format(**locals()))

			elif inseq == True:
				for file_name in export_file_names:
					out_fastq_dict[file_name] += "{line}".format(**locals())
					#with open(file_name, "a") as out_file_handle: #REMOVE W+ LATER IT WILL BE A HUGE ISSUE
					#	out_file_handle.write("{line}".format(**locals()))
			#print (line)
	
	for file_name in out_fastq_dict:
		with open(file_name, "a") as file_out_handle:
			file_out_handle.write(out_fastq_dict[file_name])

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='')

	parser.add_argument('-fastq_file', action="store", dest="fastq_file", required=True) 
	parser.add_argument('-direction', action="store", dest="direction", required=True) 
	parser.add_argument('-work_path', action="store", dest="work_path", required=True) 
	parser.add_argument('-hit_files', action="store", dest="hit_files", nargs = "+", required=True) 
	parser.add_argument('-samples', action="store", dest="samples", nargs = "+", required=True) 
	parser.add_argument('-genes', action="store", dest="genes", nargs = "+", required=True) 

	args = parser.parse_args()

	fastq_input_filename = args.fastq_file
	direction = args.direction
	hits_files = args.hit_files
	genes = args.genes
	samples = args.samples
	work_path = args.work_path

	get_blast_hits(genes, samples, hit_files, fastq_file, work_path, direction)