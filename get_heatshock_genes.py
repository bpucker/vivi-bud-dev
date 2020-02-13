annotation_file = "Vv_gene_annotation_file.txt"

candidate_gene_file = "heatshock_genes.txt"

counter = 0
with open( candidate_gene_file, "w" ) as out:
	with open( annotation_file, "r" ) as f:
		f.readline()	#remove header
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			anno = parts[1].lower()
			if "heat" in anno and "shock" in anno:
				out.write( parts[0] + '\n' )
				counter += 1
			line = f.readline()
print "number of detected genes: " + str( counter )
