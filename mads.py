### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.1 ###

import re, sys, os
from pandas import DataFrame
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
sys.path.insert(0, '/vol/python/lib64/python2.7')
sys.path.insert(0, '/vol/python/lib64/python2.7/site-packages')

# --- end of imports --- #

def load_expression( filename ):
	"""! @brief load FPKMs from given file """
	
	data = {}
	with open( filename, "r" ) as f:
		headers = f.readline().strip().split('\t')[1:]
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			exp = {}
			for idx, value in enumerate( parts[1:] ):
				header = headers[ idx ][:6]
				try:
					exp[ header ].append( float( value ) )
				except KeyError:
					exp.update( { header: [ float( value ) ] } )
			mean_exp = {}
			for key in exp.keys():
				mean_exp.update( { key: np.mean( exp[ key ] ) } )
			data.update( { parts[0]: mean_exp } )
			line = f.readline()
	return data


def construct_data_output_file( data, candidate_genes, candidate_name_mapping_table, candidate_samples ):
	"""! @brief write expression values of all candidate genes into output file """
	
	datamatrix = []
	genes = []
	for gene in candidate_genes:
		expression = [ ]
		for tissue in candidate_samples:
			expression.append( data[ gene ][ tissue ] )
		datamatrix.append( expression )
		genes.append( candidate_name_mapping_table[ gene ] )
	return genes, candidate_samples, datamatrix


def construct_heatmap( datamatrix, genes, tissues, heatmap_file ):
	"""! @brief construct heatmap from given data matrix """
	
	my_vmax = 100
	
	print "number of genes for heatmap construction: " + str( len( genes ) )
	df = DataFrame( datamatrix, index=genes[::-1], columns=tissues).round( 0 )
	
	fig, ax = plt.subplots(  )
	
	sns.heatmap( df, vmin=0, vmax= my_vmax, ax=ax, linewidths=0.3, annot=True, annot_kws={'fontsize':3}, cbar=False, fmt="g", cmap='YlGnBu' )	#binary	#cmap='YlGnBu'  = 1
	
	for idx, gene in enumerate( genes ):
		ax.text( -3.7, idx+0.6, gene, fontsize=4 )
	
	for idx, tissue in enumerate( tissues ):
		ax.text( idx+0.4, len( genes )+0.5, "20"+tissue[-2:] + "-" + tissue[2:4] + "-" + tissue[:2], rotation=90, fontsize=4 )	#NEW
	
	ax.set_yticklabels( [], rotation=0, fontsize=2 )
	ax.set_xticklabels( [] , rotation=90, fontsize=3  )
	
	ax.spines['bottom'].set_visible(False)
	ax.spines['left'].set_visible(False)
	
	ax.axes.get_yaxis().set_visible(False)
	ax.axes.get_xaxis().set_visible(False)
	
	plt.yticks( rotation=0 )
	plt.subplots_adjust( left=0.17, right=0.99, top=0.99, bottom=0.08, wspace=0.2 )
	
	plt.savefig( heatmap_file, dpi=900  )
	plt.savefig( heatmap_file.replace( ".pdf", ".svg" ) )


def filter_genes( exp_data, genes, exp_cutoff ):
	"""! @brief filter genes based on expression """
	
	final_genes = []
	for gene in genes:
		exp = exp_data[ gene ].values()
		if np.mean( exp ) >= exp_cutoff:
			final_genes.append( gene )
	return final_genes


def load_mapping_table( filename ):	#NEW
	"""! @brief load mapping table from given file """
	
	mapping_table = {}
	with open( filename, "r" ) as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			mapping_table.update( { parts[0]: parts[1] } )
			line = f.readline()
	return mapping_table


if __name__ == '__main__':
	
	data_file ="FPKMs.txt"
	prefix = "MADS/"
	candidate_gene_file = "MADS.txt"
	name_file = "names.txt"	#NEW
	
	candidate_samples = ["010616", "020616", "040616", "060616", "090616", "120616", "140616", "160616", "180616", "210616", "240616", "280616", "260716", "040816", "110816", "230816", "080916", "220916", "031116"]
	
	exp_cutoff = 1
	
	NMT = load_mapping_table( name_file )	#NEW
	
	candidate_name_mapping_table = {}
	with open( candidate_gene_file, "r" ) as f:
		content = f.read()
		genes = list( set( re.findall( "VIT_2\d+s\d+g\d+", content ) ) )
		for gene in genes:
			try:	#NEW
				candidate_name_mapping_table.update( { gene: gene + "-" + NMT[ gene ] } )	#NEW
			except KeyError:	#NEW
				print "ERROR: no name mapped - " + gene
				candidate_name_mapping_table.update( { gene: gene } )	#NEW
	candidate_genes = sorted( candidate_name_mapping_table.keys() )
	
	# --- load data --- #
	expression_data = load_expression( data_file )	
	candidate_genes = filter_genes( expression_data, candidate_genes, exp_cutoff )	#just order by gene ID
	
	# --- adjust format for heatmap construction --- #
	genes, tissues, datamatrix = construct_data_output_file( expression_data, candidate_genes, candidate_name_mapping_table, candidate_samples )
	
	
	# -- construct heatmap via seaborn --- #
	heatmap_file = prefix + "AdditionalFileG.pdf"
	construct_heatmap( datamatrix, genes, tissues, heatmap_file )
	
	print "all done!"
