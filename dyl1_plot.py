### Sarah Becker and Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.1 ###

__usage__ = """ 	python plot_triplicates_one_y_axis.py
							--exp_file <INPUT_TPMS_FILE>
							--genes_file <GENE_IDS_OF_GENES_TO_SHOW_IN_PLOT,ONE_GENE_ID_PER_LINE>
							--out <OUTPUT_DIRECTORY>
						"""

import matplotlib.pyplot as plt
import sys, os, re, math
import numpy as np
import matplotlib.patches as mpatches
from scipy import stats
import random

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
			data.update( { parts[0]: exp } )
			line = f.readline()
	return data


def load_genes( genes_file ):
	"""! @brief load genes """
	
	genes = []
	with open( genes_file, "r" ) as f:
		line = f.readline()
		while line:
			if len( line.strip() ) > 10:
				genes.append( line.strip() )
			line = f.readline()
	return genes


def load_day_lengths( day_len_file ):
	"""! @brief load day lenths per day """
	
	day_lengths = []
	with open( day_len_file, "r" ) as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			day_lengths.append( int( parts[1] ) / 60.0 )	#convert minutes to hours
			line = f.readline()
	
	return day_lengths


def plot_expression( exp_data, genes, result_file, day_lengths ):
	"""! @brief creates expression plot and adjusts settings """
	
	dates_order = ["010616", "020616", "040616", "060616", "090616", "120616", "140616", "160616", "180616", "210616", "240616", "280616", "260716", "040816", "110816", "230816", "080916", "220916", "031116"]
	labels = ["01 Jun", "", "04 Jun", "06 Jun", "09 Jun", "12 Jun", "14 Jun", "16 Jun", "18 Jun", "21 Jun", "24 Jun", "28 Jun", "26 Jul", "04 Aug", "11 Aug", "23 Aug", "08 Sep", "22 Sep", "03 Nov"]
	#02 Jun removed from labels
		
	fig, ax  = plt.subplots( figsize = (10, 5) )
	
	ticks = [ 0, 1, 3, 5, 8, 11, 13, 15, 17, 20, 23, 27, 55, 64, 71, 83, 99, 113, 155 ]
	
	# --- add gene expression --- #
	all_values = []
	colors = [ "red", "blue" ]
	name_mapping = { 'VIT_210s0003g00090': 'VviDRM1' }
	patches = []
	for k, gene in enumerate( genes ):
		mean_values_to_plot = []
		for idx, date in enumerate( dates_order ):
			mean_values_to_plot.append( np.mean( exp_data[ gene ][ date ] ) )
			ax.plot( [ ticks[ idx ], ticks[ idx ] ], [ min( exp_data[ gene ][ date ] ), max( exp_data[ gene ][ date ] ) ], color="grey", linewidth=1 )
			all_values.append( max( exp_data[ gene ][ date ] ) )
		ax.plot( ticks, mean_values_to_plot, "--o", color=colors[ k ], linewidth=1 )
		patches.append( mpatches.Patch( color=colors[ k ], label= "$\it{" + name_mapping[ gene ] +'}$ expression' ) )
	
	# --- add day lengths --- #
	ax2 = ax.twinx()
	ax2.plot( range( 156 ), day_lengths, color="black" )
	ax2.set_ylim( 6, 18 )
	ax2.set_ylabel( "day length [h]" )
	patches.append( mpatches.Patch( color='black', label='day length' ) )
	
	
	# --- add legend --- #
	ax.legend( handles=patches )
	
	# --- adjust figure properties to make it pretty --- #
	ax.set_xticks(ticks)
	ax.set_xticklabels(labels, {'rotation': 90, 'fontsize': 8})
	ax.set_xlabel("Days")
	ax.set_ylabel("transcript level [Counts Per Million]")
	
	ax.set_xlim( -0.5, 156.5 )
	ax.set_ylim( 0, max( all_values ) )
	
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	
	fig.subplots_adjust( left=0.0675, right=0.92, top=0.99, bottom=0.15 ) 
	fig.savefig( result_file, dpi = 600 )
	plt.close( "all" )


def main(arguments):
	"""! @brief gets all arguments from command line and calls functions for collecting, sorting and averaging TPMs and creating the expression plot """
	exp_file = arguments[arguments.index("--exp_file" ) + 1]
	genes_file = arguments[arguments.index("--genes_file" ) + 1]
	directory = arguments[arguments.index("--out" ) + 1]
	day_len_file = arguments[arguments.index("--days" ) + 1]
	
	if directory[-1] != "/":
		directory += "/"
	if not os.path.exists( directory ):
		os.makedirs( directory )
	
	exp_cutoff = 30
	
	result_file = directory + "VviDRM1_expression.jpg"
	
	exp_data = load_expression( exp_file )
	day_lengths = load_day_lengths( day_len_file )
	
	# --- analysing candidate genes --- #
	genes = load_genes( genes_file )
	print str( len( genes) ) + " genes subjected to analysis"
	plot_expression( exp_data, genes, result_file, day_lengths )
	

if __name__ == "__main__":
	if "--exp_file" in sys.argv and "--genes_file" in sys.argv and '--out' in sys.argv:
		main(sys.argv)
	else:
		sys.exit(__usage__)
