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


def plot_expression( exp_data, genes, result_file ):
	"""! @brief creates expression plot and adjusts settings """
	
	dates_order = ["010616", "020616", "040616", "060616", "090616", "120616", "140616", "160616", "180616", "210616", "240616", "280616", "260716", "040816", "110816", "230816", "080916", "220916", "031116"]
	labels = ["01 Jun", "02 Jun", "04 Jun", "06 Jun", "09 Jun", "12 Jun", "14 Jun", "16 Jun", "18 Jun", "21 Jun", "24 Jun", "28 Jun", "26 Jul", "04 Aug", "11 Aug", "23 Aug", "08 Sep", "22 Sep", "03 Nov"]
	temp_data=[16.19,16.46,17.5,18.74,20.47,20.43,16.98,17.32,17.23,16.59,16.7,16.08,14.98,14.25,14.19,14.04,14.53,14.4,16.01,17.48,21.38, 25.3,25.71,19.31,16.86,16.34,17.73,19.35,19.01,19.96,16.8,14.87,18.2,19.05,17.08,18.8,20.93,21.93,23.65,21.89,19.16,16.28,15,15.26,17.7, 21.02,22.35,23.1,25.08,21.24,21.6,20.64,21.54,22.95,21.78,21.36,20.4,21.03,21,20.28,18.22,17.7,21.06,18.73,18.11,18.56,19.51,20.05,17.78, 14.22,14.08,17.37,20.75,21.13,21.4,19.93,19.98,17.98,19.49,16.29,16.15,15.97,19.75,23.49,23.84,24.65,24.93,25.25,21.68,20.1,20.99,21.17,
	19.65,20.56,17.79,18.15,17.91,19.64,20.43,20.54,21.41,21.74,21.68,22.46,23.43,20.86,15.59,16.23,16.27,15.49,14.26,13.07,12.71,14.54,14.7,
	15.91,16.42,16.08,16.39,19.17,16.73,14.75,12.95,11.85,10.44,9.96,9.38,10.95,8.84,7.34,7.49,7.87,8.28,7.92,9.49,11.14,9.86,11.66,12.6,9.51,
	8.31,6.41,7.16,6.68,10.36,11.92,9.04,8.07,8.51,10.89,8.13,9.05,6.63,8.23,4.58]
	rain_data=[2.86,13.31,17.14,0,0,0,11.7,0,0,6.44,3.78,16.31,3.23,3.2,0.17,6.99,0.54,0,2.62,1.74,0,0,0,4.07,0,0,2.48,0,0,0,0.3,0,0,0,0,0,0,0,0, 6.65,4.25,0.86,0,0,0,0,0,0,0,7.82,27.22,0.17,0,0,0.44,0,0,0,0,0.23,0,0.38,0,14.5,0,0,0,0,0.28,0,0.19,0.06,0,0,0,0,0,0,0.13,5.54,0.16,0,0,0,0, 0,0,0,0,0,0,0,0,0,1.23,0,0,0,0,0,0,0,0,0,0,0.35,0.2,2.42,9.56,0,0,0,0,0,0,0,0,0.06,0,0,0,11.58,1.18,0,0,0,0,0,0.08,0,0,0,0,0,0,0,0,0.34,0.47,
	4.21,11.94,0,0.18,3.95,10.61,3.38,0,0,0,0,0,0,0,0,0]
	
	
	fig, ax  = plt.subplots( figsize = (10, 5) )
	
	ticks = [ 0, 1, 3, 5, 8, 11, 13, 15, 17, 20, 23, 27, 55, 64, 71, 83, 99, 113, 156 ]
	
	# --- add weather data to figure --- #
	ax2 = ax.twinx()
	ax3 = ax.twinx()
	ax2.plot( range( 1, 156, 1 ), temp_data, "--", color="red", alpha=0.5 )
	for idx, value in enumerate( rain_data ):
		ax3.plot( [ idx+0.2, idx+0.2 ], [ 0, rain_data[ idx ] ], color="blue", alpha=0.5 )
	
	# --- add gene expression --- #
	all_values = []
	for gene in genes:
		mean_values_to_plot = []
		for idx, date in enumerate( dates_order ):
			mean_values_to_plot.append( np.mean( exp_data[ gene ][ date ] ) )
			ax.plot( [ ticks[ idx ], ticks[ idx ] ], [ min( exp_data[ gene ][ date ] ), max( exp_data[ gene ][ date ] ) ], color="grey", linewidth=0.1 )
			all_values.append( max( exp_data[ gene ][ date ] ) )
		ax.plot( ticks, mean_values_to_plot, "--o", color="lime", linewidth=1 )
	
	# --- add legend --- #
	red_patch = mpatches.Patch( color='red', label='temperature', alpha=0.5 )
	blue_patch = mpatches.Patch( color='blue', label='rain', alpha=0.5 ) 
	green_patch = mpatches.Patch( color='lime', label='gene expression' ) 
	ax.legend( handles=[ red_patch, blue_patch, green_patch ] )
	
	# --- adjust figure properties to make it pretty --- #
	ax.set_xticks(ticks)
	ax.set_xticklabels(labels, {'rotation': 90, 'fontsize': 8})
	ax.set_xlabel("Days")
	ax.set_ylabel("gene expression [Tags Per Million assigned tags]")
	
	ax.set_xlim( -0.5, 156.5 )
	ax.set_ylim( 0, max( all_values ) )
	ax2.set_ylim( 0, max( temp_data ) )
	ax3.set_ylim( 0, max( rain_data ) )
	
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
		
	fig.subplots_adjust( left=0.0675, right=0.999, top=0.999, bottom=0.15 ) 
	fig.savefig( result_file, dpi = 300 )
	plt.close( "all" )


def filter_genes( exp_data, genes, exp_cutoff ):
	"""! @brief filter genes based on expression """
	
	final_genes = []
	for gene in genes:
		exp = [ x for sublist in exp_data[ gene ].values() for x in sublist ]
		if np.mean( exp ) >= exp_cutoff:
			final_genes.append( gene )
	return final_genes


def analyze_correlation( exp_data, genes, cor_file ):
	"""! @brief get correlation between temperature and expression values of heatshock genes """
	
	# --- merge expression data by timepoint --- #
	merged_exp_data = {}
	for gene in genes:
		exp = exp_data[ gene ]
		for key in exp.keys():
			try:
				merged_exp_data[ key ] += exp[ key ]
			except KeyError:
				merged_exp_data.update( { key: exp[ key ] } )
		
	# --- check correlation --- #
	dates_order = [ "020616", "040616", "060616", "090616", "120616", "140616", "160616", "180616", "210616", "240616", "280616", "260716", "040816", "110816", "230816", "080916", "220916", "031116"]
	ticks = [ 1, 3, 5, 8, 11, 13, 15, 17, 20, 23, 27, 55, 64, 71, 83, 99, 113, 156 ]
	temp_data=[16.19,16.46,17.5,18.74,20.47,20.43,16.98,17.32,17.23,16.59,16.7,16.08,14.98,14.25,14.19,14.04,14.53,14.4,16.01,17.48,21.38, 25.3,25.71,19.31,16.86,16.34,17.73,19.35,19.01,19.96,16.8,14.87,18.2,19.05,17.08,18.8,20.93,21.93,23.65,21.89,19.16,16.28,15,15.26,17.7, 21.02,22.35,23.1,25.08,21.24,21.6,20.64,21.54,22.95,21.78,21.36,20.4,21.03,21,20.28,18.22,17.7,21.06,18.73,18.11,18.56,19.51,20.05,17.78, 14.22,14.08,17.37,20.75,21.13,21.4,19.93,19.98,17.98,19.49,16.29,16.15,15.97,19.75,23.49,23.84,24.65,24.93,25.25,21.68,20.1,20.99,21.17,
	19.65,20.56,17.79,18.15,17.91,19.64,20.43,20.54,21.41,21.74,21.68,22.46,23.43,20.86,15.59,16.23,16.27,15.49,14.26,13.07,12.71,14.54,14.7,
	15.91,16.42,16.08,16.39,19.17,16.73,14.75,12.95,11.85,10.44,9.96,9.38,10.95,8.84,7.34,7.49,7.87,8.28,7.92,9.49,11.14,9.86,11.66,12.6,9.51,
	8.31,6.41,7.16,6.68,10.36,11.92,9.04,8.07,8.51,10.89,8.13,9.05,6.63,8.23,4.58]
	
	print len( temp_data )
	
	x_values = []	#temperature
	y_values = []	#gene expression
	for idx, date in enumerate( dates_order[1:] ):
		y_values.append( np.mean( merged_exp_data[ date ] ) )
		x_values.append( temp_data[ ticks[ idx ] ] )
	
	print stats.pearsonr( x_values, y_values )
	print stats.spearmanr( x_values, y_values )
	
	
	fig, ax = plt.subplots()
	ax.scatter( x_values, y_values, s=10, color="lime" )
	ax.set_xlabel( "temperature" )
	ax.set_ylabel( "heatshock gene expression" )
	fig.savefig( cor_file, dpi=300 )


def plot_merged_expression( exp_data, genes, merged_result_file  ):
	"""! @brief plot merged expression of all genes """
		
	counter = 0	#count genes in plot
	surviving_genes = []
	# --- merge expression data by timepoint --- #
	merged_exp_data = {}
	for gene in genes:
		exp = exp_data[ gene ]
		status = False
		for key in exp.keys():
			try:
				merged_exp_data[ key ] += exp[ key ]
				if not status:
					counter += 1
					status = True
					surviving_genes.append( gene )
			except KeyError:
				merged_exp_data.update( { key: exp[ key ] } )
	
	print "number of genes in plot: " + str( counter )
	print sorted( list( set( surviving_genes ) ) )
	
	# --- fixed data --- #
	dates_order = ["010616", "020616", "040616", "060616", "090616", "120616", "140616", "160616", "180616", "210616", "240616", "280616", "260716", "040816", "110816", "230816", "080916", "220916", "031116"]
	labels = ["01 Jun", "02 Jun", "04 Jun", "06 Jun", "09 Jun", "12 Jun", "14 Jun", "16 Jun", "18 Jun", "21 Jun", "24 Jun", "28 Jun", "26 Jul", "04 Aug", "11 Aug", "23 Aug", "08 Sep", "22 Sep", "03 Nov"]
	temp_data=[16.19,16.46,17.5,18.74,20.47,20.43,16.98,17.32,17.23,16.59,16.7,16.08,14.98,14.25,14.19,14.04,14.53,14.4,16.01,17.48,21.38, 25.3,25.71,19.31,16.86,16.34,17.73,19.35,19.01,19.96,16.8,14.87,18.2,19.05,17.08,18.8,20.93,21.93,23.65,21.89,19.16,16.28,15,15.26,17.7, 21.02,22.35,23.1,25.08,21.24,21.6,20.64,21.54,22.95,21.78,21.36,20.4,21.03,21,20.28,18.22,17.7,21.06,18.73,18.11,18.56,19.51,20.05,17.78, 14.22,14.08,17.37,20.75,21.13,21.4,19.93,19.98,17.98,19.49,16.29,16.15,15.97,19.75,23.49,23.84,24.65,24.93,25.25,21.68,20.1,20.99,21.17,
	19.65,20.56,17.79,18.15,17.91,19.64,20.43,20.54,21.41,21.74,21.68,22.46,23.43,20.86,15.59,16.23,16.27,15.49,14.26,13.07,12.71,14.54,14.7,
	15.91,16.42,16.08,16.39,19.17,16.73,14.75,12.95,11.85,10.44,9.96,9.38,10.95,8.84,7.34,7.49,7.87,8.28,7.92,9.49,11.14,9.86,11.66,12.6,9.51,
	8.31,6.41,7.16,6.68,10.36,11.92,9.04,8.07,8.51,10.89,8.13,9.05,6.63,8.23,4.58]
	# rain_data=[2.86,13.31,17.14,0,0,0,11.7,0,0,6.44,3.78,16.31,3.23,3.2,0.17,6.99,0.54,0,2.62,1.74,0,0,0,4.07,0,0,2.48,0,0,0,0.3,0,0,0,0,0,0,0,0, 6.65,4.25,0.86,0,0,0,0,0,0,0,7.82,27.22,0.17,0,0,0.44,0,0,0,0,0.23,0,0.38,0,14.5,0,0,0,0,0.28,0,0.19,0.06,0,0,0,0,0,0,0.13,5.54,0.16,0,0,0,0, 0,0,0,0,0,0,0,0,0,1.23,0,0,0,0,0,0,0,0,0,0,0.35,0.2,2.42,9.56,0,0,0,0,0,0,0,0,0.06,0,0,0,11.58,1.18,0,0,0,0,0,0.08,0,0,0,0,0,0,0,0,0.34,0.47,
	# 4.21,11.94,0,0.18,3.95,10.61,3.38,0,0,0,0,0,0,0,0,0]
	
	fig, ax  = plt.subplots( figsize = (10, 5) )
	
	ticks = [ 0, 1, 3, 5, 8, 11, 13, 15, 17, 20, 23, 27, 55, 64, 71, 83, 99, 113, 156 ]
	
	# --- add weather data to figure --- #
	ax2 = ax.twinx()
	ax3 = ax.twinx()
	ax2.plot( range( 1, 156, 1 ), temp_data, "--", color="red", alpha=0.5 )
	# for idx, value in enumerate( rain_data ):
		# ax3.plot( [ idx+0.2, idx+0.2 ], [ 0, rain_data[ idx ] ], color="blue", alpha=0.5 )
	
	
	# --- add gene expression --- #
	x_values = []	#date
	y_values = []	#gene expression
	y_dev_values = []	#standard deviation
	for idx, date in enumerate( dates_order ):
		y_values.append( np.mean( merged_exp_data[ date ] ) )
		x_values.append( dates_order.index( date ) )
		#y_dev_values.append( np.std( merged_exp_data[ date ] ) )
		#y_dev_values.append( [ sorted( merged_exp_data[ date ] )[ int( 0.25*len( merged_exp_data[ date ] ) ) ], sorted( merged_exp_data[ date ] )[ int( 0.75*len( merged_exp_data[ date ] ) ) ] ] )
	
	# all_values = []
	# for idx, date in enumerate( dates_order ):
		# ax.plot( [ ticks[ idx ], ticks[ idx ] ], y_dev_values[ idx ], color="grey", linewidth=0.5 )
		# all_values.append( max( y_dev_values[ idx ] ) )
	ax.plot( ticks, y_values, "--o", color="lime", linewidth=1 )
	
	# --- add legend --- #
	red_patch = mpatches.Patch( color='red', label='temperature', alpha=0.5 )
	#blue_patch = mpatches.Patch( color='blue', label='rain', alpha=0.5 ) 
	green_patch = mpatches.Patch( color='lime', label='heat shock gene expression' ) 
	ax.legend( handles=[ red_patch, green_patch ] )	#blue_patch
	
	# --- adjust figure properties to make it pretty --- #
	ax.set_xticks(ticks)
	ax.set_xticklabels(labels, {'rotation': 90, 'fontsize': 8})
	ax.set_xlabel("Days")
	ax.set_ylabel("gene expression [Tags Per Million assigned tags]")
	ax2.set_ylabel( "temperature [C]" )
	#ax3.set_ylabel( "precipitation [mm]" )
	
	ax.set_xlim( -0.5, 156.5 )
	ax.set_ylim( 0, max( y_values )+10 )	#all_values
	ax2.set_ylim( 0, 30 )
	ax3.set_ylim( 0, 30 )
	
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	
	# lspine = ax2.spines['right']
	# lspine.set_position( ('axes', 1.1 ) )
	
	fig.subplots_adjust( left=0.0675, right=0.85, top=0.98, bottom=0.15 ) 
	fig.savefig( merged_result_file, dpi = 300 )
	plt.close( "all" )


def merged_expression_boxplots( exp_data, genes, boxplot_result_file  ):
	"""! @brief plot merged expression of all genes """
	
	
	counter = 0	#count genes in plot
	# --- merge expression data by timepoint --- #
	merged_exp_data = {}
	for gene in genes:
		exp = exp_data[ gene ]
		status = False
		for key in exp.keys():
			try:
				merged_exp_data[ key ] += exp[ key ]
				if not status:
					counter += 1
					status = True
			except KeyError:
				merged_exp_data.update( { key: exp[ key ] } )
	
	# --- fixed data --- #
	dates_order = ["010616", "020616", "040616", "060616", "090616", "120616", "140616", "160616", "180616", "210616", "240616", "280616", "260716", "040816", "110816", "230816", "080916", "220916", "031116"]
	labels = ["01 Jun", "02 Jun", "04 Jun", "06 Jun", "09 Jun", "12 Jun", "14 Jun", "16 Jun", "18 Jun", "21 Jun", "24 Jun", "28 Jun", "26 Jul", "04 Aug", "11 Aug", "23 Aug", "08 Sep", "22 Sep", "03 Nov"]
	
	fig, ax  = plt.subplots( figsize = (10, 5) )
	
	# --- add boxplots --- #
	boxplot_values = []
	for date in dates_order:
		boxplot_values.append( merged_exp_data[ date ] )
	ax.boxplot( boxplot_values, labels=labels, showfliers=False, showmeans=True )
	
	
	# --- adjust figure properties to make it pretty --- #
	ax.set_xlabel("Days")
	ax.set_ylabel("gene expression [Tags Per Million assigned tags]")
	
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
		
	fig.subplots_adjust( left=0.0675, right=0.999, top=0.999, bottom=0.15 ) 
	fig.savefig( boxplot_result_file, dpi = 300 )
	plt.close( "all" )


def main(arguments):
	"""! @brief gets all arguments from command line and calls functions for collecting, sorting and averaging TPMs and creating the expression plot """
	exp_file = arguments[arguments.index("--exp_file" ) + 1]
	genes_file = arguments[arguments.index("--genes_file" ) + 1]
	directory = arguments[arguments.index("--out" ) + 1]
	
	if '--exp_cutoff' in arguments:
		exp_cutoff = int( arguments[arguments.index("--exp_cutoff" ) + 1] )
	else:
		exp_cutoff = 30
	
	if directory[-1] != "/":
		directory += "/"
	if not os.path.exists( directory ):
		os.makedirs( directory )
	
	result_file = directory + "expression_plot.png"
	cor_file = directory + "correlation.png"
	merged_result_file = directory + "merged_expression.png"
	boxplot_result_file = directory + "merged_expression_boxplot.png"
	
	result_file_control = directory + "expression_plot_control.png"
	cor_file_control = directory + "correlation_control.png"
	
	exp_data = load_expression( exp_file )
	
	# --- analysing candidate genes --- #
	genes = load_genes( genes_file )
	print str( len( genes) ) + " genes subjected to analysis"
	genes = filter_genes( exp_data, genes, exp_cutoff )
	print str( len( genes) ) + " genes survived filtering"
	plot_expression( exp_data, genes, result_file )
	analyze_correlation( exp_data, genes, cor_file )
	plot_merged_expression( exp_data, genes, merged_result_file  )
	merged_expression_boxplots( exp_data, genes, boxplot_result_file  )
	
	# # --- running analysis on set of random genes --- #
	# random_genes = random.sample( exp_data.keys(), k=len( genes ) )
	# random_genes = filter_genes( exp_data, random_genes, exp_cutoff )
	# print str( len( random_genes) ) + " random genes survived filtering"
	# plot_expression( exp_data, random_genes, result_file_control )
	# analyze_correlation( exp_data, random_genes, cor_file_control )
	


if __name__ == "__main__":
	if "--exp_file" in sys.argv and "--genes_file" in sys.argv and '--out' in sys.argv:
		main(sys.argv)
	else:
		sys.exit(__usage__)
