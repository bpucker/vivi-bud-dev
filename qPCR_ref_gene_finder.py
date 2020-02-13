### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.2 ###

__usage__ = """
					python qPCR_ref_gene_finder.py
					--exp <FULL_PATH_TO_EXPRESSION_FILE>
					--out <FULL_PATH_TO_OUTPUT_FILE>
					--sum_cut <CUTOFF_FOR_SUM_OVER_ALL_SAMPLES> [100]
					--num_cut <NUMBER_OF_CANDIDATES_TO_REPORT> [100]
					--anno <FULL_PATH_TO_ANNOTATION_FILE>
					
					bug reports and feature requests: bpucker@cebitec.uni-bielefeld.de
					"""
import numpy as np
from operator import itemgetter
import sys

# --- end of imports --- #

def load_expression_values( data_file ):
	"""! @brief load expression values """
	
	exp_data = {}
	with open( data_file, "r" ) as f:
		headers = f.readline().strip().split('\t')[1:]
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			exp_data.update( { parts[0]: map( float, parts[1:] ) } )
			line = f.readline()
	return exp_data


def load_annotation( annotation_file ):
	"""! @brief load annotation from given file """
	
	annotation = {}
	with open( annotation_file, "r" ) as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			annotation.update( { parts[0]: parts[1] } )
			line = f.readline()
	return annotation


def main( arguments ):
	"""! @brief runs everything """
	
	data_file = arguments[ arguments.index( '--exp' ) +1 ]
	report_file = arguments[ arguments.index( '--out' ) +1 ]
	
	if '--sum_cut' in arguments:
		exp_cutoff = int( arguments[ arguments.index( '--sum_cut' ) +1 ] )
	else:
		exp_cutoff = 100
	
	if '--num_cut' in arguments:
		num_cutoff = int( arguments[ arguments.index( '--num_cut' ) +1 ] )
	else:
		num_cutoff = 100
	
	if '--anno' in arguments:
		annotation_file = arguments[ arguments.index( '--anno' ) +1 ]
		annotation = load_annotation( annotation_file )
	else:
		annotation = {}
	
	expression = load_expression_values( data_file )
	
	gene_characteristics = []
	for gene in expression.keys():
		avg = np.median( expression[ gene ] )
		if avg > 0:
			dev = np.std( expression[ gene ] ) / avg
			if avg > exp_cutoff and min(  expression[ gene ] ) > exp_cutoff:
				gene_characteristics.append( { 'id': gene, 'avg': avg, 'dev': dev } )
	gene_characteristics.sort( key=itemgetter( 'avg' ), reverse=True )
	gene_characteristics.sort( key=itemgetter( 'dev' ) )
	
	
	with open( report_file, "w" ) as out:
		out.write( "GeneID\tNormalizedStandardDeviation\tAverageExpression\n" )
		for idx, candidate in enumerate( gene_characteristics ):
			if idx <= num_cutoff:
				try:
					out.write( "\t".join( map( str, [ candidate['id'],  candidate['dev'], candidate['avg'], annotation[ candidate['id'] ] ] ) ) + '\n' )
				except KeyError:
					out.write( "\t".join( map( str, [ candidate['id'],  candidate['dev'], candidate['avg'], "no annotation available" ] ) ) + '\n' )


if __name__ == '__main__':
	
	if '--exp' in sys.argv and '--out' in sys.argv:
		main( sys.argv )
	else:
		sys.exit( __usage__ )
