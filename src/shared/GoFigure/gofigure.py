##
## Summary Visualizations of Gene Ontology Enrichment with GO-Figure!
##
## Authors: Maarten JMF Reijnders and Robert M Waterhouse
##
## Contact: maarten.reijnders@unil.ch or robert.waterhouse@unil.ch
##

import matplotlib
matplotlib.use('Agg')
import sys
import logging
from datetime import date
from collections import defaultdict
from collections import OrderedDict
from sklearn.feature_extraction import DictVectorizer
import numpy as np
from sklearn.manifold import MDS
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
from math import log
from math import log10
import matplotlib.patches as mpatches
import os
import argparse
from adjustText import adjust_text ## Python package which rearranges plot labels to ensure minimal overlap
from matplotlib.colors import ListedColormap

def warn(*args, **kwargs):
    pass

import warnings
warnings.warn = warn

scriptPath = os.path.dirname(os.path.abspath(__file__))

## The below relative script path is only used for building the GO-Figure one-file binaries 
#scriptPath = sys.executable.replace('/gofigure','')

#scriptPath = os.path.dirname(os.path.abspath(__file__)) ## Get the absolute path of this script. Used to refer to data files used by this script.

## Long list of various arguments used by this script. Only an input file and output directory are required.
def read_arguments():
	parser = argparse.ArgumentParser(description="GOFigure!: GO enrichment visualization based on semantic similarity.",add_help=True)
	required = parser.add_argument_group('Required arguments')
	optional = parser.add_argument_group('Optional arguments')
	input_arguments = parser.add_argument_group('Optional input file and data handling arguments')
	legend_arguments = parser.add_argument_group('Optional legend arguments')
	scatter_arguments = parser.add_argument_group('Optional scatterplot arguments')
	title_arguments = parser.add_argument_group('Optional title arguments')
	output_arguments = parser.add_argument_group('Optional outfile arguments')
	
	required.add_argument('-i','--input',help="Input file in tabular format with the columns: GO ID + P-value for standard input, GO ID + P-Value + Significant for standard-plus input, TopGO output as an input, or GoStats output as an input. Can use # or ! or %% to comment out lines. Change --input_type accordingly. Default input is 'standard'. Example input files are found in the root directory.",required=True,action='store')
	
	required.add_argument('-o','--output',help='Output directory',required=True,action='store')

	## Input arguments
	input_arguments.add_argument('-a','--max_clusters',help="Maximum amount of clusters to plot (integer value). Default = 50.",required=False,action='store',default=50,type=int)
	
	input_arguments.add_argument('-j','--input_type', help="Type of input file. Use 'topgo' for standard TopGO input, 'gostats' for standard GOStats input, 'standard' for an input file where the first column is the GO ID and the second is the p-value, and 'standard-plus' for standard input but with a third column containing a user defined numerical value (for TopGO and GOStats input, this is the 'significant' value). Default = standard", required=False,action='store',default='standard',choices=['topgo','gostats','standard','standard-plus'])

	input_arguments.add_argument('-n','--ontology',help="Which ontology to use: biological process ('bpo'), molecular function ('mfo'), cellular component ('cco'), or all ontologies ('all'). Default is all.", required=False,action='store',default='all',choices=['bpo','mfo','cco','all'])

	input_arguments.add_argument('-r','--representatives',help="A list of GO terms that have priority for being a representative of a cluster. I.e. if one term in a cluster has priority, that term will always be the representative. If two terms in a cluster have priority, only those two will be considered. Input is a list of GO terms separated by a comma, such as 'GO:0000001,GO:0000002'.",required=False,action='store',type=str)

	input_arguments.add_argument('-si','--similarity_cutoff',help='Similarity cutoff to be used between GO terms, a value between 0 and 1. Default = 0.5.',required=False,action='store',default=0.5,type=float)

	input_arguments.add_argument('-v','--max_pvalue',help="Maximum p-value to consider (floating value). Default = 99999.99999.",required=False,action='store',default=99999,type=float)
	
	input_arguments.add_argument('-so','--sort_by',help="Which values to use for sorting the clusters: 'pval' (p-value) or 'user' (user-defined value) or 'user-descending' (user-defined value descending). Default = pval.",required=False,action='store',default='pval',type=str,choices=['pval','user','user-descending'])

	input_arguments.add_argument('-su','--sum_user',help="To sum the user-defined values (column 3) for each member of a cluster. Either 'True' or 'False'. Default = False.",required=False,default=False,action='store',type=bool,choices=[True,False])

	input_arguments.add_argument('-to','--top_level',help="Set top level GO terms for clusters as given by the GO DAG (see https://www.ebi.ac.uk/QuickGO). Top level GO terms have to be given separated by comma's, without spaces. E.g.: 'GO:000001,GO:000008'.",required=False,action='store',type=str)

	input_arguments.add_argument('-nc','--name_changes',help="A list of GO terms that will be used as representative of a cluster specifically for naming purposes, but not for internal calculations. This is opposed to the'--representatives option, which provides GO terms to be used as representatives of a cluster both internally and externally. This specific option allows for the renaming of clusters without recalculating the clusters when there is a need to reproduce the original figure. Input is a list of GO terms separated by a comma, such as 'GO:0000001,GO:0000002'.",required=False,action='store',type=str)

	## Scatterplot layout arguments
	scatter_arguments.add_argument('-b','--opacity',help="Set opacity for the clusters, floating point from 0 to 1. Default = 1", required=False,action='store',default=1.0,type=float)
	
	scatter_arguments.add_argument('-c','--colours',help="Color GO clusters based on p-value ('pval'), log10 p-value ('log10-pval'), number of GO terms that are a member of the cluster ('members'), frequency of the GO term in the GOA database ('frequency'), a unique colour per cluster ('unique'), or a user defined value ('user'). Default = 'log10-pval'.",required=False,action='store',default='log10-pval',choices=['pval','log10-pval','user','members','frequency','unique'])
	
	scatter_arguments.add_argument('-p','--palette',help='Set to any color brewer palette available for MatPlotLib (https://matplotlib.org/3.1.1/tutorials/colors/colormaps.html). Default = plasma',required=False,action='store',default='plasma')

	scatter_arguments.add_argument('-rs','--random_state',help="Set random state for the calculation of the X and Y label. Useful if you want the figures to be exactly the same. Needs to be an integer or None. Specifying 'None' will create a different orientation of the plot every time. Default = 1",required=False,action='store',default=1)
	scatter_arguments.add_argument('-s','--size',help="Size of GO clusters based on how many GO terms are a member of the cluster ('members'), frequency in GOA database ('frequency'), p-value where smaller = larger size ('pval'), the user defined value ('user'), or a fixed integer for every cluster. Default = 'members'.",required=False,action='store',default='members')

	scatter_arguments.add_argument('-sr','--size_range',help="Set cluster size range to 'small', 'medium', or 'big'. Default = medium.", required=False,action='store',default='medium',choices=['small','medium','big'])

	scatter_arguments.add_argument('-u','--cluster_labels',help="Label clusters numbered based on the sorting option ('numbered'), based on the representative GO term ('go'), based on the representative GO term with arrows ('go-arrows') based on the representative GO term name ('description'), or based on the representative GO term name with arrows ('description-numbered'). Default = numbered.",required=False,action='store',default='numbered',choices=['numbered','go','go-arrows','description','description-numbered'])
	
	scatter_arguments.add_argument('-x','--xlabel',help="X-axis label. Default = 'Semantic space X'",required=False,action='store',default='Semantic space X',type=str)

	scatter_arguments.add_argument('-y','--ylabel',help="Y-axis label. Default = 'Semantic space Y'",required=False,action='store',default='Semantic space Y',type=str)

	scatter_arguments.add_argument('-z','--colour_label',help="Colourbar label. Default is the description of the input used for --colours",required=False,action='store',type=str)

	## Legend layout arguments
	legend_arguments.add_argument('-d','--legend_position',help="Position the legend at the bottom left ('left'), bottom right ('right'), or bottom center ('center'). Default = center.",required=False,action='store',default='center',choices=['left','right','center'])
	
	legend_arguments.add_argument('-e','--description_limit',help="Integer character limit of GO term descriptions in the legend. Default = 35.",required=False,action='store',default=35,type=int)
	
	legend_arguments.add_argument('-f','--font_size',help="Font size of the legend. 'xx-small', 'x-small', 'small', 'medium', 'large' , 'x-large', or 'xx-large'. Default = medium",required=False,action='store',default='medium',choices=['xx-small','x-small','small','medium','large','x-large','xx-large'])

	legend_arguments.add_argument('-g','--legend_columns',help="Legend as a single ('single') column or double ('double') column. Default = double.",required=False,action='store',default='double',choices=['single','double','triple'])

	legend_arguments.add_argument('-l','--legend',help="Option to show GO terms and descriptions in the legend ('full'), GO term only 'go', description only ('description'), or no legend ('none'). Default = description.",required=False,action='store',default='description',choices=['full','go','description','none'])

	legend_arguments.add_argument('-m','--max_label',help='Maximum labels to display in the legend. Default = 10.',required=False,action='store',default=10,type=int)
	
	## Output files arguments
	output_arguments.add_argument('-dp','--dpi',help="Set DPI for outfiles. Default = 400",required=False,action='store',default=400,type=int)
	
	output_arguments.add_argument('-q','--file_type',help="Image file type. 'png', 'pdf', 'svg', or 'tiff'. Default = 'png'",required=False,action='store',default='png',choices=['png','pdf','tiff','svg'])

	output_arguments.add_argument('-w','--outfile_appendix',help="What to add to the outfiles after 'biological_process', 'molecular_function', or 'cellular_component'. By default it will add the value given for --outdir, replacing '/' with '_'.",required=False,action='store',default='')

	## Figure title arguments
	title_arguments.add_argument('-k','--title_size',help="Set title size. 'xx-small', 'x-small', 'small', 'medium', 'large' , 'x-large', or 'xx-large'. Default = medium",required=False,action='store',default='medium',choices=['xx-small','x-small','small','medium','large','x-large','xx-large'])

	title_arguments.add_argument('-t','--title',help="Set figure title. Has to be between single or double quotation marks.",required=False,action='store')
	
	return(parser)

## Creating logfiles. Logger_std is for writing to the console, logger is for writing to a file.
def create_logger(output_dir_path,appendix):
	logger = logging.getLogger('maartenviz')
	logger_std = logging.getLogger('maartenviz_std')
	logger.setLevel(logging.INFO)
	logger_std.setLevel(logging.INFO)

	if appendix == '':
		error_file_handler = logging.FileHandler(output_dir_path+'/logfile_'+output_dir_path.replace('/','_')+'.txt',mode='w')
	else:
		error_file_handler = logging.FileHandler(output_dir_path+'/logfile_'+appendix+'.txt',mode='w')
	error_console_handler = logging.StreamHandler()
	error_console_handler.setLevel(logging.CRITICAL)
	logger.addHandler(error_file_handler)
	logger.addHandler(error_console_handler)

	console_handler = logging.StreamHandler()

	logger_std.addHandler(console_handler)
	logger_std.propagate = False
	logger.propagate=False
	return(logger,logger_std)

## Further handling of the command line input arguments
def argparse_extended_error_handling(parser):
	required = parser.parse_args()
	optional = parser.parse_args()
	input_file_path = required.input
	output_dir_path = required.output
	max_label = optional.max_label
	palette = optional.palette
	colours = optional.colours
	size = optional.size
	description_limit = optional.description_limit
	legend = optional.legend
	legend_position = optional.legend_position
	legend_columns = optional.legend_columns
	max_pvalue = optional.max_pvalue
	max_clusters = optional.max_clusters
	cluster_labels = optional.cluster_labels
	similarity_cutoff = optional.similarity_cutoff
	input_type = optional.input_type
	ontology = optional.ontology
	opacity = optional.opacity
	size_range = optional.size_range
	font_size = optional.font_size
	title = optional.title
	title_size = optional.title_size
	dpi = optional.dpi
	file_type = optional.file_type
	outfile_appendix = optional.outfile_appendix
	xlabel = optional.xlabel
	ylabel = optional.ylabel
	zlabel = optional.colour_label
	priorities = optional.representatives
	random_state = optional.random_state
	sort_by = optional.sort_by
	sum_user = optional.sum_user
	top_level = optional.top_level
	name_changes = optional.name_changes

	## Initialize logger instances
	logger,logger_std = create_logger(output_dir_path,outfile_appendix)

	## Convert priorities input to list of GO terms
	if priorities == None:
		priorities = []
	else:
		priorities = priorities.split(',')

	## Check if priorities are valid GO terms
	for priority in priorities:
		if not priority.startswith('GO:') or not len(priority) == 10:
			logger.critical('gofigure.py: error: argument -r/--representatives contains a GO term that is not in valid GO term format: '+priority)
			exit()
	
	try:
		random_state = int(random_state)
	except:
		if random_state.lower() == 'none':
			random_state = None
		else:
			logger.critical("gofigure.py: error: argument -rs/--random_state: invalid choice: '"+str(random_state)+"' (requires an integer or 'None'")
			exit()

	## Opacity is a number from 0 to 1
	if opacity < 0 or opacity > 1:
		logger.critical("gofigure.py: error: argument -o/--opacity: invalid choice: '"+str(opacity)+"' (requires a floating point from 0 to 1)")
		exit()
	if not size in ['user','pval','frequency','members']:
		try:
			size = int(size)
		except:
			logger.critical("gofigure.py: error: argument -s/--size: invalid choice: '"+size+"' (choose from 'members', 'frequency', 'pval', 'significant', or an integer)")
			exit()

	## Similarity cutoff is a number from 0 to 1
	if similarity_cutoff < 0 or similarity_cutoff > 1:
		logger.critical("gofigure.py: error: argument -c/--colours: invalid choice: '"+str(similarity_cutoff)+"' (choose any value ranging from 0 to 1)")
		exit()

	## See if the given colour palette is valid
	try:
		cmap = ListedColormap(sns.color_palette(palette))
	except:
		logger.critical("gofigure.py: error: argument -p/--palette: invalid choice: '"+palette+"' (for example, pick one of https://matplotlib.org/3.1.1/tutorials/colors/colormaps.html)")
		exit()

	## Standard input is incompatible with user option for size and colours
	if input_type == 'standard':
		if colours == 'user':
			logger.critical("gofigure.py: error: -it/--input_type standard and -c/--colour user are incompatible")
			exit()
		if size == 'user':
			logger.critical("gofigure.py: error: -it/--input_type standard and -s/--size user are incompatible")
			exit()

	## Two bools to see if there is valid TopGO or GOStats input. One line of valid input, bool = true
	gostats_bool = False
	topgo_bool = False
	line_count = 1 ## Keeps track of which line the error is caught
	
	## Various input checks for all input types. I.e. column numbers, if there are GO terms found, etc.
	for line in open(input_file_path):
		if not line.startswith('!') and not line.startswith('#') and not line.startswith('%') and not line.strip() == '':
			ssline = line.strip().split('\t')
			if input_type.lower() == 'standard' or input_type.lower() == 'standard-plus' and 'GO:' in line:
				if not ssline[0].startswith('GO:') and not len(ssline[0]) == 10:
					logger.critical("gofigure.py: input error line "+str(line_count)+": standard or standard-plus input first column has to be a GO term. Please see 'example_data/standard_example_input.tsv' or 'example_data/standard-plus_example_input.tsv' for compatible standard or standard-plus input")
					exit()
				try:
					float(ssline[1].replace('<','').replace(',','.'))
				except:
					logger.critical("gofigure.py: input error line "+str(line_count)+": standard or standard-plus input second column is not a number. Needs to be p-value. Please see 'example_data/standard_example_input.tsv' or 'example_data/standard-plus_example_input.tsv' for compatible standard or standard-plus input")
					exit()
				if input_type.lower() == 'standard':
					if len(ssline) != 2:
						logger.critical("gofigure.py: input error: --input type given is 'standard', but column number is not 2. Please see 'example_data/standard_example_input.tsv' for correct standard input.")
						exit()
				if input_type.lower() == 'standard-plus':
					if len(ssline) != 3:
						logger.critical("gofigure.py: input error: --input_type given is 'standard-plus', but column number is not 3. Please see 'example_data/standard-plus_example_input.tsv' for correct standard-plus input.")
						exit()
					try:
						float(ssline[2].replace('<',''))
					except:
						logger.critical("gofigure.py: input error line "+str(line_count)+": standard input requires third column to be a user defined numerical value, but is not a number. Please see 'example_data/standard-plus_example_input.tsv' for compatible standard-plus input")
						exit()
			elif input_type.lower() == 'topgo':
				if len(ssline) == 7 and 'GO:' in line:
					if not ssline[1].startswith('GO:') and not len(ssline[1]) == 10:
						logger.critical("gofigure.py: input error line "+str(line_count)+": TopGO input second column needs to be a GO term. Pleas see 'example_data/topgo_example_input.tsv' for correct TopGO input.")
						exit()
					try:
						float(ssline[6].replace('<',''))
					except:
						logger.critical("gofigure.py: input error line "+str(line_count)+": TopGO input 7th column is not a number. Needs to be p-value. Please see 'example_data/topgo_example_input.tsv' for compatible TopGO input")
						exit()
					try:
						float(ssline[4].replace('<',''))
					except:
						logger.critical("gofigure.py: input error line "+str(line_count)+": TopGO input 5th column is not a number. This column needs to be the significantly enriched GO terms, but is not a number. has to be integer (no decimals). Please see 'example_data/topgo_example_input.tsv' for compatible TopGO input")
						exit()
					topgo_bool = True
				
			elif input_type.lower() == 'gostats':
				if not len(ssline) > 1:
					ssline = ' '.join(line.split()).split(' ')
				if len(ssline) == 7 and 'GO:' in line:
					if not ssline[1].startswith('GO:') and not len(ssline[1]) == 10:
						logger.critical("gofigure.py: input error line "+str(line_count)+": Column 2 is not a GO term.. Please see 'example_data/gostats_example_input.tsv' for compatible GOStats input")
						exit()
					try:
						float(ssline[2].replace('<',''))
					except:
						logger.critical("gofigure.py: input error line "+str(line_count)+": Column 3 is not a number. Needs to be a p-value. Please see 'example_data/gostats_example_input.tsv' for compatible GOStats input")
						exit()
					try:
						int(ssline[5].replace('<',''))
					except:
						logger.critical("gofigure.py: input error line "+str(line_count)+": Column 6 is not a number. Needs to be the number of significantly enriched GO terms. Has to be an integer (no decimals). Please see 'example_data/gostats_example_input.tsv' for compatible GOStats input")
						exit()
					gostats_bool = True
		
		line_count += 1

	## If TopGO or GOStats bools are still false, throw an error
	if input_type.lower() == 'topgo' and topgo_bool == False:
		logger.critical("gofigure.py: input error: No TopGO input found. Checked by looking for a row with 7 columns and 'GO:' in it. Pleas see 'example_data/topgo_example_input.tsv' for correct TopGO input")
		exit()
	if input_type.lower() == 'gostats' and gostats_bool == False:
		logger.critical("gofigure.py: input error: No GOStats input found. Checked by looking for a row with 7 columns and 'GO:' in it. Please see 'example_data/gostats_example_input.tsv' for compatible GOStats input")
		exit()

	return(input_file_path,output_dir_path,max_label,palette,colours,size,description_limit,legend,legend_position,max_pvalue,max_clusters,legend_columns,cluster_labels,similarity_cutoff,input_type,ontology,opacity,size_range,font_size,title,title_size,dpi,file_type,outfile_appendix,xlabel,ylabel,zlabel,priorities,logger,logger_std,random_state,sort_by,sum_user,top_level,name_changes)

## Function to process the input files. Either TopGO, GOStats, standard (GO term + P-value columns), or standard-plus (GO term + P-value + user defined value)
def process_input(input_file_path,input_type):
	infile_dict = defaultdict(list)
	infile = open(input_file_path,'r')
	infile.readline()
	go_set = set()
	line_count = 1
	if input_type.lower() == 'topgo':
		for line in infile.readlines():
			if not line.startswith('#') and not line.startswith('!') and not line.startswith('%') and not line.strip() == '' and 'GO:' in line:
				goid,term,description,annotated,user,expected,pvalue = line.strip().split('\t')
				if term in go_set:
					logger.error('Found duplicate GO term ID: '+term+' at line '+str(line_count)+'. Please fix duplicates before continuing.')
					exit()
				go_set.add(term)
				pvalue = pvalue.replace('<','').replace(',','.')
				infile_dict[term] = [term,user,pvalue]
			line_count += 1
	elif input_type.lower() == 'gostats':
		for line in infile.readlines():
			if not line.startswith('#') and not line.startswith('!') and not line.startswith('%') and not line.strip() == '' and 'GO:' in line:
				ssline = line.strip().split('\t')
				if not len(ssline) > 1:
					ssline = ' '.join(line.split()).split(' ')
				goid,term,pvalue,odds_ratio,expected,user,annotated = ssline
				if term in go_set:
					logger.error('Found duplicate GO term ID: '+term+' at line '+str(line_count)+'. Please fix duplicates before continuing.')
					exit()
				go_set.add(term)
				pvalue = pvalue.replace('<','').replace(',','.')
				infile_dict[term] = [term,user,pvalue]
			line_count += 1
	elif input_type.lower() == 'standard':
		for line in infile.readlines():
			if not line.startswith('#') and not line.startswith('!') and not line.startswith('%') and not line.strip() == '' and 'GO:' in line:
				term,pvalue = line.strip().split('\t')
				if term in go_set:
					logger.error('Found duplicate GO term ID: '+goid+' at line '+str(line_count)+'. Please fix duplicates before continuing.')
					exit()
				go_set.add(term)
				pvalue = pvalue.replace('<','').replace(',','.')
				infile_dict[term] = [term,pvalue]
			line_count += 1
	elif input_type.lower() == 'standard-plus':
		for line in infile.readlines():
			if not line.startswith('#') and not line.startswith('!') and not line.startswith('%') and not line.strip() == '' and 'GO:' in line:
				term,pvalue,user = line.strip().split('\t')
				if term in go_set:
					logger.error('Found duplicate GO term ID: '+term+' at line '+str(line_count)+'. Please fix duplicates before continuing.')
					exit()
				go_set.add(term)
				pvalue = pvalue.replace('<','').replace(',','.')
				infile_dict[term] = [term,user,pvalue]
			line_count += 1
	return(infile_dict)

## Instead of calculating the IC contents on the spot, read them from the IC file. This means less input file requirements and is slightly faster.  The file is created using scripts/ics.py.
def read_IC(scriptPath):
	ic_dict = defaultdict(float)
	frequency_dict = defaultdict(float)
	for line in open(scriptPath+'/data/ic.tsv'):
		if not line.startswith('#'):
			go,ic,frequency = line.strip().split('\t')
			ic_dict[go] = float(ic)
			frequency_dict[go] = float(frequency)
		elif line.startswith('# GOA version used: '):
			logger.info('GOA version used to calculate IC and frequencies: '+line.strip().split(': ')[1]+'\n')
	return(ic_dict,frequency_dict)

## Read the GO OBO file to match GO terms to their names and namespaces (description_dict and namespace_dict, respectively)
def read_OBO(scriptPath):
	description_dict = defaultdict(str)
	namespace_dict = defaultdict(str)
	obsolete_dict = defaultdict(set)
	alt_dict = defaultdict(str)
	for line in open(scriptPath+'/data/go.obo'):
		if line.startswith('id: GO:'):
			go = line.strip().split('id: ')[1]
			obsolete_bool = False
		elif line.startswith('alt_id: GO:'):
			alt_id = line.strip().split('alt_id: ')[1]
			alt_dict[alt_id] = go
		elif line.startswith('name: '):
			name = line.strip().split('name: ')[1]
			description_dict[go] = name
		elif line.startswith('namespace: '):
			namespace = line.strip().split('namespace: ')[1]
			namespace_dict[go] = namespace
		elif line.lower().startswith('is_obsolete: true'):
			obsolete_bool = True
			del description_dict[go]
			del namespace_dict[go]
		elif line.startswith('consider: GO:') and obsolete_bool == True:
			consider_go = line.strip().split('consider: ')[1]
			obsolete_dict[go].add(consider_go)
		elif line.startswith('data-version: '):
			obo_version = line.strip().split('data-version: ')[1]
			logger.info('go.obo version used: '+obo_version+'\n')
	return(description_dict,namespace_dict,obsolete_dict,alt_dict)

## Read child - parent relations, used for calculating semantic similarities
def read_parents_children(scriptPath):
	parents_dict = defaultdict(list)
	children_dict = defaultdict(list)
	for line in open(scriptPath+'/data/relations_full.tsv'):
		if not line.startswith('#'):
			go,parent = line.strip().split('\t')
			parents_dict[go].append(parent)
			children_dict[parent].append(go)
		elif line.startswith('#go.obo version used: '):
			logger.info('go.obo version used to create GO relations: '+line.strip().split(': ')[1]+'\n')
	return(parents_dict,children_dict)

## Calculates semantic similarities using Lin's formula
def calc_sem_sim(go1,go2,parents_dict,ic_dict):
	go1_parents = parents_dict[go1]
	go2_parents = parents_dict[go2]
	go1_IC = ic_dict[go1]
	go2_IC = ic_dict[go2]
	parent_IC = 0
	sem_sim = 0
	if go1 == go2: ## If both GO terms are the same, semantic similarity = 1
		sem_sim = 1
	elif go1 in go2_parents: ## If go1 is a parent of go2, the parent term used in the sem_sim equation is go1
		parent_IC = go1_IC
	elif go2 in go1_parents: ## If go2 is a parent of go1, the parent term used in the sem_sim equation is go2
		parent_IC = go2_IC
	else:
		parent_IC = 0
		closest_parent = ''
		for parent in go1_parents: ## If neither GO is a parent of the other, calculate common parent term with the highest IC content (meaning closest parent term)
			if parent in go2_parents:
				if ic_dict[parent] > parent_IC:
					parent_IC = ic_dict[parent]
					closest_parent = parent
	if not go1 == go2:
		if not go1_IC == 0 and not go2_IC == 0: ## If one of the GO terms has no IC, something is wrong. Check if the GO term is available in the data files
			sem_sim = (2*parent_IC)/(go1_IC+go2_IC) ## Calculate semantic similarity
	return(sem_sim)

## Creates 
def create_GO_dict(input_dict,input_type,namespace_input,namespace_dict,ic_dict,frequency_dict,max_pval,alt_dict):
	go_dict = defaultdict(list)
	conversion_bool = False
	if input_type == 'standard':
		for key,values in input_dict.items():
			go,pval = values
			pval = float(pval)
			if not go in namespace_dict:
				alt_dict_namespace = False
				if go in alt_dict:
					alt_go = go
					go = alt_dict[go]
					if namespace_dict[go] == namespace_input:
						if conversion_bool == False:
							logger.info("The following GO terms are converted or deleted because either they are an alt_id and not a primary id in the go.obo, because they are obsolete in the go.obo, or because they can't be found in the go.obo: \n\nOriginal GO\tNew GO(s)\tReason")
							conversion_bool = True
						logger.info(alt_go+'\t'+go+'\talt_id')
						alt_dict_namespace == True
				elif not go in alt_dict or alt_dict_namespace == False:
					if conversion_bool == False:
						logger.info("The following GO terms are converted or deleted because either they are an alt_id and not a primary id in the go.obo, because they are obsolete in the go.obo, or because they can't be found in the go.obo: \n\nOriginal GO\tNew GO(s)\tReason")
						conversion_bool = True
					logger.info('\t'+go+'\tremoved')
			if namespace_dict[go] == namespace_input and max_pval >= pval:
				if go in obsolete_dict:
					if conversion_bool == False:
						logger.info("The following GO terms are converted or deleted because either they are an alt_id and not a primary id in the go.obo, because they are obsolete in the go.obo, or because they can't be found in the go.obo: \n\nOriginal GO\tNew GO(s)\tReason")
						conversion_bool = True
					for new_go in obsolete_dict[go]:
						ic = float(ic_dict[new_go])
						frequency = float(frequency_dict[new_go])
						go_dict[new_go] = [pval,0,ic,frequency]
					logger.info(go+'\t'+str(obsolete_dict[go]).replace('{','').replace("'",'').replace('}','')+'\tObsolete')
				else:
					if namespace_dict[go] == namespace_input and max_pval >= float(pval):
						ic = float(ic_dict[go])
						frequency = float(frequency_dict[go])
						go_dict[go] = [pval,0,ic,frequency]
	else:
		for key,values in input_dict.items():
			go,user,pval = values
			pval = float(pval)
			if not go in namespace_dict:
				alt_dict_namespace = False
				if go in alt_dict:
					alt_go = go
					go = alt_dict[go]
					if namespace_dict[go] == namespace_input:
						if conversion_bool == False:
							logger.info("The following GO terms are converted or deleted because either they are an alt_id and not a primary id in the go.obo, because they are obsolete in the go.obo, or because they can't be found in the go.obo: \n\nOriginal GO\tNew GO(s)\tReason")
							conversion_bool = True
						logger.info(alt_go+'\t'+go+'\talt_id')
						alt_dict_namespace == True
				elif not go in alt_dict or alt_dict_namespace == False:
					if conversion_bool == False:
						logger.info("The following GO terms are converted or deleted because either they are an alt_id and not a primary id in the go.obo, because they are obsolete in the go.obo, or because they can't be found in the go.obo: \n\nOriginal GO\tNew GO(s)\tReason")
						conversion_bool = True
					logger.info(go+'\tNo replacement\tremoved')
			if namespace_dict[go] == namespace_input and max_pval >= pval:
				pval = float(pval)
				user = float(user)
				if go in obsolete_dict:
					if conversion_bool == False:
						logger.info("The following GO terms are converted or deleted because either they are an alt_id and not a primary id in the go.obo, because they are obsolete in the go.obo, or because they can't be found in the go.obo: \n\nOriginal GO\tNew GO(s)\tReason")
						conversion_bool = True
					for new_go in obsolete_dict[go]:
						ic = float(ic_dict[new_go])
						frequency = float(frequency_dict[new_go])
						go_dict[new_go] = [pval,user,ic,frequency]
					logger.info(go+'\t'+str(obsolete_dict[go]).replace('{','').replace("'",'').replace('}','')+'\tObsolete')
				elif not go in namespace_dict:
					alt_dict_namespace = False
					if go in alt_dict:
						alt_go = go
						go = alt_dict[go]
						if namespace_dict[go] == namespace_input:
							if conversion_bool == False:
								logger.info("The following GO terms are converted or deleted because either they are an alt_id and not a primary id in the go.obo, because they are obsolete in the go.obo, or because they can't be found in the go.obo: \n\nOriginal GO\tNew GO(s)\tReason")
								conversion_bool = True
							logger.info(alt_go+'\t'+go+'\talt_id')
							alt_dict_namespace == True
					elif not go in alt_dict or alt_dict_namespace == False:
						if conversion_bool == False:
							logger.info("The following GO terms are converted or deleted because either they are an alt_id and not a primary id in the go.obo, because they are obsolete in the go.obo, or because they can't be found in the go.obo: \n\nOriginal GO\tNew GO(s)\tReason")
							conversion_bool = True
						logger.info(go+'\tNo replacement\tremoved')
				else:
					ic = float(ic_dict[go])
					frequency = float(frequency_dict[go])
					go_dict[go] = [pval,user,ic,frequency]
	if 'GO:0008150' in go_dict:
		del go_dict['GO:0008150']
	if 'GO:0003674' in go_dict:
		del go_dict['GO:0003674']
	if 'GO:0005575' in go_dict:
		del go_dict['GO:0005575']
	return(go_dict)

## Creates clusters of GO terms, and chooses a representative term based on the workflow outlined below
def create_clusters(go_dict,parents_dict,children_dict,ic_dict,similarity_cutoff,priorities):
	go_dict = OrderedDict(go_dict.items())
	if priorities == None:
		priorities = []
	bool = True
	sem_sim_dict = defaultdict(list)
	prev_set = set()

	if top_level == None:
		top_level_terms = []
	else:
		top_level_terms = top_level.split(',')
	top_level_dict = defaultdict(float)
	for top_level_term in top_level_terms:
		if top_level_term in go_dict:
			pval,user,ic,frequency = go_dict[top_level_term]
			top_level_dict[top_level_term] = ic
	while bool == True:
		bool2 = True
		golist = list(OrderedDict(go_dict))
		for go in golist: ## Loop over every GO term
			if go in go_dict:
				pval,user,ic,frequency = go_dict[go]
				highest_sem_sim = 0
				highest_sem_sim_GO = ''
				highest_sem_sim_pval = 100
				highest_sem_sim_user = 0
				highest_sem_sim_frequency = 1
				highest_sem_sim_IC = 0
				highest_top_level_ic = 0
				most_informative_top_level_term = ''
				for top_level_term in top_level_terms:
					if go in children_dict[top_level_term]:
						top_level_ic = top_level_dict[top_level_term]
						if top_level_ic > highest_top_level_ic:
							highest_top_level_ic = top_level_dict[top_level_term]
							most_informative_top_level_term = top_level_term
				for go2 in golist:
					go2_top_level_term = go2
					if go in top_level_terms:
						go_top_level_bool = True
					else:
						go_top_level_bool = False
					go2_top_level_bool = False
					if go2 in top_level_terms:
						go2_top_level_bool = True
					else:
						for member in sem_sim_dict[go2]:
							if member[0] in top_level_terms:
								go2_top_level_bool = True
								go2_top_level_term = member[0]
					if not most_informative_top_level_term == '' and not go2_top_level_term == most_informative_top_level_term:
						continue
				#for go2,values2 in go_dict.items(): ## Loop over every GO term within the loop so they can be compared
					if not go2 == go and go2 in go_dict: ## Exclude the same GO term
							if (go_top_level_bool == False and go2_top_level_bool == False) or (go_top_level_bool == True and go2_top_level_bool == False and go2 in children_dict[go]) or (go_top_level_bool == False and go2_top_level_bool == True and go in children_dict[go2_top_level_term]):
								pval2,user2,ic2,frequency2 = go_dict[go2]
								sem_sim = calc_sem_sim(go,go2,parents_dict,ic_dict)
								if sem_sim > highest_sem_sim: ## Calculate the GO term which has the highest semantic similarity to the original GO term and adjust all values accordingly
									highest_sem_sim = sem_sim
									highest_sem_sim_GO = go2
									highest_sem_sim_pval = pval2
									highest_sem_sim_user = user2
									highest_sem_sim_frequency = frequency2
									highest_sem_sim_IC = ic2
				if highest_sem_sim >= similarity_cutoff: ## If the highest semantic similarity is higher than the cutoff proceed with grouping the GO terms and choosing the representative
					go_parents = parents_dict[go]
					highest_sem_sim_GO_parents = parents_dict[highest_sem_sim_GO]
					go1_bool = False
					go2_bool = False
					if go in priorities and highest_sem_sim_GO not in priorities:
						sem_sim_dict[go].append([go,pval,user,ic,frequency])
						sem_sim_dict[go].append([highest_sem_sim_GO,highest_sem_sim_pval,highest_sem_sim_user,highest_sem_sim_IC,highest_sem_sim_frequency])
						go2_bool = True
					elif highest_sem_sim_GO in priorities and not go in priorities:
						sem_sim_dict[highest_sem_sim_GO].append([go,pval,user,ic,frequency])
						sem_sim_dict[highest_sem_sim_GO].append([highest_sem_sim_GO,highest_sem_sim_pval,highest_sem_sim_user,highest_sem_sim_IC,highest_sem_sim_frequency])
						go1_bool = True
					elif frequency > 0.05 and not highest_sem_sim_frequency > 0.05: ## If the frequency of the original GO term in the entire GOA database is higher than 0.05, exclude this GO term as a potential representative term
						sem_sim_dict[highest_sem_sim_GO].append([go,pval,user,ic,frequency])
						sem_sim_dict[highest_sem_sim_GO].append([highest_sem_sim_GO,highest_sem_sim_pval,highest_sem_sim_user,highest_sem_sim_IC,highest_sem_sim_frequency])
						go1_bool = True
					elif highest_sem_sim_frequency > 0.05 and not frequency > 0.05: ## If the frequency of the second GO term in the entire GOA database is higher than 0.05, exclude this GO term as a potential representative term
						sem_sim_dict[go].append([go,pval,user,ic,frequency])
						sem_sim_dict[go].append([highest_sem_sim_GO,highest_sem_sim_pval,highest_sem_sim_user,highest_sem_sim_IC,highest_sem_sim_frequency])
						go2_bool = True
					elif highest_sem_sim_pval/pval < 0.5: ## If the p-value of the original term is significantly higher than the 2nd term, the 2nd term is the representative
						sem_sim_dict[highest_sem_sim_GO].append([go,pval,user,ic,frequency])
						sem_sim_dict[highest_sem_sim_GO].append([highest_sem_sim_GO,highest_sem_sim_pval,highest_sem_sim_user,highest_sem_sim_IC,highest_sem_sim_frequency])
						go1_bool = True
					elif pval/highest_sem_sim_pval < 0.5: ## If the p-value of the 2nd term is significantly higher than the original term, the original term is the representative
						sem_sim_dict[go].append([go,pval,user,ic,frequency])
						sem_sim_dict[go].append([highest_sem_sim_GO,highest_sem_sim_pval,highest_sem_sim_user,highest_sem_sim_IC,highest_sem_sim_frequency])
						go2_bool = True
					elif go in highest_sem_sim_GO_parents: ## If the original GO term is a parent of the 2nd GO term, the 2nd GO term is the representative
						sem_sim_dict[highest_sem_sim_GO].append([go,pval,user,ic,frequency])
						sem_sim_dict[highest_sem_sim_GO].append([highest_sem_sim_GO,highest_sem_sim_pval,highest_sem_sim_user,highest_sem_sim_IC,highest_sem_sim_frequency])
						go1_bool = True
					elif highest_sem_sim_GO in go_parents: ## If the 2nd GO term is a parent of the original GO term, the original GO term is the representative
						sem_sim_dict[go].append([go,pval,user,ic,frequency])
						sem_sim_dict[go].append([highest_sem_sim_GO,highest_sem_sim_pval,highest_sem_sim_user,highest_sem_sim_IC,highest_sem_sim_frequency])
						go2_bool = True
					else: ## If none of the above reaches to a representative, the original GO term is the representative
						sem_sim_dict[go].append([go,pval,user,ic,frequency])
						sem_sim_dict[go].append([highest_sem_sim_GO,highest_sem_sim_pval,highest_sem_sim_user,highest_sem_sim_IC,highest_sem_sim_frequency])
						go2_bool = True
					if go1_bool == True: ## If the original term is not chosen as a representative and already has an entry in the sem_sim_dict, add its values to the 2nd GO term in the sem_sim_dict and delete the old entry
						if go in sem_sim_dict:
							for val in sem_sim_dict[go]:
								if not val in sem_sim_dict[highest_sem_sim_GO]:
									sem_sim_dict[highest_sem_sim_GO].append(val)
							del sem_sim_dict[go]
						del go_dict[go]
						
					if go2_bool == True:
						if highest_sem_sim_GO in sem_sim_dict:
							for val in sem_sim_dict[highest_sem_sim_GO]:
								if not val in sem_sim_dict[go]:
									sem_sim_dict[go].append(val)
							del sem_sim_dict[highest_sem_sim_GO]
						del go_dict[highest_sem_sim_GO]
				
				else: ## If no go term is found with a semantic similarity >= cutoff, make a new entry in the sem_sim_dict
					sem_sim_dict[go].append([go,pval,user,ic,frequency])
		for go in golist: ## This loop makes sure the above workflow is looped until no changes are made, by checking if the go_dict keys are present as a key in sem_sim_dict. Any key not present in sem_sim_dict after every loop is deleted. Once all go_dict keys == sem_sim_dict keys, bool2 == true and therefore bool == false and while loop ends
			if not go in sem_sim_dict:
				bool2 = False
				if go in go_dict:
					del go_dict[go]
		if bool2 == True:
			bool = False
	return(sem_sim_dict)

## Creates a cluster dict containing just the bare information needed to create the dataframe for the scatterplot
def create_clusterdict(sem_sim_dict,description_dict):
	cluster_dict = defaultdict(list)
	for key,values in sem_sim_dict.items():
		dup_set = set()
		for value in values:
			go = value[0]
			user = value[2]
			if not go in dup_set:
				dup_set.add(go)
				name = description_dict[go]
				key_name = description_dict[key]
				cluster_dict[key].append([go,float(user)])
	return(cluster_dict)

## Creates the data frame for the scatterplot
def create_df(cluster_dict,go_dict,description_dict,colours,size,description_limit,legend,cluster_labels,random_state_input,semsim_dict):
	if not random_state_input == None:
		random_state_input = int(random_state_input)
	order = list(cluster_dict.keys()) ## Keep static order of df
	sem_array = np.zeros((len(order),len(order)))

	for i in range(0,len(order)): ## Create a sem_array containing semantic similarities between each representative of each cluster
		go1 = order[i]
		for k in range(0,len(order)):
			go2 = order[k]
			sem = calc_sem_sim(go1,go2,parents_dict,ic_dict)
			sem_array[i,k] = sem

	embedding = MDS(n_components=2,dissimilarity='euclidean',random_state=random_state_input) ## Two dimensional array based on semantic similarities between all representative GO terms
	if len(sem_array) == 0:
		return(pd.DataFrame())

	sem_array_transformed = embedding.fit_transform(sem_array)
	df = pd.DataFrame(data=sem_array_transformed, index=order, columns=['x','y']) ## Create the dataframe starting with the x and y semantic similarities

	names = df.index.values

	frequency_list = [] ## Frequency of GO in GOA database
	colour_criterium_list = [] ## Values which will be used to colour clusters
	pval_list = [] ## Pvalues
	description_list = [] ## Names of GO terms
	count_list = [] ## Counts of GO terms in cluster
	legend_list = [] ## DF holding the legend descriptions
	size_criterium_list = [] ## Values which will be used to size the clusters in the scatterplot
	go_list = [] ## List of GO terms part of the cluster
	count = 1 ## Index number
	user_defined_values_list = [] ## User defined value
	
	if name_changes == None:
		name_changes_list = []
	else:
		name_changes_list = name_changes.split(',')

	for go in names: ## Loop over each cluster and add values to above lists
		frequency = go_dict[go][2]
		if frequency == 0: ## If it somehow bugs out, make 0.001 so it wont throw error
			frequency = 0.001
		colours = colours.lower()
		colour_criterium = log10(float(go_dict[go][0]))
		pval = float(go_dict[go][0])
		if sum_user == False:
			user = go_dict[go][1]
		pval_list.append(pval)
		goterms = cluster_dict[go]
		tmpGOList = []
		new_go = ''
		new_description = ''
		for goterm in goterms: ## Add GO terms and user defined values
			tmpGOList.append(goterm[0])
			if goterm[0] in name_changes_list:
				new_go = goterm[0]
				new_description = description_dict[goterm[0]]		
		if sum_user == False:
			user = go_dict[go][1]
		elif sum_user == True:
			user = 0
			for goterm in goterms:
				user += go_dict[goterm[0]][1]
		go_list.append(tmpGOList)
		user_defined_values_list.append(user)
		if colours == 'pval': ## Decide which criteria to use for colouring of clusters
			colour_criterium = float(go_dict[go][0])
		elif colours == 'log10-pval':
			colour_criterium = log10(float(go_dict[go][0]))
		elif colours == 'user':
			colour_criterium = user
		elif colours == 'members':
			colour_criterium = 1+len(goterms)
		elif colours == 'frequency':
			colour_criterium = -log(frequency)
		if size == 'user': ## Decide which criteria to use for sizeing of clusters
			size_criterium = user
		elif size == 'pval':
			size_criterium = -log10(pval)
		elif size == 'frequency':
			size_criterium = log(frequency)
		elif size == 'members':
			size_criterium = log(1+len(goterms))
		else:
			size_criterium = size
		frequency_list.append(log(frequency))
		colour_criterium_list.append(colour_criterium)
		size_criterium_list.append(size_criterium)
		description = description_dict[go]
		if not new_description == '':
			description = new_description
		if not new_go == '':
			names = np.where(names == go, new_go, names)
			semsim_dict[new_go] = semsim_dict.pop(go)
			go = new_go
		if len(description) > description_limit: ## If a description exceeds the description character limit, truncate it
			description = description[0:description_limit]+'...'
		description_list.append(description)
		count_list.append(count)
		if legend.lower() == 'description': ## Decide how to fill the legend
			legend_list.append('. '+description)
		elif legend.lower() == 'go':
			legend_list.append('. '+go)
		elif legend.lower() == 'full':
			legend_list.append('. '+go+' '+description)
		elif legend.lower() == 'none':
			legend_list.append('. ')
		elif legend.lower() == 'exhaustive':
			legend_list.append('. '+go+' '+description)
		count += 1
	df['frequency'] = frequency_list
	df['colour'] = colour_criterium_list
	df['size'] = size_criterium_list
	df['pval'] = pval_list
	df['description'] = description_list
	df['legend'] = legend_list
	df['members'] = go_list
	df['representative'] = names
	df['user'] = user_defined_values_list
	if sort_by == 'user-descending':
		df = df.sort_values(by='user',ascending=False)	
	else:
		df = df.sort_values(by=sort_by) ## Sorts the DF based on P-value (lowest first)
	
	df['dotCount'] = count_list
	df = df.set_index('dotCount')
	for index, row in df.iterrows():
		if legend == 'exhaustive': ## Add member GO terms to exhaustive legend
			count = 0
			legend_description = ''
			for go in row['members']:
				description = description_dict[go]
				if len(description) > description_limit:
					description = description[0:description_limit]+'...'
				if count == 0:
					if cluster_labels.lower() == 'numbered':
						legend_description = str(index)+'. '+go+' '+description
					elif cluster_labels.lower() == 'go':
						legend_description = str(index)+'. '+go+' '+description
					count += 1
				else:
					legend_description = legend_description+'\n    '+go+' '+description
			df.loc[index,'legend'] = legend_description
		else: ## Else, fill legend normally
			if cluster_labels.lower() == 'numbered':
				df.loc[index,'legend'] = str(index)+row['legend']
			elif cluster_labels.lower() == 'go' or cluster_labels.lower() == 'go-arrows' or cluster_labels.lower() == 'description-numbered' or cluster_labels.lower() == 'description':
				df.loc[index,'legend'] = str(index)+row['legend']
	return(df,semsim_dict)

def scatterplot(df,output_dir_path,ontology,max_labels,palette,colour,size,legend,legend_position,max_clusters,legend_columns,cluster_labels,opacity,size_range,font_size,title,title_size,dpi,file_type,outfile_appendix,xlabel,ylabel,zlabel):
	cmap = ListedColormap(sns.color_palette(palette)) ## Colour palette for clusters
	df = df.head(max_clusters) ## If a maximum amount of clusters is given, truncate dataframe
	colour_criterium = df['colour']
	n = df.size
	if legend_columns == 'single': ## Determine legend column number
		colnumber = 1
	elif legend_columns == 'double':
		colnumber = 2
	elif legend_columns == 'triple':
		colnumber = 3
	if size_range == 'small': ## Determine size range of clusters
		size_range = (200,2000)
	elif size_range == 'medium':
		size_range = (400,4000)
	elif size_range == 'big':
		size_range = (600,6000)
	if str.isdigit(str(size)):  ## If size range is changed using fixed integer, use that
		size = int(size)*200
		size_range = (size,size)
	if colour.lower() == 'log10-pval' or colour.lower() == 'user' or colour.lower() == 'pval' or colour.lower() == 'members': ## This and below if-statement to colour clusters
		norm = plt.Normalize(colour_criterium.min(), colour_criterium.max())
		sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
		ax = sns.scatterplot(x='x',y='y',size='size',sizes=size_range,hue='colour',palette=cmap,edgecolor='black',data=df.iloc[::-1],alpha=0.5,linewidth=1)
	elif colour.lower() == 'unique':
		ax = sns.scatterplot(x='x',y='y',size='size',sizes=size_range,hue='dotCount',palette=cmap,edgecolor='black',data=df.iloc[::-1],alpha=0.5,linewidth=1)
	if legend_position.lower() == 'left': ## This and below if-statement to set legend position
		legend_position = 'right'
	elif legend_position.lower() == 'right':
		legend_position = 'left'
	empty_handle = matplotlib.patches.Rectangle((0,0), 1, 1, fill=False, edgecolor='none',visible=False) ## Below statements fix the legend
	handle_list = [empty_handle]*max_labels
	legend_ax = ax.legend(handle_list,df.legend.iloc[0:max_labels],bbox_to_anchor=(0.5, -0.15),loc='upper '+legend_position.lower(),handlelength=0, handletextpad=0, ncol=colnumber, frameon=False, fontsize=font_size)
	for item in legend_ax.legendHandles:
		item.set_visible(False)
	if legend == 'none': ## No legend if excluded in options
		ax.get_legend().remove()
	
	if cluster_labels == 'numbered': ## Couple of if-statements to determine how clusters are labeled
		if len(df.index) >= max_labels:
			for line in range(1,max_labels+1):
				ax.text(df.x[line], df.y[line], str(line), horizontalalignment='center', verticalalignment='center', size='small', color='black', weight='semibold', alpha=opacity)
		else:
			for line in range(1,len(df.index)+1):
				ax.text(df.x[line], df.y[line], str(line), horizontalalignment='center', verticalalignment='center', size='small', color='black', weight='semibold', alpha=opacity)
	elif cluster_labels == 'go' or cluster_labels == 'go-arrows' or cluster_labels == 'description' or cluster_labels == 'description-numbered':
		texts = []
		if cluster_labels == 'description-numbered':
			if len(df.index) >= max_labels:
				for line in range(1,max_labels+1):
					texts.append(ax.text(df.x[line], df.y[line], str(line)+'. '+df.description[line], horizontalalignment='center', verticalalignment='center', size='small', color='black', weight='semibold', alpha=opacity))
			else:
				for line in range(1,len(df.index)+1):
					texts.append(ax.text(df.x[line], df.y[line], str(line)+'. '+df.description[line], horizontalalignment='center', verticalalignment='center', size='small', color='black', weight='semibold',alpha=opacity))
		elif cluster_labels == 'description':
			if len(df.index) >= max_labels:
				for line in range(1,max_labels+1):
					texts.append(ax.text(df.x[line], df.y[line], df.description[line], horizontalalignment='center', verticalalignment='center', size='small', color='black', weight='semibold', alpha=opacity))
			else:
				for line in range(1,len(df.index)+1):
					texts.append(ax.text(df.x[line], df.y[line], df.description[line], horizontalalignment='center', verticalalignment='center', size='small', color='black', weight='semibold',alpha=opacity))	
		else:
			if len(df.index) >= max_labels:
				for line in range(1,max_labels+1):
					texts.append(ax.text(df.x[line], df.y[line], str(line)+'. '+df.representative[line], horizontalalignment='center', verticalalignment='center', size='small', color='black', weight='semibold', alpha=opacity))
			else:
				for line in range(1,len(df.index)+1):
					texts.append(ax.text(df.x[line], df.y[line], str(line)+'. '+df.representative[line], horizontalalignment='center', verticalalignment='center', size='small', color='black', weight='semibold',alpha=opacity))
		
		if cluster_labels == 'go-arrows' or cluster_labels == 'description' or cluster_labels == 'description-numbered':
			adjust_text(texts, arrowprops=dict(arrowstyle='->', color='red', alpha=0.5)) ## Python package to make sure labels overlap as little as possible
		else:
			adjust_text(texts)
	ax.set_xlabel(xlabel)
	ax.set_ylabel(ylabel)
	if not colour.lower() == 'unique': ## Change label on the right based on values used for colour determination
		kwargs = {'format': '%.0f'}
		colorbar = ax.figure.colorbar(sm,alpha=0.5,**kwargs)
		if zlabel == None:
			if colour.lower() == 'user':
				colorbar.set_label(zlabel)
				
			elif colour.lower() == 'members':
				colorbar.set_label('Unique GO terms per cluster')
			elif colour.lower() == 'frequency':
				colorbar.set_label('Frequency of GO term in the GOA database')
			elif colour.lower() == 'user':
				colorbar.set_label('Significantly enriched GO terms per cluster')
			elif colour.lower() == 'pval':
				colorbar.set_label('pvalue')
			elif colour.lower() == 'log10-pval':
				colorbar.set_label('log10 pvalue')
		else:
			colorbar.set_label(zlabel)
		colorbar.ax.tick_params(size=0)
	if not title == None:
		ax.set_title(title,fontsize=title_size)
	fig = ax.get_figure()
	if outfile_appendix == '':
		outfile_name = ontology+'_'+output_dir_path.replace('/','_')+'.'+file_type
	else:
		outfile_name = ontology+'_'+outfile_appendix+'.'+file_type
	fig.savefig(output_dir_path+'/'+outfile_name,dpi=dpi,bbox_extra_artists=(legend_ax,),bbox_inches='tight')
	plt.close()
	logger_std.info("Output for "+ontology.replace('_',' ')+"found in: "+output_dir_path+"/"+outfile_name)

def output_table(df,output_dir_path,ontology,appendix): ## Output raw table file based on the dataframe
	if appendix == '':
		df.to_csv(output_dir_path+'/'+ontology+'_'+output_dir_path.replace('/','_')+'.tsv',sep='\t')
	else:
		df.to_csv(output_dir_path+'/'+ontology+'_'+appendix+'.tsv',sep='\t')

def create_cluster_table(semsim_dict,description_dict,output_dir_path,ontology,appendix):
	if appendix == '':
		outfile = open(output_dir_path+'/'+ontology+'_full_table_'+output_dir_path.replace('/','_')+'.tsv','w')
	else:
		outfile = open(output_dir_path+'/'+ontology+'_full_table_'+appendix+'.tsv','w')
	outfile.write('Cluster representative\tCluster member\tCluster member description\tMember P-value\tMember user defined value\tMember IC\tMember frequency')
	for representative,members in semsim_dict.items():
		member_set = set()
		for member in members:
			go,pval,user_val,ic,frequency = member
			description = description_dict[go]
			if not go in member_set:
				outfile.write('\n'+representative+'\t'+go+'\t'+description+'\t'+"{:.3e}".format(pval)+'\t'+"{:.3e}".format(user_val)+'\t'+"{:.3e}".format(ic)+'\t'+"{:.3e}".format(frequency))
			member_set.add(go)
	outfile.close()

parser = read_arguments()
args = parser.parse_args()
outdir = args.output
if not os.path.isdir(outdir):
        os.mkdir(outdir)

input_file_path,output_dir_path,max_labels,palette,colours,size,description_limit,legend,legend_position,max_pvalue,max_clusters,legend_columns,cluster_labels,similarity_cutoff,input_type,ontology,opacity,size_range,font_size,title,title_size,dpi,file_type,outfile_appendix,xlabel,ylabel,zlabel,priorities,logger,logger_std,random_state,sort_by,sum_user,top_level,name_changes = argparse_extended_error_handling(parser)

logger.info('GO-Figure! version: 1.0.1\n')
logger.info('Date run: '+str(date.today())+'\n')
logger.info('Command line input: '+str(sys.argv)+'\n')

logger_std.info('Starting calculations')

input_dict = process_input(input_file_path,input_type)
ic_dict,frequency_dict = read_IC(scriptPath)
description_dict,namespace_dict,obsolete_dict,alt_dict = read_OBO(scriptPath)
parents_dict,children_dict = read_parents_children(scriptPath)

if ontology == 'all': ## Which ontologies to loop over
	ontology_list = ['biological_process','molecular_function','cellular_component']
elif ontology == 'bpo':
	ontology = 'biological_process'
	ontology_list = [ontology]
elif ontology == 'mfo':
	ontology = 'molecular_function'
	ontology_list = [ontology]
elif ontology == 'cco':
	ontology = 'cellular_component'
	ontology_list = [ontology]

for ontology in ontology_list:
	logger_std.info('\nCalculating '+ontology.replace('_',' '))
	logger.info('\nCalculating '+ontology.replace('_',' ')+'\n')
	go_dict = create_GO_dict(input_dict,input_type,ontology,namespace_dict,ic_dict,frequency_dict,max_pvalue,alt_dict)
	semsim_dict = create_clusters(go_dict,parents_dict,children_dict,ic_dict,similarity_cutoff,priorities)
	goset = set()
	for key,values in semsim_dict.items():
		goset.add(key)
		for value in values:
			goset.add(value[0])
	cluster_dict = create_clusterdict(semsim_dict,description_dict)
	if not len(semsim_dict) == 0:
		df,semsim_dict = create_df(cluster_dict,go_dict,description_dict,colours,size,description_limit,legend,cluster_labels,random_state,semsim_dict)
		scatterplot(df,output_dir_path,ontology,max_labels,palette,colours,size,legend,legend_position,max_clusters,legend_columns,cluster_labels,opacity,size_range,font_size,title,title_size,dpi,file_type,outfile_appendix,xlabel,ylabel,zlabel)
		output_table(df,output_dir_path,ontology,outfile_appendix)
		create_cluster_table(semsim_dict,description_dict,output_dir_path,ontology,outfile_appendix)
	else:
		logger_std.info('No GO terms for ontology found: '+ontology)
		logger.info('No GO terms for ontology found: '+ontology+'\n')

logger_std.info("Finished!")
