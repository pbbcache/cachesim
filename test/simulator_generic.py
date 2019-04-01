#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
# simulator_generic.py
#
##############################################################################
#
# Copyright (c) 2019 Juan Carlos Saez<jcsaezal@ucm.es>, Jorge Casas <jorcasas@ucm.es>, Adrian Garcia Garcia<adriagar@ucm.es>
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
from simulator_core import *
from simulator_results import *
from common_charts import *
import time
import sys
import argparse

## Turn list of lists [key=val] into dict
## For now it only recognizes None, Strings, Int and Booleans as values
def parse_option_list(optlist):
	optdict={}

	if not optlist:
		return optdict

	for opt in optlist:
		item=opt[0]
		tok=item.split("=")
		if len(tok)==2:
			key=tok[0]
			val=tok[1]

			if val=="False":
				val=False
			elif val=="True":
				val=True
			elif val=="None":
				val=None
			else:
				## Try int (otherwise it will be a string)
				try:
					val=int(tok[1])
				except ValueError:
					val=tok[1]
		else:
			## Assume bool
			key=item
			val=True
		## Add to dict
		optdict[key]=val
	return optdict



def parse_range(str_range,nr_workloads):
	if not str_range or str_range=="-" or str_range=="*":
		return [i for i in range(nr_workloads)]
	else:
		## Available formats
		## 1-
		## 1-3
		## -4
		## 1,2,3-5
		subranges=str_range.split(",")
		items=[]
		for subrange in subranges:
			if  subrange[-1]=="-":
				left_limit=int(subrange[0:len(subrange)-1])
				items.extend([i for i in range(left_limit-1,nr_workloads)])
			elif subrange[0]=="-":
				right_limit=int(subrange[1:len(subrange)])
				items.extend([i for i in range(0,right_limit)])			
			elif "-" in subrange:
				limits=subrange.split("-")
				left_limit=int(limits[0])
				right_limit=int(limits[1])

				assert left_limit>=1 and left_limit<=nr_workloads
				assert right_limit>=1 and right_limit<=nr_workloads
				items.extend([i for i in range(left_limit-1,right_limit)])
			else: ## Number
				num=int(subrange)
				assert num>=1 and num<=nr_workloads
				items.append(num-1)

		return items
## Prepare parser
parser = argparse.ArgumentParser(description='Test main file for simulator (UCP, Whirlpool, etc.)')
parser.add_argument("-s","--summary-file",default="../data/volta2_turbo_off.300W.csv",help='File with the data collected offline for every benchmark')
parser.add_argument("-b","--max-bandwidth",type=float,default=float('Inf'),help='Max bandwidth observed for the target machine')
parser.add_argument("-p","--parallel",action='store_true',help='Enable parallel search for the optimal')
parser.add_argument("-H","--harness",action='store_true',help='Use harness workload file as input')
parser.add_argument("-C","--generate-chart",action='store_true',help='Enable the STP-vs-UNF chart generator')
parser.add_argument("-m","--bw-model",default="simple",help='Select BW model to be applied (%s)' % str(glb_available_bw_models))
parser.add_argument("-a","--algorithms",default="ucp",help='Comma-separated algorithms (%s)' % str(glb_available_bw_models))
parser.add_argument("-f","--format",default="full",help='Format style for output')
parser.add_argument("-d","--debugging",action='store_true',help='Enable debugging mode (assume normal input files) and keep workload 0')
parser.add_argument("-r","--use-range",default="-",help='Pick selected workloads only by specifying a range')
parser.add_argument("-O","--option",action='append',nargs=1,metavar="key=val", help='Use generic algorithm-specific options')
parser.add_argument('wfile', metavar='WORKLOAD_FILE', help='a workload file')	

args=parser.parse_args(sys.argv[1:])
optdict=parse_option_list(args.option)

select_bw_model(args.bw_model)
## Add bw_model to global dict
optdict["bw_model"]=args.bw_model

max_bandwidth=args.max_bandwidth

if args.harness:
	benchs=parse_harness_file(args.wfile)
	(workloads,nr_ways)=get_workloads_table_from_list_gen(args.summary_file,[benchs],separator=',')

else:
	## For testing with normal input file
	(workloads,nr_ways)=get_workloads_table_from_csv_gen(args.summary_file,args.wfile,separator=',')


workloads_name=args.wfile.split("/")[-1].rsplit(".", 1)[0]

## Schedulers 
algorithms=args.algorithms.split(",")


## get workload range
workload_range=parse_range(args.use_range,len(workloads))

## Structure to generate charts
data = {}
markerSpecs= {}
markers_available=['.','^','s','P','*','H','x']
nr_markers=len(markers_available)

## Verify and prepare data for charts...
for idx_alg,algorithm in enumerate(algorithms):
	if  not algorithm in glb_algorithms_supported:
		print >> sys.stderr, "No such partitioning algorithm: ", algorithm
		exit(1)
	
	## Prepare data and markers
	if args.generate_chart:
		data[algorithm]=[]
		markerSpecs[algorithm]=markers_available[idx_alg % nr_markers]

show_header=True
for idx in workload_range:
	workload=workloads[idx]
	i=idx+1
	wname="W%d" % i
	for alg_idx,algorithm in enumerate(algorithms):
		if ("print_times" in optdict) and optdict["print_times"]:
			print >> sys.stderr, "********* Workload %s - Algorithm %s ***********" %  (wname,algorithm)
		sol_data=apply_part_algorithm(algorithm,workload,nr_ways,max_bandwidth,parallel=args.parallel,uoptions=optdict,workload_name=wname)
		if args.format=="harness":
			sim_print_sol_masks(algorithm,i,sol_data,max_bandwidth)
		elif args.format=="harness-debussy":
			sim_print_sol_masks_debussy(algorithm,i,sol_data,max_bandwidth)
		elif args.format=="quiet":
			pass
		elif args.format=="table":
			sim_print_sol_table(algorithm,i,sol_data,max_bandwidth,print_header=show_header)
		else:
			sim_print_sol_simple(algorithm,i,sol_data,max_bandwidth,print_header=(alg_idx==0))

		show_header=False
		
		## Prepare stuff
		if args.generate_chart:
			## REtrieve metrics
			metrics=compute_basic_metrics(sol_data,max_bandwidth)
			data[algorithm].append(LabeledPoint(wname, metrics["stp"], metrics["unfairness"]))


if args.generate_chart:
	fsize=16
	rcParams['text.usetex'] = True #Let TeX do the typsetting
	rcParams['text.latex.preamble'] = [r'\usepackage{sansmath}', r'\sansmath'] #Force sans-serif math mode (for axes labels)
	rcParams['font.family'] = 'sans-serif' # ... for regular text
	rcParams['font.sans-serif'] = 'Helvetica' #'Helvetica, Avant Garde, Computer Modern Sans serif'
	rcParams['xtick.labelsize'] = fsize
	rcParams['ytick.labelsize'] = fsize
	rcParams['legend.fontsize'] = fsize
	rcParams['grid.linewidth']= 1.0
	figSize=[9.0,8.0]
	print "Generating chart as",workloads_name+".pdf ...",
	scatter2DGen(data,markerSpecs,"STP","Unfairness",0,figSize,filename=workloads_name+".pdf",windowTitle=workloads_name,legendLocation='upper left',axes_labelsize=fsize)
	print "Done"
	#Uncomment to display the chart windows
	#plt.show()

