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
from __future__ import print_function
from simulator_core import *
import simulator_core
from simulator_results import *
#from common_charts import *
import time
import sys
import argparse
import ipyparallel as ipp
import pandas as pd

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

def get_user_defined_assignment(user_file):
	ufd = open(user_file,"r")
	assignment = ufd.read()

	workloads_assignment=[]
	for line in assignment.split("\n"):
		if ";" in line:
			cluster_ids,cluster_ways = line.split(";")
			# Cluster id stores to which clust each app belongs
			cluster_id = [int(clid) for clid in cluster_ids.split(",")]
			cluster_way = [int(clw) for clw in cluster_ways.split(",")]

			per_clust_masks = get_partition_masks(cluster_way)
			per_app_masks = [None]*len(cluster_id)
			per_app_ways = [None]*len(cluster_id)
			for i in range(len(cluster_id)):
				per_app_ways[i] = cluster_way[cluster_id[i]]
				per_app_masks[i] = per_clust_masks[cluster_id[i]]

			workloads_assignment.append((per_app_ways,per_app_masks,cluster_id))

	ufd.close()
	return workloads_assignment

def get_user_defined_assignment(user_file):
	ufd = open(user_file,"r")
	assignment = ufd.read()

	workloads_assignment=[]
	for line in assignment.split("\n"):
		if ";" in line:
			cluster_ids,cluster_ways = line.split(";")
			# Cluster id stores to which clust each app belongs
			cluster_id = [int(clid) for clid in cluster_ids.split(",")]
			cluster_way = [int(clw) for clw in cluster_ways.split(",")]

			per_clust_masks = get_partition_masks(cluster_way)
			per_app_masks = [None]*len(cluster_id)
			per_app_ways = [None]*len(cluster_id)
			for i in range(len(cluster_id)):
				per_app_ways[i] = cluster_way[cluster_id[i]]
				per_app_masks[i] = per_clust_masks[cluster_id[i]]

			workloads_assignment.append((per_app_ways,per_app_masks,cluster_id))

	ufd.close()
	return workloads_assignment

def print_mapping_app(df,w,alg_name,mode):
	subdf=df[(df["W#"]==w) & (df["Algorithm"]==alg_name)]
	if len(subdf)==0:
		print("Algorithm %s not found in file:" % alg_name)
		exit(1)	
  
	benchids=subdf["BenchID"].unique().tolist()
	ways={}
	masks={}
	llc_idv={}
	cluster_llc=[]
	group_count={}
	cluster_ids=[]
	mapping_spec=[]
	app_masks=[]
	app_names=[]
	cluster_zero=False
	old=False
	
	for benchid in benchids:
		if (mode==1 and benchid==len(benchids)) or (mode==2 and benchid==0):
			continue
		dfs=subdf[subdf["BenchID"]==benchid]
		app=dfs["Name"].values[0]
		if mode==1:
			raw=(dfs["Mask/NR_WAYS"].values[0])
			way_mask=raw.split("(")[0]
			nr_ways=int(raw.split("(")[1][:-1])
			#print(raw,"-->",nr_ways)
			way_mask=(dfs["Mask/NR_WAYS"].values[0]).split("(")[0]
			llc_id=int(dfs["LLCID"].values[0])
			cluster_id=int(dfs["Cluster"].values[0])
		else:
			way_mask=dfs[dfs["Property"]=="way_mask"].values[0]
			nr_ways=int(dfs[dfs["Property"]=="way_count"].values[0])
			llc_id=int(dfs[dfs["Property"]=="llc_id"].values[0])
			cluster_id=int(dfs[dfs["Property"]=="cluster_id"].values[0])

		if llc_id in group_count:
			count=group_count[llc_id]
			group_count[llc_id]=count + 1 
		else:
			count=0
			group_count[llc_id]=1	

		mapping_spec.append("G%d/%d" % (llc_id,count))
   
		## For debugging
		#print(app,llc_id,cluster_id,way_mask)
		cid=cluster_id
		ways[cid]=nr_ways
		masks[cid]=way_mask
		llc_idv[cid]=llc_id
		cluster_ids.append(cid)
		cluster_llc.append(llc_id)
		app_masks.append(way_mask)
		app_names.append(app)

		if cid==0:
			cluster_zero=True

	if cluster_zero:
		start=0
		end=len(masks)
	else:
		start=1
		end=len(masks)+1
  
	if old:
		## Combine pairs idx_cluster-app/llc_group;
		cluster_spec=list([ "%s/%d" % (masks[i],llc_idv[i])  for i in range (start,end)])		
		cluster_spec_str=','.join(cluster_spec)
		id_clusters=','.join(map(str, cluster_ids))
		id_groups=','.join(map(str, cluster_llc))
		app_masks_str=','.join(app_masks)
		print("%s;%s;%s;%s" % (app_masks_str,id_groups,id_clusters,cluster_spec_str))			
	else:
		cluster_per_app=','.join(map(str, cluster_ids))

		cluster_spec=list([ "%d/%s/%d" % (i,masks[i],llc_idv[i])  for i in range (start,end)])		
		cluster_spec_str=','.join(cluster_spec)
		mapping_spec_str=','.join(mapping_spec)
		app_names_str=','.join(app_names)
		print("%s;%s;%s;%s" % (app_names_str,cluster_per_app,cluster_spec_str,mapping_spec_str))			

 
def parse_opt_switch(opt_spec):
	correspondence={"max": True, "min": False, "maximize": True , "minimize": False}
	aliases={"unf": "unf_stp","gmipc":"gmean_ipc","aipc":"aggregate_ipc","hspeed":"hmean_speedup","m1":"m1_metric","STP":"stp","cov":"cov_unfairness_metric","jain":"jain_fairness","antt":"antt_metric" }
	tokens=opt_spec.split(":")

	if len(tokens)!=3 or tokens[2] not in correspondence.keys():
		print("Wrong format of argument in opt option")
		exit(1)

	## Apply alias
	if tokens[1] in aliases:
		funname=aliases[tokens[1]]
	else:
		funname=tokens[1]

	f=getattr(simulator_core,funname,None)
	if f is None:
		print("Could not find function name %s in simulator_core" % tokens[1])
		exit(1)		

	return (f,correspondence[tokens[2]],"%s-%s-%s" % (tokens[0],tokens[1],tokens[2]))

def parse_partitioning_algorithms(alg_spec):
	algorithms=alg_spec.split(",")
	alg_spec=[]
	for alg in algorithms:
		index=alg.find(":")
		if index==-1:
			alg_spec.append((alg,[unfairness,False,alg]))
		else:
			alg_spec.append((alg[0:index],parse_opt_switch(alg)))
	return alg_spec


def list_algorithms():
	names=get_algorithm_names()
	print("**Partitioning Algorithms**")
	for name in names:
		print(name)


## Prepare parser
parser = argparse.ArgumentParser(description='Test main file for simulator (UCP, Whirlpool, etc.)')
parser.add_argument("-s","--summary-file",default="../data/volta2_turbo_off.300W.csv",help='File with the data collected offline for every benchmark')
parser.add_argument("-b","--max-bandwidth",type=float,default=float('Inf'),help='Max bandwidth observed for the target machine')
parser.add_argument("-p","--parallel",action='store_true',help='Enable parallel search for the optimal')
parser.add_argument("-P","--parsim",action='store_true',help='Launch all simulations in parallel')
parser.add_argument("-H","--harness",action='store_true',help='Use harness workload file as input')
parser.add_argument("-C","--generate-chart",action='store_true',help='Enable the STP-vs-UNF chart generator')
parser.add_argument("-W","--show-window",action='store_true',help='Show matplotlib window with chart (must be used along with -C)')
parser.add_argument("-L","--list-algorithms",action='store_true',help='List algorithms implemented in the simulator.')
parser.add_argument("-m","--bw-model",default="simple",help='Select BW model to be applied (%s)' % str(glb_available_bw_models))
parser.add_argument("-a","--algorithms",default="ucp",help='Comma-separated algorithms (%s)' % str(glb_algorithms_supported))
parser.add_argument("-f","--format",default="table",help='Format style for output')
parser.add_argument("-d","--debugging",action='store_true',help='Enable debugging mode (assume normal input files) and keep workload 0')
parser.add_argument("-r","--use-range",default="-",help='Pick selected workloads only by specifying a range')
parser.add_argument("-w","--force-ways",type=int,default=0,help='Force number of ways of the cache (it must be smaller or equal than the maximum number of ways inferred from input file)')
parser.add_argument("-t","--topology",default="uma",help='Select topology to be considered (%s)' % str(glb_available_topologies))
parser.add_argument("-O","--option",action='append',nargs=1,metavar="key=val", help='Use generic algorithm-specific options')
parser.add_argument('wfile', metavar='WORKLOAD_FILE', help='a workload file')
parser.add_argument("-o","--offline",default='-', help='Use simulator output file to generate data for static mapping')

## Special case for listing algorithms (add a dummy argument if needed)
if len(sys.argv)==2 and sys.argv[1]=="-L":
	sys.argv.append("foo")

args=parser.parse_args(sys.argv[1:])
initialize_simulator(args)
optdict=parse_option_list(args.option)

## Add bw_model and topology to global dict
optdict["bw_model"]=args.bw_model
optdict["topology"]=args.topology
if "user" in args.algorithms:
	if "user_file" in optdict:
		optdict["user_assignment"] = get_user_defined_assignment(optdict["user_file"])
	else:
		print("Selected user defined schedule but no file provided with (e.g., -O user_file=data/user_assignment.csv)")
		exit(1)

max_bandwidth=args.max_bandwidth

if  args.wfile=="-L" or args.list_algorithms:
	list_algorithms()
	exit(0)	


if args.harness:
	benchs=parse_harness_file(args.wfile)
	(workloads,nr_ways)=get_workloads_table_from_list_gen(args.summary_file,[benchs],separator=',')

else:
	## For testing with normal input file
	(workloads,nr_ways)=get_workloads_table_from_csv_gen(args.summary_file,args.wfile,separator=',')

if args.force_ways>0 and args.force_ways<=2*nr_ways:
	nr_ways=args.force_ways

workloads_name=args.wfile.split("/")[-1].rsplit(".", 1)[0]

## Parse algorithms
alg_spec=parse_partitioning_algorithms(args.algorithms)

## get workload range
workload_range=parse_range(args.use_range,len(workloads))

## Process offline mode here
if args.offline != "-":
	try:
		df=pd.read_csv(args.offline, delim_whitespace=True)
		## Detect format
		if "Algorithm" in df.columns:
			mode=1
		else:
			mode=2
			df.rename(columns={"workload":"W#","scheme":"Algorithm","id":"BenchID","app":"Name","property":"Property","value":"Value"},inplace=True)
	except:
		print("Could not open or process simulator's log file:", args.offline)
		exit(1)
	## Go ahead
	for idx in workload_range:
		workload=workloads[idx]
		i=idx+1
		wname="W%d" % i
		for alg_idx,alg in enumerate(alg_spec):
			algorithm=alg[0]
			alg_name=alg[1][2]
			opt_spec=alg[1]
			print_mapping_app(df,wname,algorithm,mode)

	exit(0)

## Structure to generate charts
data = {}
markerSpecs= {}
markers_available=['.','^','s','P','*','H','x']
nr_markers=len(markers_available)

## Verify and prepare data for charts...
for idx_alg,alg in enumerate(alg_spec):
	algorithm=alg[0]
	alg_name=alg[1][2]
	if not algorithm in glb_algorithms_supported:
		print("No such partitioning algorithm: ", algorithm, file=sys.stderr)
		exit(1)
	
	## Prepare data and markers
	if args.generate_chart:
		data[alg_name]=[]
		markerSpecs[alg_name]=markers_available[idx_alg % nr_markers]


## Launch simulations in parallel if the user requested it 
if args.parsim:
	results=launch_simulations_in_parallel(workloads,alg_spec,workload_range,nr_ways,max_bandwidth, optdict,args.bw_model)


show_header=True
for idx in workload_range:
	workload=workloads[idx]
	i=idx+1
	wname="W%d" % i
	for alg_idx,alg in enumerate(alg_spec):
		algorithm=alg[0]
		alg_name=alg[1][2]
		opt_spec=alg[1]
		if ("print_times" in optdict) and optdict["print_times"]:
			print("********* Workload %s - Algorithm %s ***********" %  (wname,alg_name), file=sys.stderr)
		if  args.parsim:
			sol_data=results[(idx,alg_idx)]
		else:
			sol_data=apply_part_algorithm(algorithm,workload,nr_ways,max_bandwidth,parallel=args.parallel,debugging=args.debugging,uoptions=optdict,workload_name=wname,opt_spec=opt_spec)
		if args.format=="harness":
			sim_print_sol_masks(alg_name,i,sol_data,max_bandwidth)
		elif args.format=="harness-debussy":
			sim_print_sol_masks_debussy(alg_name,i,sol_data,max_bandwidth)
		elif args.format=="quiet":
			pass
		elif args.format=="table":
			sim_print_sol_table(alg_name,i,sol_data,max_bandwidth,print_header=show_header,user_options=optdict)
		elif args.format=="dataframe" or args.format=="df":
			sim_print_sol_dataframe(alg_name,i,sol_data,max_bandwidth,print_header=show_header,user_options=optdict)			
		elif args.format=="cluster":
			sim_print_cluster_info(alg_name,i,sol_data,max_bandwidth,print_masks=False)
		elif args.format=="harness-cluster":
			sim_print_cluster_info(alg_name,i,sol_data,max_bandwidth,print_masks=True)
		else:
			sim_print_sol_simple(alg_name,i,sol_data,max_bandwidth,print_header=(alg_idx==0))

		show_header=False
		
		## Prepare stuff
		if args.generate_chart:
			## Retrieve metrics
			metrics=compute_basic_metrics(sol_data,max_bandwidth)
			data[alg_name].append(LabeledPoint(wname, metrics["stp"], metrics["unfairness"]))


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
	print("Generating chart as",workloads_name+".pdf ...",
	scatter2DGen(data,markerSpecs,"STP","Unfairness",0,figSize,filename=workloads_name+".pdf",windowTitle=workloads_name,legendLocation='upper left',axes_labelsize=fsize,show_window=args.show_window))
	print("Done")
	#Uncomment to display the chart windows
	#plt.show()

