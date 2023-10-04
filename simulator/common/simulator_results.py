# -*- coding: utf-8 -*-

#
# simulator_results.py
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
from opt_clustering import *
## Must be commented out in public version
from simulator_exploration import *
import numpy as np
import pandas as pd

## Efficient function to determine all metrics from slowdown vector...
## Returns a dictionary with the collected metrics....
def compute_basic_metrics(sol_data,max_bandwidth=float('Inf')):
	metrics={}
	(sol_spec,statistics)=sol_data
	(workload,per_app_ways,per_app_masks,cluster_id)=sol_spec
	slowdown_vector=get_slowdown_vector(workload,per_app_ways,max_bandwidth)
	slowdown_vector_np=np.array(slowdown_vector, dtype='float')
	slowdown_vector_nobw=get_slowdown_vector(workload,per_app_ways)
	nr_apps=len(slowdown_vector)
	predicted_ipc_vector=[]
	ipc_alone_vector=[]
 
	metrics["unfairness"] = max(slowdown_vector) / min(slowdown_vector)
	metrics["antt"] =sum(slowdown_vector) / nr_apps

	stp=0.0
	gmean_ipc=1.0
	gmean_speedup=1.0
	m1 = 0
	hmean_speedup = 0
	aggregate_ipc=0

	for i,slowdown in enumerate(slowdown_vector):
		predicted_ipc=workload[i].get_ipc_alone()/slowdown    
		stp += 1.0/slowdown
		gmean_ipc *= predicted_ipc
		gmean_speedup *= 1.0/slowdown

		aggregate_ipc += predicted_ipc;
		predicted_ipc_vector.append(predicted_ipc)
		ipc_alone_vector.append(workload[i].get_ipc_alone())
		hmean_speedup += slowdown ## slowdown is the inverse of speedup

		for j, slowdown2 in enumerate(slowdown_vector):
			if (j!=i):
				m1 += abs(slowdown2 - slowdown)

	gmean_ipc=round(gmean_ipc**(1.0/nr_apps),3)
	gmean_speedup=round(gmean_speedup**(1.0/nr_apps),3)
 
	metrics["stp"] = round(stp,3) 
	metrics["slowdown"]= slowdown_vector
	metrics["slowdown_no_bw"]= slowdown_vector_nobw
	metrics["gmean_speedup"] = gmean_speedup
	metrics["gmean_ipc"]= gmean_ipc
	cov=np.std(slowdown_vector)/np.mean(slowdown_vector)
	metrics["unfairness-cov"]=cov
	metrics["jain-fairness"]=1.0/(1.0+cov**2)
	metrics["predicted_ipc"]=predicted_ipc_vector
	metrics["ipc_alone"]=ipc_alone_vector
	metrics["m1_metric"] = round(m1, 3)
	metrics["hmean_speedup"] = len(slowdown_vector)/hmean_speedup
	metrics["aggregate_ipc"] = aggregate_ipc
 
	return metrics



def mirror_mask(mask,nr_ways):
	low_bit=0
	count=0
	found=False
	hex_mask=int(mask,16)
	half_idx=nr_ways/2

## Find first bit
	for i in range(nr_ways):
		is_one=((hex_mask & 1<<i)!=0)
		if not found and is_one:
			low_bit=i
			count=1
			found=1
		elif found:
			if is_one:
				count+=1
			else:
				break ## Nothing else to parse

	high_bit=low_bit+count-1


	if high_bit<=nr_ways/2: ## Beyond half (move right)
		new_high=nr_ways-1-low_bit
		new_low=new_high-count+1
	else:
		new_low=nr_ways-1-high_bit

	return hex(((1<<count)-1)<<new_low)


def fix_intel_bug_aux(ways,masks,cluster_ids):
	## Calculate nr_ways
	nr_ways=0
	available_masks=masks.values()

	for way_count in ways.values():
		nr_ways+=way_count

	## Build conflicting mask
	conflicting_mask=hex(1<<(nr_ways-1))
	first_mask=hex(1<<0)

	## Mask found
	if (first_mask not in available_masks) and (conflicting_mask in available_masks):

		## Specular way
		return [mirror_mask(mask,nr_ways) for mask in available_masks]

	else:
		## Nothing to do
		return masks
		

def fix_intel_bug(per_app_ways,per_app_masks,cluster_id):
	per_app_ways = list(map(int, per_app_ways))
	## Calculate way count 
	ways={}

	for idx in range(len(per_app_ways)):
		cid=cluster_id[idx]
		ways[cid]=per_app_ways[idx]

	## Calculate nr_ways
	nr_ways=0

	for way_count in ways.values():
		nr_ways+=way_count

	## Build conflicting mask
	conflicting_mask=hex(1<<(nr_ways-1))
	first_mask=hex(1<<0)

	## Mask found
	if (first_mask not in per_app_masks) and (conflicting_mask in per_app_masks):
		## Replace value
		for idx,mask in enumerate(per_app_masks):
			new_mask=mirror_mask(mask,nr_ways)
			## Patch per_app_masks directly
			per_app_masks[idx]=new_mask

# Function to pretty-print a workload
def print_workload(workload,id):
	print("****WORKLOAD #%d****" % id)
	bnames=[]
	for benchmark in workload:
		bnames.append(benchmark.name)
	print(bnames)

def sim_print_sol_simple(sched_title,idx,sol_data,max_bandwidth=float('Inf'),print_header=False):
	metrics=compute_basic_metrics(sol_data,max_bandwidth)
	(sol_spec,statistics)=sol_data
	(workload,per_app_ways,per_app_masks,cluster_id)=sol_spec

	## Determine slowdown with and widthout BW
	sbw=metrics["slowdown"]
	snbw=metrics["slowdown_no_bw"]

	if print_header:
		print_workload(workload,idx)

	print("==== Solution for %s ===" % sched_title)
	i=0
	for benchmark in workload:
		print("%s --> (%d ways,%f,%f)" % (benchmark.name,per_app_ways[i],sbw[i],snbw[i]))
		i=i+1
	print("--------------")
	print("--------------")
	print("STP=%f" % metrics["stp"])
	print("Unf=%f" % metrics["unfairness"])
	print("--------------")
	print("=======================")

def sim_print_sol_table(sched_title,idx,sol_data,max_bandwidth=float('Inf'),print_header=False,**kwargs):

	kwargs.setdefault("user_options",None)
	kwargs.setdefault("workload_prefix","W")
	kwargs.setdefault("alg_suffix","")
	kwargs.setdefault("use_csv",False)
	kwargs.setdefault("mapping_output",False)

	if "user_options" in kwargs:
		## Merge dict with user options... [Overriding defaults...]
		kwargs.update(kwargs["user_options"])

	workload_prefix=kwargs["workload_prefix"]
	alg_suffix=kwargs["alg_suffix"]
	use_csv=kwargs["use_csv"]
	mapping_output=kwargs["mapping_output"]

	metrics=compute_basic_metrics(sol_data,max_bandwidth)
	(sol_spec,statistics)=sol_data
	(workload,per_app_ways,per_app_masks,cluster_id)=sol_spec
	fix_intel_bug(per_app_ways,per_app_masks,cluster_id)


	## Determine slowdown with and widthout BW
	sbw=metrics["slowdown"]
	snbw=metrics["slowdown_no_bw"]

	if use_csv:
		if mapping_output:
			formatted_string="%s,%s,%s,%s,%s,%s,%s,%s,%s,%s"
		else:
			formatted_string="%s,%s,%s,%s,%s,%s,%s,%s,%s"
	else:
		if mapping_output:
			formatted_string="%-Xs %-Ys %-9s %-17s %-17s %-13s %-15s %-17s %-17s %-17s"
		else:
			formatted_string="%-Xs %-Ys %-9s %-17s %-13s %-15s %-17s %-17s %-17s"
		formatted_string=formatted_string.replace("X",str(len(workload_prefix)+3)).replace("Y",str(len(alg_suffix)+19))

	workload_name="%s%i" % (workload_prefix,idx)
	alg_name="%s%s" % (sched_title,alg_suffix)

	if print_header:
		if  mapping_output:
			print(formatted_string % ("W#","Algorithm","BenchID","Name","Cluster","LLCID","Mask/NR_WAYS","SlowdownNB/ANTT","STP","Slowdown"))
		else:
			print(formatted_string % ("W#","Algorithm","BenchID","Name","Cluster","Mask/NR_WAYS","SlowdownNB/ANTT","STP","Slowdown"))

	for i,benchmark in enumerate(workload):	
		if  mapping_output:
			print(formatted_string % (workload_name,alg_name,i+1,benchmark.name,cluster_id[i], benchmark.llc_id, per_app_masks[i]+"("+str(per_app_ways[i])+")",round(snbw[i],5),round(1.0/sbw[i],5),round(sbw[i],5)))
		else:
			print(formatted_string % (workload_name,alg_name,i+1,benchmark.name,cluster_id[i],per_app_masks[i]+"("+str(round(per_app_ways[i]))+")",round(snbw[i],5),round(1.0/sbw[i],5),round(sbw[i],5)))

	time_string= "%.6fs" % statistics["sim_time"]
	if  mapping_output:
		print(formatted_string % (workload_name,alg_name,len(workload)+1,"OVERALL","#"+str(statistics["total_branches"]),-1,time_string,round(metrics["antt"],5),round(metrics["stp"],5),round(metrics["unfairness"],5)))
		print(formatted_string % (workload_name,alg_name,-1,len(workload)+1,"OVERALL","-","-",round(metrics["gmean_speedup"],5),round(metrics["gmean_ipc"],5),round(metrics["unfairness-cov"],5)))

	else:
		print(formatted_string % (workload_name,alg_name,len(workload)+1,"OVERALL","#"+str(statistics["total_branches"]),time_string,round(metrics["antt"],5),round(metrics["stp"],5),round(metrics["unfairness"],5)))
		print(formatted_string % (workload_name,alg_name,len(workload)+1,"OVERALL","-","-",round(metrics["gmean_speedup"],5),round(metrics["gmean_ipc"],5),round(metrics["unfairness-cov"],5)))


	if "times" in statistics:
		times=statistics["times"]
		for i,t in enumerate(times):
			print("W%i,%s,%d,%d,%.6f" % (idx,sched_title,len(workload),i,t), file=sys.stderr)


def sim_print_sol_dataframe(sched_title,idx,sol_data,max_bandwidth=float('Inf'),print_header=False,**kwargs):

	kwargs.setdefault("user_options",None)
	kwargs.setdefault("workload_prefix","W")
	kwargs.setdefault("alg_suffix","")
	kwargs.setdefault("use_csv",False)

	if "user_options" in kwargs:
		## Merge dict with user options... [Overriding defaults...]
		kwargs.update(kwargs["user_options"])

	workload_prefix=kwargs["workload_prefix"]
	alg_suffix=kwargs["alg_suffix"]
	use_csv=kwargs["use_csv"]

	metrics=compute_basic_metrics(sol_data,max_bandwidth)
	(sol_spec,statistics)=sol_data
	(workload,per_app_ways,per_app_masks,cluster_id)=sol_spec
	fix_intel_bug(per_app_ways,per_app_masks,cluster_id)


	## Determine slowdown with and widthout BW
	sbw=metrics["slowdown"]
	snbw=metrics["slowdown_no_bw"]
	pred_ipcv=metrics["predicted_ipc"]
	alone_ipcv=metrics["ipc_alone"] 

# ("W#","Algorithm","BenchID","Name","Cluster","Mask/NR_WAYS","SlowdownNB/ANTT","STP","Slowdown")

	if use_csv:
		formatted_string="%s,%s,%s,%s,%s,%s"
	else:
		formatted_string="%-Xs %-Ys %-9s %-17s %-17s %-17s"
		formatted_string=formatted_string.replace("X",str(len(workload_prefix)+len("workload")+3)).replace("Y",str(len(alg_suffix)+19))

	workload_name="%s%i" % (workload_prefix,idx)
	alg_name="%s%s" % (sched_title,alg_suffix)

	if print_header:
		print(formatted_string % ("workload","scheme","id","app","property","value"))

	for i,benchmark in enumerate(workload):	
		print(formatted_string % (workload_name,alg_name,i+1,benchmark.name,"cluster_id",cluster_id[i]))
		print(formatted_string % (workload_name,alg_name,i+1,benchmark.name,"llc_id",benchmark.llc_id))
		print(formatted_string % (workload_name,alg_name,i+1,benchmark.name,"way_count",per_app_ways[i]))
		print(formatted_string % (workload_name,alg_name,i+1,benchmark.name,"way_mask",per_app_masks[i]))
		print(formatted_string % (workload_name,alg_name,i+1,benchmark.name,"slowdown",round(sbw[i],5)))
		print(formatted_string % (workload_name,alg_name,i+1,benchmark.name,"slowdown_nb",round(snbw[i],5)))
		print(formatted_string % (workload_name,alg_name,i+1,benchmark.name,"ipc",round(pred_ipcv[i],5)))
		print(formatted_string % (workload_name,alg_name,i+1,benchmark.name,"ipc_alone",round(alone_ipcv[i],5)))  

	## Overall workload output
	time_string= "%.6f" % statistics["sim_time"]

	print(formatted_string % (workload_name,alg_name,0,"OVERALL","STP",round(metrics["stp"],5)))
	print(formatted_string % (workload_name,alg_name,0,"OVERALL","ANTT",round(metrics["antt"],5)))
	print(formatted_string % (workload_name,alg_name,0,"OVERALL","Unfairness",round(metrics["unfairness"],5)))
	print(formatted_string % (workload_name,alg_name,0,"OVERALL","UnfairnessCoV",round(metrics["unfairness-cov"],5)))
	print(formatted_string % (workload_name,alg_name,0,"OVERALL","JainFairness",round(metrics["jain-fairness"],5)))
	print(formatted_string % (workload_name,alg_name,0,"OVERALL","GmeanSpeedup",round(metrics["gmean_speedup"],5)))
	print(formatted_string % (workload_name,alg_name,0,"OVERALL","GmeanIPC",round(metrics["gmean_ipc"],5)))
	print(formatted_string % (workload_name,alg_name,0,"OVERALL","AggregateIPC",round(metrics["aggregate_ipc"],5)))
	print(formatted_string % (workload_name,alg_name,0,"OVERALL","HmeanSpeedup",round(metrics["hmean_speedup"],5)))
	print(formatted_string % (workload_name,alg_name,0,"OVERALL","M1",round(metrics["m1_metric"],5)))
	print(formatted_string % (workload_name,alg_name,0,"OVERALL","sim_time_s",time_string))
	print(formatted_string % (workload_name,alg_name,0,"OVERALL","explored_sols",str(statistics["total_branches"])))

	if "times" in statistics:
		times=statistics["times"]
		for i,t in enumerate(times):
			print(formatted_string % (workload_name,alg_name,0,"OVERALL","sim_time%d" % (i,t)))



def sim_print_sol_masks(sched_title,idx,sol_data,max_bandwidth=float('Inf'),print_header=False):
	(sol_spec,statistics)=sol_data
	(workload,per_app_ways,per_app_masks,cluster_id)=sol_spec
	fix_intel_bug(per_app_ways,per_app_masks,cluster_id)
	for mask in per_app_masks:
		print(mask,)
	print

def sim_print_sol_masks_debussy(sched_title,idx,sol_data,max_bandwidth=float('Inf'),print_header=False):
	banned_benchmarks = ["omnetpp06","povray06","leslie3d06","hmmer06","gamess06","perlbench06","vortex00","gromacs06","gobmk06","h264ref06","sjeng06"]
	(sol_spec,statistics)=sol_data
	(workload,per_app_ways,per_app_masks,cluster_id)=sol_spec

	if len(cluster_id) == len(per_app_masks):
		# Transform cluster_id to app_id vector
		cluster_id = [i for i in range(len(cluster_id))]
	
	cl_ixs_problematics = []
	cl_i_max = -1
	ways_max = -1
	cluster_ways = [ 0 for _ in range(len(cluster_id))]
	cluster_apps_ixs = [ [] for _ in range(len(cluster_id))]

	for i,app in enumerate(workload):
		cluster_apps_ixs[cluster_id[i]].append(i)
		cluster_ways[cluster_id[i]] = per_app_ways[i]
		if per_app_ways[i] > ways_max:
			ways_max = per_app_ways[i] 
			cl_i_max = cluster_id[i]
		if app.name in banned_benchmarks and per_app_ways[i] == 1:
			if len(cl_ixs_problematics) < 2:
				if cluster_id[i] not in cl_ixs_problematics:
					cl_ixs_problematics.append(cluster_id[i]) 

	if cl_ixs_problematics:
		# Assign first problematic clust
		ways_offset = 1 
		cl_id = cl_ixs_problematics[0]
		for app_id in cluster_apps_ixs[cl_id]:
			per_app_masks[app_id] = "0x3"

		# Assign max clust
		for i_max in cluster_apps_ixs[cl_i_max]:
			app_ways = per_app_ways[i_max]
			per_app_masks[i_max] = hex((1<<app_ways)-1 <<ways_offset)

		ways_offset+=app_ways 

		# Assign second problematic clust
		if len(cl_ixs_problematics) == 2:
			cl_id = cl_ixs_problematics[1]
			for app_id in cluster_apps_ixs[cl_id]:
				per_app_masks[app_id] = hex((1<<2)-1 << ways_offset-1)
			ways_offset += 1

		# Assign the rest of applications taking into account the offset
		for cl_idx in range(len(cluster_id)):
			if cl_idx not in cl_ixs_problematics and cl_idx != cl_i_max:
				cl_ways = cluster_ways[cl_idx]
				for app_id in cluster_apps_ixs[cl_idx]:
					per_app_masks[app_id] = hex((1<<cl_ways)-1 << ways_offset)
				ways_offset+=cl_ways


	for mask in per_app_masks:
		print(mask,)
	print





def sim_print_cluster_info(sched_title,idx,sol_data,max_bandwidth=float('Inf'),print_masks=False):
	(sol_spec,statistics)=sol_data
	(workload,per_app_ways,per_app_masks,cluster_id)=sol_spec

	## Issue if mask end bit


	#print workload_data
	ways={}
	masks={}
	llc_id={}
	cluster_llc=[]

	cluster_zero=False
	llc_id_matters=False

	for idx,benchmark in enumerate(workload):	
		cid=cluster_id[idx]
		ways[cid]=per_app_ways[idx]
		masks[cid]=per_app_masks[idx]
		llc_id[cid]=benchmark.llc_id	
		cluster_llc.append(benchmark.llc_id)

		if cid==0:
			cluster_zero=True

		if not llc_id_matters and benchmark.llc_id!=-1:
			llc_id_matters=True


	if print_masks:
		if sched_title!="dunn":
			masks=fix_intel_bug_aux(ways,masks,masks.keys())
		target=masks
	else:
		target=ways


	## Deal with the issue of gaps during cluster assignment
	if sched_title=="dunn":
		sorted_cluster_ids=sorted(target.keys())
		cluster_spec=','.join([str(target[i]) for i in sorted_cluster_ids])
		cluster_map=[sorted_cluster_ids.index(i) for i in cluster_id]
		print("%s;%s" % (','.join(map(str, cluster_map)),cluster_spec))
	else:
		start=0 if cluster_zero else 1

		if not llc_id_matters:
			## Print out things
			cluster_spec=','.join([str(target[i]) for i in range(start,len(target))])
			print("%s;%s" % (','.join(map(str, cluster_id)),cluster_spec))
		else:
			## Combine pairs idx_cluster-app/llc_group;
			cluster_spec=[ "%s/%d" % (target[i],llc_id[i])  for i in range (start,len(target)) ]		
			cluster_spec_str=','.join(cluster_spec)
			id_clusters=','.join(map(str, cluster_id))
			id_groups=','.join(map(str, cluster_llc))
			print("%s;%s;%s" % (id_clusters,id_groups,cluster_spec_str))			



def initialize_simulator(args):
	select_bw_model(args.bw_model)
	select_topology(args.topology)
	## Must be commented out in public version
	add_experimental_algorithms()

## Parameters
# ------------
# algoritm-name (string), 
# workload: list of apps in the workload 
# nr_ways: total number of ways on the machine
# max_bandwidth: maximum bandwidth of the platform...
## Return value
# ------------ 
# Tuple with two mandatory items
# 1st item [solution]: (nr_ways_per_app, cat_masks, cluster_id_app)
# 2nd item [statistics dict]: (time to reach solution, total branches explored)
def apply_part_algorithm(algorithm,workload,nr_ways,max_bandwidth=float('Inf'),parallel=False,debugging=False,uoptions={},workload_name="W",opt_spec=(unfairness,"False")):
	solution=None
	total_branches=1
	stats={}
	metadata = {}
	iterations=1
	time_array=[]
	time_model=False

	if "iterations" in uoptions:
		iterations=uoptions["iterations"]
		if iterations<=1:
			iterations=1

	if "time_model" in uoptions:
		time_model=uoptions["time_model"]

	for i in range(iterations):
		## Keep track of sim time
		start = time.time()
		if algorithm=="opt-stp":
			(cost,solution,total_branches,metadata)=get_optimal_schedule(workload,throughput,True,nr_ways,max_bandwidth,multiprocessing=parallel,user_options=uoptions)
		elif algorithm=="opt-unf":
			(cost,solution,total_branches,metadata)=get_optimal_schedule(workload,unfairness_max_throughput,False,nr_ways,max_bandwidth,multiprocessing=parallel,user_options=uoptions)
		elif algorithm=="opt-stp-s":
			(cost,solution,total_branches,metadata)=get_optimal_schedule(workload,throughput,True,nr_ways,max_bandwidth,multiprocessing=False,user_options=uoptions)
		elif algorithm=="opt-unf-s":
			(cost,solution,total_branches,metadata)=get_optimal_schedule(workload,unfairness_max_throughput,False,nr_ways,max_bandwidth,multiprocessing=False,user_options=uoptions)
		elif algorithm=="opt-stp-bf":
			(cost,solution,total_branches,metadata)=get_optimal_schedule_bf(workload,throughput,True,nr_ways,max_bandwidth,multiprocessing=parallel,user_options=uoptions)
		elif algorithm=="opt-unf-bf":
			(cost,solution,total_branches,metadata)=get_optimal_schedule_bf(workload,unfairness_max_throughput,False,nr_ways,max_bandwidth,multiprocessing=parallel,user_options=uoptions)
		elif algorithm=="opt-bf":
			(cost,solution,total_branches,metadata)=get_optimal_schedule_bf(workload,opt_spec[0],opt_spec[1],nr_ways,max_bandwidth,multiprocessing=parallel,user_options=uoptions)
		elif algorithm=="optc-stp":
			(patched_workload,per_app_ways,per_app_masks,cluster_id,total_branches,metadata)=get_optimal_clustering(workload,throughput,True,nr_ways,max_bandwidth,multiprocessing,user_options=uoptions)
		elif algorithm=="optc-unf":
			(patched_workload,per_app_ways,per_app_masks,cluster_id,total_branches,metadata)=get_optimal_clustering(workload,unfairness_max_throughput,False,nr_ways,max_bandwidth,multiprocessing=parallel,user_options=uoptions)
		elif algorithm=="optc-stp-bf":
			(patched_workload,per_app_ways,per_app_masks,cluster_id,total_branches,metadata)=get_optimal_clustering_bf(workload,throughput,True,nr_ways,max_bandwidth,multiprocessing,user_options=uoptions)
		elif algorithm=="optc-unf-bf":
			(patched_workload,per_app_ways,per_app_masks,cluster_id,total_branches,metadata)=get_optimal_clustering_bf(workload,unfairness_max_throughput,False,nr_ways,max_bandwidth,multiprocessing=parallel,user_options=uoptions)
		elif algorithm=="optc":
			(patched_workload,per_app_ways,per_app_masks,cluster_id,total_branches,metadata)=get_optimal_clustering_bf(workload,opt_spec[0],opt_spec[1],nr_ways,max_bandwidth,multiprocessing=parallel,user_options=uoptions)
		elif algorithm=="opt-map-stp":
			uopt_patched=dict(uoptions)
			uopt_patched["opt_mapping"]=True
			(patched_workload,per_app_ways,per_app_masks,cluster_id,total_branches,metadata)=get_optimal_clustering_bf(workload,throughput,True,nr_ways,max_bandwidth,multiprocessing,user_options=uopt_patched)
		elif algorithm=="opt-map-unf":
			uopt_patched=dict(uoptions)
			uopt_patched["opt_mapping"]=True
			(patched_workload,per_app_ways,per_app_masks,cluster_id,total_branches,metadata)=get_optimal_clustering_bf(workload,unfairness_max_throughput,False,nr_ways,max_bandwidth,multiprocessing=parallel,user_options=uopt_patched)
		elif algorithm=="opt-map-cov":
			uopt_patched=dict(uoptions)
			uopt_patched["opt_mapping"]=True
			(patched_workload,per_app_ways,per_app_masks,cluster_id,total_branches,metadata)=get_optimal_clustering_bf(workload,cov_unfairness_metric,False,nr_ways,max_bandwidth,multiprocessing=parallel,user_options=uopt_patched)
		elif algorithm=="opt-map-fair":
			uopt_patched=dict(uoptions)
			uopt_patched["opt_mapping"]=True
			(patched_workload,per_app_ways,per_app_masks,cluster_id,total_branches,metadata)=get_optimal_clustering_bf(workload,fairness_metric,True,nr_ways,max_bandwidth,multiprocessing=parallel,user_options=uopt_patched)
		elif algorithm=="opt-map-antt":
			uopt_patched=dict(uoptions)
			uopt_patched["opt_mapping"]=True
			(patched_workload,per_app_ways,per_app_masks,cluster_id,total_branches,metadata)=get_optimal_clustering_bf(workload,antt_metric,False,nr_ways,max_bandwidth,multiprocessing=parallel,user_options=uopt_patched)			
		elif algorithm=="opt-map":
			uopt_patched=dict(uoptions)
			uopt_patched["opt_mapping"]=True
			(patched_workload,per_app_ways,per_app_masks,cluster_id,total_branches,metadata)=get_optimal_clustering_bf(workload,opt_spec[0],opt_spec[1],nr_ways,max_bandwidth,multiprocessing=parallel,user_options=uopt_patched)	
		elif algorithm=="yu-petrov":
			solution=get_schedule_yu_petrov(workload,nr_ways)
		elif algorithm=="ucp":
			if "ucp_metric" in uoptions:
				solution=get_schedule_UCP_gen(workload,nr_ways,metric=uoptions["ucp_metric"])
			else:
				solution=get_schedule_UCP(workload,nr_ways)
		elif algorithm=="ucp-slowdown":
			solution=get_schedule_UCP_gen(workload,nr_ways,metric="slowdown")
		elif algorithm=="kpart":
			(patched_workload,per_app_ways,per_app_masks,cluster_id)=get_kpart_best_gen(workload,nr_ways,max_bandwidth,debugging=debugging,variant=0,user_options=uoptions)
		elif algorithm=="kpart-opt":
			(patched_workload,per_app_ways,per_app_masks,cluster_id)=get_kpart_best_gen(workload,nr_ways,max_bandwidth,debugging=debugging,variant=1,user_options=uoptions)
		elif algorithm=="whirlpool":
			(patched_workload,per_app_ways,per_app_masks,cluster_id)=get_schedule_whirlpool_float(workload, nr_ways)
		elif algorithm=="whirlpool-c":
			(patched_workload,per_app_ways,per_app_masks,cluster_id)=get_schedule_whirlpool_float(workload, nr_ways,metric="llcmpkc")
		elif algorithm == "user":
			(patched_workload,per_app_ways,per_app_masks,cluster_id) = get_user_assignment(workload,nr_ways,workload_name,user_options=uoptions)
		elif algorithm=="equal-part":
			solution=get_equal_llc_schedule(len(workload),nr_ways)
		else:
			## Must be commented out in public version
			(patched_workload,per_app_ways,per_app_masks,cluster_id)=invoke_extra_algorithm(algorithm,workload,nr_ways,max_bandwidth,parallel,debugging,uoptions,opt_spec)
			#print "algorithm not valid: %s" % algorithm
			#exit(1)		

		if time_model:
			unfairness_max_throughput(workload, solution, max_bandwidth)

		end = time.time()
		sim_time=round(end - start,6)
		time_array.append(sim_time)

		## Hack to enable successful clean up of communications in remote engines
		if "broadcast" in uoptions and  uoptions["broadcast"]:
			time.sleep(0.5)

	if "trace" in metadata:
		trace_fd = open("paraver_trace_%s_%s.prv" % (workload_name,algorithm),"w+")
		trace_fd.write(metadata["header"])
		for line in metadata["trace"]:
			trace_fd.write(":".join(line)+"\n")
		trace_fd.close()


	## Build solution based on return values
	if not solution:
		solution_spec=(patched_workload,per_app_ways,per_app_masks,cluster_id)
	else:
		solution_spec=(workload,solution,get_partition_masks(solution),[i for i in range(1,len(workload)+1)])


	if iterations>1:
		stats["times"]=time_array
		stats["sim_time"]=np.array(time_array).mean()
	else:
		stats["sim_time"]=sim_time
	stats["total_branches"]=total_branches
	return (solution_spec,stats)



## function to set up global data
gbl_sim_dict={}
def set_global_simulation_properties(key_vals):
	global gbl_sim_dict
	gbl_sim_dict.update(key_vals)
	## Establish BW model
	if "bw_model" in key_vals:
		select_bw_model(key_vals["bw_model"])

def get_global_simulation_properties(foo=None):
	global gbl_sim_dict
	return gbl_sim_dict 

def execute_simulation(simulation):
	## simulation is just a pair (idx_workload, idx_algorithm )
	(idx_workload, idx_algorithm)=simulation
	params=get_global_simulation_properties()
	wname="W%d" % (idx_workload+1)
	workloads=params["workloads"]
	uoptions=params["uoptions"]
	algorithms=params["algorithms"]
	nr_ways=params["nr_ways"]
	max_bandwidth=params["max_bandwidth"]

	alg=algorithms[idx_algorithm]
	algorithm=alg[0]
	opt_spec=alg[1]

	return apply_part_algorithm(algorithm,workloads[idx_workload],nr_ways,max_bandwidth,parallel=False,debugging=False,uoptions=uoptions,workload_name=wname,opt_spec=opt_spec)



def launch_simulations_in_parallel(workloads,algorithms,workload_range,nr_ways,max_bandwidth, uoptions,bw_model):
	parallel=False ## Forced
	debugging=False #True

	# Establish connection
	rc = ipp.Client(timeout=30000) 
	nr_engines=len(rc.ids)
	# Set up global simulation parameters in engines
	dview = rc[:]
	dict_globals={"workloads":workloads,
			"max_bandwidth": max_bandwidth,
			"algorithms": algorithms,
			"nr_ways": nr_ways,
			"uoptions": uoptions,
			"bw_model": bw_model}

	if debugging:
		set_global_simulation_properties(dict_globals)
	
	# Copy global data
	ret=dview.apply_sync(lambda x: set_global_simulation_properties(x),dict_globals)

	## Generate simulation set 
	simulations=[]
	for widx in workload_range:
		for alg_idx,algorithm in enumerate(algorithms):
			simulations.append((widx,alg_idx))

	## Execute simulations in parallel
	lview = rc.load_balanced_view(block=True)

	if debugging:
		results=[r for r in map(lambda simdata: execute_simulation(simdata), simulations)]
	else:
		results=lview.map(lambda simdata: execute_simulation(simdata), simulations)

	## Turn results into a dictionary
	resdict={}
	for i,simdata in enumerate(simulations):
		resdict[simdata]=results[i]

	## Return dict with simulation results
	return resdict
