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
from simulator_core import *
from opt_clustering import *
import numpy as np

## Efficient function to determine all metrics from slowdown vector...
## Returns a dictionary with the collected metrics....
def compute_basic_metrics(sol_data,max_bandwidth=float('Inf')):
	metrics={}
	(sol_spec,statistics)=sol_data
	(workload,per_app_ways,per_app_masks,cluster_id)=sol_spec
	slowdown_vector=get_slowdown_vector(workload,per_app_ways,max_bandwidth)
	slowdown_vector_nobw=get_slowdown_vector(workload,per_app_ways)

	metrics["unfairness"] = max(slowdown_vector) / min(slowdown_vector)
	metrics["antt"] =sum(slowdown_vector) / len(slowdown_vector)
	stp=0
	for slowdown in slowdown_vector:
		stp += 1/slowdown

	metrics["stp"] = stp 
	metrics["slowdown"]= slowdown_vector
	metrics["slowdown_no_bw"]= slowdown_vector_nobw

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
	print "****WORKLOAD #%d****" % id
	bnames=[]
	for benchmark in workload:
		bnames.append(benchmark.name)
	print bnames

def sim_print_sol_simple(sched_title,idx,sol_data,max_bandwidth=float('Inf'),print_header=False):
	metrics=compute_basic_metrics(sol_data,max_bandwidth)
	(sol_spec,statistics)=sol_data
	(workload,per_app_ways,per_app_masks,cluster_id)=sol_spec

	## Determine slowdown with and widthout BW
	sbw=metrics["slowdown"]
	snbw=metrics["slowdown_no_bw"]

	if print_header:
		print_workload(workload,idx)

	print "==== Solution for %s ===" % sched_title
	i=0
	for benchmark in workload:
		print "%s --> (%d ways,%f,%f)" % (benchmark.name,per_app_ways[i],sbw[i],snbw[i])
		i=i+1
	print "--------------"
	print "--------------"
	print "STP=%f" % metrics["stp"]
	print "Unf=%f" % metrics["unfairness"]
	print "--------------"
	print "======================="

def sim_print_sol_table(sched_title,idx,sol_data,max_bandwidth=float('Inf'),print_header=False):

	metrics=compute_basic_metrics(sol_data,max_bandwidth)
	(sol_spec,statistics)=sol_data
	(workload,per_app_ways,per_app_masks,cluster_id)=sol_spec
	fix_intel_bug(per_app_ways,per_app_masks,cluster_id)


	## Determine slowdown with and widthout BW
	sbw=metrics["slowdown"]
	snbw=metrics["slowdown_no_bw"]

	formatted_string="%-4s %-19s %-7s %-20s %-13s %-15s %-17s %-17s %-17s"
	if print_header:
		print formatted_string % ("W#","Algorithm","BenchID","Name","Cluster","Mask/NR_WAYS","SlowdownNB/ANTT","STP","Slowdown")

	for i,benchmark in enumerate(workload):
		print formatted_string % ("W%i" % idx,sched_title,i+1,benchmark.name,cluster_id[i],per_app_masks[i]+"("+str(per_app_ways[i])+")",round(snbw[i],5),round(1.0/sbw[i],5),round(sbw[i],5))

	time_string= "%.6fs" % statistics["sim_time"]
	print formatted_string % ("W%i" % idx,sched_title,len(workload)+1,"OVERALL","#"+str(statistics["total_branches"]),time_string,round(metrics["antt"],5),round(metrics["stp"],5),round(metrics["unfairness"],5))

	if "times" in statistics:
		times=statistics["times"]
		for i,t in enumerate(times):
			print >> sys.stderr, "W%i,%s,%d,%d,%.6f" % (idx,sched_title,len(workload),i,t)

def sim_print_sol_masks(sched_title,idx,sol_data,max_bandwidth=float('Inf'),print_header=False):
	(sol_spec,statistics)=sol_data
	(workload,per_app_ways,per_app_masks,cluster_id)=sol_spec
	fix_intel_bug(per_app_ways,per_app_masks,cluster_id)
	for mask in per_app_masks:
		print mask,
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
		print mask,
	print





def sim_print_cluster_info(sched_title,idx,sol_data,max_bandwidth=float('Inf'),print_masks=False):
	(sol_spec,statistics)=sol_data
	(workload,per_app_ways,per_app_masks,cluster_id)=sol_spec

	## Issue if mask end bit


	#print workload_data
	ways={}
	masks={}

	cluster_zero=False
	for idx in range(len(workload)):
		cid=cluster_id[idx]
		ways[cid]=per_app_ways[idx]
		masks[cid]=per_app_masks[idx]

		if cid==0:
			cluster_zero=True


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
		print "%s;%s" % (','.join(map(str, cluster_map)),cluster_spec) 
	else:
		## Print out things
		if cluster_zero:
			cluster_spec=','.join([str(target[i]) for i in range(0,len(target))])
		else:
			cluster_spec=','.join([str(target[i]) for i in range(1,len(target)+1)])
	
		print "%s;%s" % (','.join(map(str, cluster_id)),cluster_spec) 


def initialize_simulator(bw_model):
	select_bw_model(bw_model)

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
def apply_part_algorithm(algorithm,workload,nr_ways,max_bandwidth,parallel=False,debugging=False,uoptions={},workload_name="W"):
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
		elif algorithm=="optc-stp":
			(patched_workload,per_app_ways,per_app_masks,cluster_id,total_branches,metadata)=get_optimal_clustering(workload,throughput,True,nr_ways,max_bandwidth,multiprocessing,user_options=uoptions)
		elif algorithm=="optc-unf":
			(patched_workload,per_app_ways,per_app_masks,cluster_id,total_branches,metadata)=get_optimal_clustering(workload,unfairness_max_throughput,False,nr_ways,max_bandwidth,multiprocessing=parallel,user_options=uoptions)
		elif algorithm=="optc-stp-bf":
			(patched_workload,per_app_ways,per_app_masks,cluster_id,total_branches,metadata)=get_optimal_clustering_bf(workload,throughput,True,nr_ways,max_bandwidth,multiprocessing,user_options=uoptions)
		elif algorithm=="optc-unf-bf":
			(patched_workload,per_app_ways,per_app_masks,cluster_id,total_branches,metadata)=get_optimal_clustering_bf(workload,unfairness_max_throughput,False,nr_ways,max_bandwidth,multiprocessing=parallel,user_options=uoptions)
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
			print "algorithm not valid: %s" % algorithm
			exit(1)		

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

	return apply_part_algorithm(algorithms[idx_algorithm],workloads[idx_workload],nr_ways,max_bandwidth,parallel=False,debugging=False,uoptions=uoptions,workload_name=wname)



def launch_simulations_in_parallel(workloads,algorithms,workload_range,nr_ways,max_bandwidth, uoptions,bw_model):
	parallel=False ## Forced
	debugging=False

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
		results=map(lambda simdata: execute_simulation(simdata), simulations)
	else:
		results=lview.map(lambda simdata: execute_simulation(simdata), simulations)

	## Turn results into a dictionary
	resdict={}
	for i,simdata in enumerate(simulations):
		resdict[simdata]=results[i]

	## Return dict with simulation results
	return resdict
