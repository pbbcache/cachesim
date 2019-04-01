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


#Function to pretty-print a workload 
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

	if "per_app_centroids" in statistics:
		cluster_id = statistics["per_app_centroids"]

	## Determine slowdown with and widthout BW
	sbw=metrics["slowdown"]
	snbw=metrics["slowdown_no_bw"]

	formatted_string="%-4s %-19s %-7s %-20s %-13s %-15s %-17s %-17s %-17s"
	if print_header:
		print formatted_string % ("W#","Algorithm","BenchID","Name","Cluster","Mask/NR_WAYS","SlowdownNB/ANTT","STP","Slowdown")

	for i,benchmark in enumerate(workload):
		print formatted_string % ("W%i" % idx,sched_title,i+1,benchmark.name,cluster_id[i],per_app_masks[i]+"("+str(per_app_ways[i])+")",round(snbw[i],5),round(1.0/sbw[i],5),round(sbw[i],5))

	print formatted_string % ("W%i" % idx,sched_title,len(workload)+1,"OVERALL","#"+str(statistics["total_branches"]),statistics["sim_time"]+"s",round(metrics["antt"],5),round(metrics["stp"],5),round(metrics["unfairness"],5))



def sim_print_sol_masks(sched_title,idx,sol_data,max_bandwidth=float('Inf'),print_header=False):
	(sol_spec,statistics)=sol_data
	(workload,per_app_ways,per_app_masks,cluster_id)=sol_spec
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
		#Assign first problematic clust
		ways_offset = 1 
		cl_id = cl_ixs_problematics[0]
		for app_id in cluster_apps_ixs[cl_id]:
			per_app_masks[app_id] = "0x3"

		#Assign max clust
		for i_max in cluster_apps_ixs[cl_i_max]:
			app_ways = per_app_ways[i_max]
			per_app_masks[i_max] = hex((1<<app_ways)-1 <<ways_offset)

		ways_offset+=app_ways 

		#Assign second problematic clust
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


glb_algorithms_supported=["opt-stp","yu-petrof","opt-unf","ucp","ucp-slowdown","kpart","dunn","whirlpool","whirlpool-c","equal-part","optc-stp","optc-unf","opt-unf-np","opt-stp-np","yu-petrof"]

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
	elif algorithm=="optc-stp":
		(patched_workload,per_app_ways,per_app_masks,cluster_id,total_branches,metadata)=get_optimal_clustering(workload,throughput,True,nr_ways,max_bandwidth,multiprocessing,user_options=uoptions)
	elif algorithm=="optc-unf":
		(patched_workload,per_app_ways,per_app_masks,cluster_id,total_branches,metadata)=get_optimal_clustering(workload,unfairness_max_throughput,False,nr_ways,max_bandwidth,multiprocessing=parallel,user_options=uoptions)
	elif algorithm=="yu-petrof":
		solution=get_schedule_yu_petrof(workload,nr_ways)
	elif algorithm=="ucp":
		solution=get_schedule_UCP(workload,nr_ways)
	elif algorithm=="ucp-slowdown":
		solution=get_schedule_UCP_gen(workload,nr_ways,metric="slowdown")
	elif algorithm=="kpart":
		(patched_workload,per_app_ways,per_app_masks,cluster_id)=get_kpart_best_gen(workload,nr_ways,max_bandwidth,debugging)
	elif algorithm=="dunn":
		patched_workload=workload
		(per_app_ways,per_app_masks,per_app_centroids,cluster_id)=get_schedule_dunn(workload, nr_ways)
		stats["per_app_centroids"]=per_app_centroids
	elif algorithm=="whirlpool":
		(patched_workload,per_app_ways,per_app_masks,cluster_id)=get_schedule_whirlpool_float(workload, nr_ways)
	elif algorithm=="whirlpool-c":
		(patched_workload,per_app_ways,per_app_masks,cluster_id)=get_schedule_whirlpool_float(workload, nr_ways,metric="llcmpkc")
	elif algorithm == "stock-linux":
		solution = [ hex(((1<<nr_ways)-1)<<0) for _ in  range(len(workload))]
	elif algorithm=="equal-part":
		solution=get_equal_llc_schedule(len(workload),nr_ways)
	else:
		print "algorithm not valid: %s" % algorithm
		exit(1)
	end = time.time()
	sim_time=str(round(end - start,4))

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

	stats["sim_time"]=sim_time
	stats["total_branches"]=total_branches
	return (solution_spec,stats)
