#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
# simulator_core.py
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
import operator
import csv
import sys
import pandas as pd
import numpy as np
import sympy as sp
import multiprocessing
from multiprocessing.pool import ThreadPool
import ipyparallel as ipp
from subprocess import call
from sklearn.cluster import KMeans
import time
import functools
import datetime
import pytz

# Best 3, 6 (split / and get last)
def parse_harness_file(filename):
	inputfile=open(filename, "r")
	benchmarks=[] ## Empty dict
	## Skip header	
	line=inputfile.readline()
	i=0
	while line!="":
		if i%3==2:
			row=line.split("/")
			benchmarks.append(row[-1].strip())
		i=i+1
		line=inputfile.readline()
	return benchmarks
	

################################################################################
##             BANDWITH SLOWDOWN MODEL [Morad et al, JPDC'16]                 ##
################################################################################

## 
# Function to predict the bandwidth alone
# Params: Total bandwidth used by the workload and application actual bandwidth.
# Return: Predicted bandwidth alone.
def estimate_bw_alone(total_bw, bw_shared):
	return (bw_shared*(bw_shared+1-total_bw))/(bw_shared*(bw_shared+1-total_bw)+1-total_bw)

## 
# Approximates the slowdown due to collisions on the memory bus.  
# Params: Total bandwidth used by the workload and application actual bandwidth.
# Return: Predicted slowdown
def predict_slowdown(total_bw,bw_shared):
	bw_alone=estimate_bw_alone(total_bw,bw_shared)
	return bw_alone/bw_shared


## 
# Determines the theoretical total memory bandwidth to guarantee that
# an application suffers a slowdown matching the target passed as a parameter  
# Params: Target slowdown we want to reach, actual app BW (when running with the others)
#         current Total bandwidth 
# Return: Total estimated BW to guarantee target.
def get_total_bandwidth(target_slowdown,bw_shared,cur_total_bw):
	bw_alone=estimate_bw_alone(cur_total_bw,bw_shared)
	new_bw_shared=bw_alone/target_slowdown
	new_total=1-((new_bw_shared**2)*(1/bw_alone -1))/(new_bw_shared*(1-1/bw_alone) + 1)

	return new_total

## 
# Given the bandwidth (BW) solo utilization with respect to the theoretical maximum 
# for n applications, this function predicts the actual aggregate BW 
# and per-application BW for each application when running together.  
# Params: BW alone vector (consisting of normalized values [0..1])
# Return: Collection with the values of the variables
#	[0]-> Total BW
#	[1]-> BW for app 0
#	[2]-> BW for app 1
# 	....
def estimate_bw_shared(bw_alone_vector):
	n=len(bw_alone_vector)
	equations=[]
	U=sp.Symbol('U')
	variables=[U]
	guess=[sum(bw_alone_vector)]
	## To build equations dynamically
	eq=-U

	for i in range(n):
		bi=bw_alone_vector[i]
		ui=sp.Symbol('u%d' % i)
		variables.append(ui)
		coeff=(1.0-(1.0/bi))
		equations.append((ui**2)*coeff+ui*(1-U)*coeff+1-U)
		eq+=ui
		guess.append(bi)

	equations.append(eq)#"-U+u0+u1+u2")#)eq)

	try:
		return sp.nsolve(equations,variables,guess) #verify=False)
	except:
		return guess



def apply_bw_model(bw_alone_vector,slowdown_vector,max_bandwidth):
	bw_alone_vnorm=map(lambda x: x/max_bandwidth,bw_alone_vector)

	## Invoke solver to determine approximate solution
	solver_solution=estimate_bw_shared(bw_alone_vnorm)
	
	## Get max predicted BW from solution
	predicted_total_bw=solver_solution[0]

	# Get BW shared for each app
	predicted_bw_shared=solver_solution[1:]

	overall_slowdown=[]
	
	## Determine actual slowdown by considering BW-related slowdown 
	for i,bw_alone in enumerate(bw_alone_vnorm): 
		## Simple model (Slowdown is proportional to BW reduction)
		bw_slowdown_factor=bw_alone/predicted_bw_shared[i]
		## Calculate performance reduction 
		overall_slowdown.append(slowdown_vector[i]*bw_slowdown_factor)

	return overall_slowdown

## Global simulator parameter
glb_active_bw_model=0
glb_available_bw_models=["simple","stalls"]

def select_bw_model(str_model):
	global glb_active_bw_model
	if str_model in glb_available_bw_models:
		glb_active_bw_model=glb_available_bw_models.index(str_model)
	else:
		print >> sys.stderr, "No such BW model: ", str_model
		exit(1)

def apply_bw_model_gen(benchmark_list,ways_assigned,max_bandwidth):
	slowdown_vector = []
	bw_alone_vector = []

	global glb_active_bw_model
	#print glb_active_bw_model
	#exit(0)
	if glb_active_bw_model==1: ## Stall_based
		total_stalls_alone_v=[]
		mem_access_stalls_alone_v=[]

		if "stalls_l3_miss" in benchmark_list[0].properties.columns:
			stall_event="stalls_l3_miss"
		else:
			stall_event="stalls_ldm_pending"

	## Retrieve data from properties as required by the model...
	for i,benchmark in enumerate(benchmark_list):
		slowdown_vector.append(benchmark.properties["slowdown"][ways_assigned[i]])
		bw_alone_vector.append(benchmark.properties["bandwidth_mbps"][ways_assigned[i]])
		if glb_active_bw_model==1:
			s_total=benchmark.properties["stalls_total"][ways_assigned[i]] 
			total_stalls_alone_v.append(s_total)
			## Probar a meter los store stalls
			mem_stalls=benchmark.properties[stall_event][ways_assigned[i]]+benchmark.properties["stalls_store_buf"][ways_assigned[i]]
			if mem_stalls>s_total:
				mem_stalls=s_total ## Truncate
			mem_access_stalls_alone_v.append(mem_stalls)

	bw_alone_vnorm=map(lambda x: x/max_bandwidth,bw_alone_vector)

	## Invoke solver to determine approximate solution
	solver_solution=estimate_bw_shared(bw_alone_vnorm)
	
	## Get max predicted BW from solution
	predicted_total_bw=solver_solution[0]

	# Get BW shared for each app
	predicted_bw_shared=solver_solution[1:]

	overall_slowdown=[]
	
	## Determine actual slowdown by considering BW-related slowdown 
	for i,bw_alone in enumerate(bw_alone_vnorm): 
		## Simple model (Slowdown is proportional to BW reduction)
		bw_slowdown_factor=bw_alone/predicted_bw_shared[i]

		## Alter bw_slowdown_factor according to stall-based model
		if glb_active_bw_model==1:
			s_alone=total_stalls_alone_v[i]
			m_alone=mem_access_stalls_alone_v[i]
			## bw_slowdown_factor is Old value provided by Morad's model
			## Determine new bw_slowdown factor
			old_morad=bw_slowdown_factor
			bw_slowdown_factor=(s_alone+m_alone*(old_morad-1))/s_alone
			#assert(old_morad<bw_slowdown_factor)

			#bw_slowdown_factor=s_shared/s_alone

		## Calculate performance reduction 
		overall_slowdown.append(slowdown_vector[i]*bw_slowdown_factor)

	return overall_slowdown


def determine_bw_shared(bw_alone_vector,max_bandwidth):
	bw_alone_vnorm=map(lambda x: x/max_bandwidth,bw_alone_vector)

	## Invoke solver to determine approximate solution
	solver_solution=estimate_bw_shared(bw_alone_vnorm)

	# Return Agg BW and BW shared for each app
	return (solver_solution[0]*max_bandwidth,map(lambda x: x* max_bandwidth,solver_solution[1:]))		




################################################################################
##                                 COST FUNCTIONS                             ##
################################################################################

## Sequential application with information on performance, 
## LLCMPKI and memory bandwidth for different cache sizes  
class App:
	def __init__(self, name, properties, parent=None):
		self.name = name
		self.properties=properties
		self.scaling_factor=1.
		## For kpart
		if parent:
			self.original_app=parent.original_app
		else:
			self.original_app=self		


	def __repr__(self):
		return repr((self.name))

	def show(self):
		print "==== %s ====" % self.name 
		print self.properties

	def get_metric_table(self,metric):
		return self.properties.sort_index()[metric]

	def get_speedup_table(self,nr_ways_fair):
		if nr_ways_fair==0:
			nr_ways_fair=1
		ipc_ref=self.original_app.properties["ipc"][nr_ways_fair]
		return 	(self.properties["ipc"]*(self.scaling_factor/ipc_ref)).sort_index()

	## Scaled ways is a python list of floats 
	## (proportional ways in cluster the application belongs to )
	def build_scaled_properties(self,scaled_ways):
		## Work with original app!!!
		float_df = pd.DataFrame()
		bdf=self.original_app.properties.drop("BENCH",axis=1)
		# Build zero interpolation
		projected_values = bdf.ix[1]+(bdf.ix[1]-bdf.ix[2])
		projected_values = projected_values.apply(lambda x: 0.000001 if x<=0 else x)
		for cur_ways in range(len(scaled_ways)):
			f_ways = scaled_ways[cur_ways]
			i_act = int(f_ways)
			i_next = i_act + 1
			frac = f_ways % 1
			if f_ways < 1:
				new_row = (1-frac)*projected_values + frac*bdf.loc[i_next]
			else: # Weighted average
				new_row = (1-frac)*bdf.loc[i_act] + frac*bdf.loc[i_next]
			new_row["float_ways"] = f_ways
			float_df=float_df.append(new_row,ignore_index=True)
		float_df.index+=1
		return float_df

def get_slowdown_vector_old(benchmark_list, solution, max_bandwidth=float('Inf')):
	throughput_value = 0
	slowdown_vector = []
	bw_alone_vector = []

	for i,benchmark in enumerate(benchmark_list):
		if i < len(solution): # i can be greater than or equal to len(solution) if the solution is partial (for bounding)
			slowdown_vector.append(benchmark.properties["slowdown"][solution[i]])
			bw_alone_vector.append(benchmark.properties["bandwidth_mbps"][solution[i]])

	## Apply BW model if necessary
	if max_bandwidth != float('Inf'):	
		## Replace slowdown vector ...
		return apply_bw_model(bw_alone_vector,slowdown_vector,max_bandwidth)
	else:
		return slowdown_vector	

def get_slowdown_vector(benchmark_list, solution, max_bandwidth=float('Inf')):
	## Apply BW model if necessary
	if max_bandwidth != float('Inf'):	
		## Replace slowdown vector ...
		return apply_bw_model_gen(benchmark_list[0:len(solution)],solution,max_bandwidth)
	else:
		slowdown_vector = []
		for i,benchmark in enumerate(benchmark_list):
			if i < len(solution): # i can be greater than or equal to len(solution) if the solution is partial (for bounding)
				slowdown_vector.append(benchmark.properties["slowdown"][solution[i]])
		return slowdown_vector	


def normalize_solution(benchmark_set, solution, remaining_ways=0):
	final_solution=list(solution)
	clustering=not(type(benchmark_set) is list)
	nr_apps_no_solution=0

	## Point to right stuff	
	if clustering:
		cluster_list=benchmark_set[1] ## Use cluster rather than that
		apps=benchmark_set[0]
	else:
		cluster_list=benchmark_set
		apps=benchmark_set

	len_partial_solution=len(solution)

	## Complete ideal solution (if necessary)
	if remaining_ways > 0:
		# Solution is partial, it is necessary to complete it (ideal solution for bounding)
		remaining_apps=len(cluster_list) - len(solution)
		for i in range(remaining_apps):
			final_solution.append(remaining_ways-(remaining_apps-1))
			if clustering:
				nr_apps_no_solution+=len(cluster_list[len_partial_solution+i])
			else:
				nr_apps_no_solution+=1

	## Obtain plain version of the solution for evaluation
	if clustering: ## Plain version of the solution
		plain_solution=[]
		for i,cluster in enumerate(cluster_list):
			cluster_ways=final_solution[i]
			for app in cluster:
				plain_solution.append(cluster_ways)

		## Overwrite variable...
		final_solution=plain_solution

	return (apps,final_solution,len(apps)-nr_apps_no_solution)
		


# Throughput function.
# Params: Benchmarks data, solution dictionary, max bandwidth and remaining ways (to know if it is a partial solution, used for bounding)
# Return: Throughput value calculated.
def throughput(benchmark_set, solution, max_bandwidth=float('Inf'), remaining_ways=0):
	##Keep track of 
	(benchmark_list,final_solution,apps_in_partial_solution)=normalize_solution(benchmark_set, solution,remaining_ways)
	slowdown_vector=get_slowdown_vector(benchmark_list,final_solution,max_bandwidth)

	## Calculate throughput
	throughput_value=0
	for slowdown in slowdown_vector:
		throughput_value += 1/slowdown

	del slowdown_vector
	del final_solution

	return throughput_value

# Unfairness function.
# Params: Benchmarks data, solution dictionary, max bandwidth and remaining ways (to know if it is a partial solution, used for bounding)
# Return: Unfairness value calculated.
def unfairness(benchmark_set, solution, max_bandwidth=float('Inf'), remaining_ways=0):
	(benchmark_list,final_solution,apps_in_partial_solution)=normalize_solution(benchmark_set, solution,remaining_ways) 
	slowdown_vector=get_slowdown_vector(benchmark_list,final_solution,max_bandwidth)

	## If remaining_ways > 0 -> theoretical unfairness does not consider 
	# greater values than the current max, nor bigger than the current min...
	## If remaining_ways for the rest, we could be making the slowdown smaller (as a result of making min smaller)

	## Pick min that cannot grow with new slowdown values...
	if remaining_ways>0:
		## Let's assume that min slowdown unfairness will not get any smaller
		unfairness_value = max(slowdown_vector) / min(slowdown_vector[0:len(solution)])
	else:
		unfairness_value = max(slowdown_vector) / min(slowdown_vector)
		
	
	del slowdown_vector
	return unfairness_value


# ANTT METRIC function.
# Params: Benchmarks data, solution dictionary, max bandwidth and remaining ways (to know if it is a partial solution, used for bounding)
# Return: Unfairness value calculated.
def antt_metric(benchmark_set, solution, max_bandwidth=float('Inf'), remaining_ways=0):
	(benchmark_list,final_solution,apps_in_partial_solution)=normalize_solution(benchmark_set, solution,remaining_ways) 
	slowdown_vector=get_slowdown_vector(benchmark_list,final_solution,max_bandwidth)

	antt_value = sum(slowdown_vector) / len(slowdown_vector)
	del slowdown_vector
	return antt_value

# Unfairness function maximizing throughput.
# Params: Benchmarks data, solution dictionary, max bandwidth and remaining ways (to know if it is a partial solution, used for bounding)
# Return: Unfairness and -throughput tuple calculated.
def unfairness_max_throughput(benchmark_set, solution, max_bandwidth=float('Inf'), remaining_ways=0):
	(benchmark_list,final_solution,apps_in_partial_solution)=normalize_solution(benchmark_set, solution,remaining_ways) 
	throughput_value = 0	
	
	slowdown_vector=get_slowdown_vector(benchmark_list,final_solution,max_bandwidth)

	## Calculate throughput
	for slowdown in slowdown_vector:
		throughput_value += 1/slowdown		

	## Pick min that cannot grow with new slowdown values...
	if remaining_ways>0:
		## Let's assume that min slowdown unfairness will not get any smaller
		unfairness_value = max(slowdown_vector) / min(slowdown_vector[0:len(solution)])
	else:
		unfairness_value = max(slowdown_vector) / min(slowdown_vector)

	del slowdown_vector
	return (unfairness_value,-throughput_value)


def fairness_metric(benchmark_set, solution, max_bandwidth=float('Inf'), remaining_ways=0):
	(benchmark_list,final_solution,apps_in_partial_solution)=normalize_solution(benchmark_set, solution,remaining_ways) 
	slowdown_vector=get_slowdown_vector(benchmark_list,final_solution,max_bandwidth)		
	fairness=1-np.std(slowdown_vector)/np.mean(slowdown_vector)
	del slowdown_vector
	return fairness


################################################################################
##                            FILE READING FUNCTIONS                          ##
################################################################################

# Constructs a memory structure of the useful data of CSV files
# Params: Path to benchmark summary CSV file and workloads CSV file. Separator for summary file
# Return: Structure with all useful data loaded.
def get_workloads_table_from_csv_gen(summary_csv_path, workloads_csv_path,separator=None):
	# Read workloads_csv and save it to structure in memory to increase performance
	workloads_table = []
	workloads = []
	workloads_csv_file = open(workloads_csv_path)
	workloads_csv_entry = csv.reader(workloads_csv_file)

	for workload in workloads_csv_entry:
		workloads.append(workload)
		workloads_table.append([]) ## Empty list of benchmarks or applications

	workloads_csv_file.close()
	del workloads_csv_file

	# At this point workloads and workloads_table structures are created successfully
	# Now, read summary stats and save benchmark data from benchmarks in workloads to structure in memory (workloads_table)
	
	# Get DataFrame with all CSV data
	if separator:
		summary_dframe = pd.read_csv(summary_csv_path,sep=separator)
	else:
		summary_dframe = pd.read_csv(summary_csv_path,sep='\t')
	
	# Get list of benchmarks
	benchmarks_list = summary_dframe["BENCH"].unique()
	benchmark_dict={}

	# Get info for each benchmark in the global table
	for benchmark_name in benchmarks_list:
		# Get benchmark_name's DataFrame using NR_WAYS as index
		benchmark_dframe = summary_dframe.loc[summary_dframe["BENCH"] == benchmark_name].set_index("NR_WAYS")
		
		## Calculate max IPS and slowdown and save on benchmark_dframe column
		ipc_vals=benchmark_dframe[["ipc"]].values
		ipc_max = max(ipc_vals)[0]
		benchmark_dframe["slowdown"] = ipc_max/benchmark_dframe["ipc"]
		nr_ways=len(ipc_vals)

		# Keep only useful data
		# THIS HAS BEEN REMOVED TO MAKE THE CODE WORK WITH ANY KIND OF EVENTS...
		#benchmark_dframe = benchmark_dframe[["bandwidth_mbps","llcrpki","llcmpki","slowdown"]]

		benchmark_dict[benchmark_name]=benchmark_dframe


	## Traverse workloads first
	i=0
	for workload in workloads:
		for benchmark_name in workload:
			if not benchmark_name in benchmarks_list:
				print "Can't find benchmark name (%s) in global table" % benchmark_name
				return None
			## Get right Dframe
			benchmark_dframe=benchmark_dict[benchmark_name]
			## Safe to continue	
			workloads_table[i].append(App(benchmark_name,benchmark_dframe))
		i += 1		

	return (workloads_table,nr_ways)


def get_workloads_table_from_csv(summary_csv_path, workloads_csv_path,separator=None):
	return get_workloads_table_from_csv_gen(summary_csv_path, workloads_csv_path,separator)[0]

# Constructs a memory structure of the useful data of CSV files
# Params: Path to benchmark summary CSV file and list of workloads
# Return: Structure with all useful data loaded.
def get_workloads_table_from_list_gen(summary_csv_path, workloads, separator=None):
	# Read workloads_csv and save it to structure in memory to increase performance
	workloads_table = []

	for workload in workloads:
		workloads_table.append([]) ## Empty list of benchmarks or applications

	# At this point workloads and workloads_table structures are created successfully
	# Now, read summary stats and save benchmark data from benchmarks in workloads to structure in memory (workloads_table)
	
	# Get DataFrame with all CSV data
	if separator:
		summary_dframe = pd.read_csv(summary_csv_path,sep=separator)
	else:
		summary_dframe = pd.read_csv(summary_csv_path)
	
	# Get list of benchmarks
	benchmarks_list = summary_dframe["BENCH"].unique()
	benchmark_dict={}

	# Get info for each benchmark in the global table
	for benchmark_name in benchmarks_list:
		# Get benchmark_name's DataFrame using NR_WAYS as index
		benchmark_dframe = summary_dframe.loc[summary_dframe["BENCH"] == benchmark_name].set_index("NR_WAYS")
		
		## Calculate max IPS and slowdown and save on benchmark_dframe column
		ipc_vals=benchmark_dframe[["ipc"]].values
		ipc_max = max(ipc_vals)[0]
		nr_ways=len(ipc_vals)
		benchmark_dframe["slowdown"] = ipc_max/benchmark_dframe["ipc"]


		# Keep only useful data
		# THIS HAS BEEN REMOVED TO MAKE THE CODE WORK WITH ANY KIND OF EVENTS...
		#benchmark_dframe = benchmark_dframe[["bandwidth_mbps","llcrpki","llcmpki","slowdown"]]

		benchmark_dict[benchmark_name]=benchmark_dframe


	## Traverse workloads first
	i=0
	for workload in workloads:
		for benchmark_name in workload:
			if not benchmark_name in benchmarks_list:
				print "Can't find benchmark name (%s) in global table" % benchmark_name
				return None
			## Get right Dframe
			benchmark_dframe=benchmark_dict[benchmark_name]
			## Safe to continue	
			workloads_table[i].append(App(benchmark_name,benchmark_dframe))
		i += 1		

	return (workloads_table,nr_ways)


def get_workloads_table_from_list(summary_csv_path, workloads, separator=None):
		return get_workloads_table_from_list_gen(summary_csv_path, workloads)[0]

################################################################################
##                               PRIVATE FUNCTIONS                            ##
################################################################################

# Private recursive function that calculates best solution for a apps data table using a cost function.
# Params: Apps data dictionary, cost function, operator to compare, slots available, list of candidate slots to assign, 
# 	, partial solution increased in each recursive call and activate bound or not.
# Return: Best cost, best solution and explored solutions.
def get_optimal_schedule_aux(benchmark_set, cost_function, max_bandwidth, default_cost, op, nr_ways, partial_solution, bound):
	
	if type(benchmark_set) is list:
		cluster_list=benchmark_set
	else:
		cluster_list=benchmark_set[1] ## Use cluster rather than that

	n_clusters=len(cluster_list)
	level=len(partial_solution)
	#cluster_selected = cluster_list[level]
	

	if level == n_clusters - 1: # Trivial Base case
		partial_solution.append(nr_ways) # Being the last app, we give you all the available 'slots'
		best_cost = cost_function(benchmark_set, partial_solution, max_bandwidth)
		best_solution = partial_solution
		total_branches=1
	elif nr_ways == n_clusters - level: ## Optimization, one way for each app....
		##Easy ...
		for i in range(nr_ways):
			partial_solution.append(1)
		best_cost = cost_function(benchmark_set, partial_solution, max_bandwidth)
		best_solution = partial_solution
		total_branches = 1	
	else: # Recursive case
		## Note that the cost value can be a tuple!!
		best_cost = default_cost 
		best_solution = None
		total_branches = 0 
		#0.0 if op == operator.gt else float('Inf')
		
		for nr_ways_assigned in range(1, nr_ways - (n_clusters - level - 1) + 1):
			new_partial_solution = list(partial_solution)
			new_partial_solution.append(nr_ways_assigned)
			remaining_ways = nr_ways - nr_ways_assigned

			ideal_cost=cost_function(benchmark_set, new_partial_solution, max_bandwidth, remaining_ways)

			## Bound
			if bound and op(best_cost, ideal_cost):
				total_branches+=1 ## Pruning counts
				continue
				
			(cost, solution, branches) = get_optimal_schedule_aux(benchmark_set, cost_function, max_bandwidth, best_cost, op, 
				remaining_ways, new_partial_solution, bound)

			total_branches += branches
			
			## Good solution so far is worse than default value	
			if solution is None:
				continue

			if op(cost, best_cost):
				best_cost = cost
				best_solution = solution
			else:
				del solution

	return (best_cost, best_solution, total_branches)

# Marginal utility function
def marginal_utility(from_way, to_way, mpki_from, mpki_to):
	return (mpki_to - mpki_from) / (from_way - to_way)


# Private function that calculates maximal marginal used by lookahead algorithm.
# Params: Apps data dictionary, benchmark name, slots assigned to benchmark, available ways and total slots.
# Return: Maximal marginal value and slots to get this value for that benchmark.
def get_max_mu_gen(curve, current_benchmark_ways, available_ways, max_nr_ways):
	current_benchmark_misses = curve[current_benchmark_ways]
	max_mu = 0
	mu_ways = -1
	exploration_ways = min(available_ways, max_nr_ways - current_benchmark_ways) 

	for i in range(1, exploration_ways + 1):
		new_assigned_ways = current_benchmark_ways + i
		new_misses = curve[new_assigned_ways]
		mu = marginal_utility(current_benchmark_ways, new_assigned_ways, current_benchmark_misses, new_misses)
		
		if mu > max_mu:
			max_mu = mu
			mu_ways = new_assigned_ways

	return (max_mu, mu_ways)

# Private function that calculates maximal marginal used by lookahead algorithm.
# Params: Apps data dictionary, benchmark name, slots assigned to benchmark, available ways and total slots.
# Return: Maximal marginal value and slots to get this value for that benchmark.
def get_max_mu(benchmark, current_benchmark_ways, available_ways, max_nr_ways, use_bandwidth=False):
	if use_bandwidth:
		llcmpki_table=benchmark.properties["bandwidth_mbps"]
	else:
		llcmpki_table=benchmark.properties["llcmpki"]
	current_benchmark_misses = llcmpki_table[current_benchmark_ways]
	max_mu = 0
	mu_ways = -1
	exploration_ways = min(available_ways, max_nr_ways - current_benchmark_ways) 

	for i in range(1, exploration_ways + 1):
		new_assigned_ways = current_benchmark_ways + i
		new_misses = llcmpki_table[new_assigned_ways]
		mu = marginal_utility(current_benchmark_ways, new_assigned_ways, current_benchmark_misses, new_misses)
		
		if mu > max_mu:
			max_mu = mu
			mu_ways = new_assigned_ways

	return (max_mu, mu_ways)


# Same parameters as lookahead_algorithm_gen but using an array of curves as a parameter
def lookahead_algorithm_gen(curves, nr_ways):
	solution = []
	available_ways = nr_ways
	there_is_maximal_marginal = True
	
	# Each app needs at least one 'slot'
	for curve in curves:
		solution.append(1)
		available_ways -= 1

	while available_ways > 0 and there_is_maximal_marginal:
		global_max_mu = 0
		global_mu_ways= -1

		# Find maximal marginal utility across applications 
		for idx,curve in enumerate(curves):
			(max_mu, mu_ways) = get_max_mu_gen(curve, solution[idx], available_ways, nr_ways)
			if max_mu > global_max_mu:
				global_max_mu = max_mu
				selected_idx = idx
				global_mu_ways = mu_ways

		if global_max_mu == 0:
			there_is_maximal_marginal = False
		else:
			# Update partitions for selected benchmark and update available ways
			available_ways -= global_mu_ways - solution[selected_idx]
			solution[selected_idx] = global_mu_ways

	# If there are available ways yet, allocate remaining ways fairly
	if available_ways > 0:
		n_apps=len(curves)
		nr_fair_ways_per_benchmark = available_ways / n_apps
		remaining_ways = available_ways % n_apps
		i = 0

		while i < n_apps and available_ways > 0:
			nr_extra_ways = nr_fair_ways_per_benchmark
			
			if remaining_ways > 0:
				nr_extra_ways += 1
				remaining_ways -= 1
			
			solution[i] += nr_extra_ways
			available_ways -= nr_extra_ways		
			i += 1 
		
	return solution


def number_of_solutions_partitioning_gen(nr_ways,nr_partitions,count_intermediate):
	assert nr_partitions<=nr_ways

	if nr_partitions==1 or nr_ways==nr_partitions:
		return 1		
	elif nr_ways==nr_partitions+1: ## Just one app can take the extra way...
		return nr_partitions
	else:
		## Recursive case
		total_solutions=1 if count_intermediate else 0

		for i in range(1,nr_ways-(nr_partitions-1)+1):
			total_solutions+=number_of_solutions_partitioning_gen(nr_ways-i,nr_partitions-1,count_intermediate)

		return total_solutions 

def repr_node(node):
	out_str = '"['
	for s in node:
		out_str += str(s) + ","
	out_str=out_str.rstrip(",")
	out_str += ']"'
	return out_str

def generate_dot_tree(nr_ways,nr_partitions,edges,padre):
	assert nr_partitions<=nr_ways

	if nr_partitions==1:
		node=list(padre)
		node.append(nr_ways)
		edges.append(repr_node(padre)+' -> '+repr_node(node))
		return 1
	elif nr_ways==nr_partitions: ## Just one app can take the extra way...
		node=list(padre)
		for i in range(nr_partitions):
			node.append(1)
		edges.append(repr_node(padre)+' -> '+repr_node(node))
		return nr_partitions
	else:
		## Recursive case
		total_solutions = 0
		for i in range(1,nr_ways-(nr_partitions-1)+1):
			node=list(padre)
			node.append(i)
			if padre:
				edges.append(repr_node(padre)+' -> '+repr_node(node))
			else:
				edges.append('"[]"'+' -> '+repr_node(node))

			total_solutions+=generate_dot_tree(nr_ways-i,nr_partitions-1,edges,node)

		return total_solutions

def write_dot_tree(nr_ways,n_apps):
	edges=[]
	generate_dot_tree(nr_ways,n_apps,edges,[])
	fd = open("tree.dot","w+")
	fd.write('digraph {\n')
	fd.write('node [fontname="courier", shape="plaintext" fontsize=24,  margin=0, width=0, height=0 , fontcolor=black];\n')

	for e in edges:
		fd.write(e+'\n')
	fd.write('}\n')
	fd.close()

def number_of_solutions_partitioning(nr_ways,nr_partitions):
	return number_of_solutions_partitioning_gen(nr_ways,nr_partitions,False)



## Bell triangle
##https://en.wikipedia.org/wiki/Bell_number
def number_of_solutions_clustering(nr_items_in_set):
	sol=[[1]]

#	print [1]

	for i in range(1,nr_items_in_set):
		prev_row=sol[i-1]
		new_row=[prev_row[-1]]
		prev_val=prev_row[-1]
		for j in range(0,len(sol[i-1])):
			new_val=prev_val+prev_row[j]
			new_row.append(new_val)
			prev_val=new_val

		sol.append(new_row)
#		print new_row

	return sol[nr_items_in_set-1][-1]

def number_of_solutions_partitioning_dp_gen(nr_ways,nr_partitions,count_intermediate):
	assert nr_partitions<=nr_ways

	sol=[]
	## Empty matrices
	for i in range(0,nr_ways):
		sol.append([0]*nr_partitions)


	for w in range(1,nr_ways+1):
		for p in range(1,nr_partitions+1):
			if p>w:
				pass ## Does not make sense
			elif p==1 or p==w:
				sol[w-1][p-1]=1
			elif w==p+1:
				sol[w-1][p-1]=p
			else:
				## Recursive case
				total_solutions=1 if count_intermediate else 0
				for i in range(1,w-(p-1)+1):
					total_solutions+=sol[w-i-1][p-2]
				sol[w-1][p-1]=total_solutions
	return sol


def number_of_solutions_partitioning_dp(nr_ways,nr_partitions):
	return number_of_solutions_partitioning_dp_gen(nr_ways,nr_partitions,False)

## https://stackoverflow.com/questions/19368375/set-partitions-in-python
def generate_possible_clusters_nofilter(collection):
	if len(collection) == 1:
		yield [ collection ]
		return

	first = collection[0]
	for smaller in generate_possible_clusters_nofilter(collection[1:]):
	# insert `first` in each of the subpartition's subsets
		for n, subset in enumerate(smaller):
			yield smaller[:n] + [[ first ] + subset]  + smaller[n+1:]
		# put `first` in its own subset 
		yield [[ first ]] + smaller


def generate_possible_clusters(collection,nr_ways):
	if len(collection) == 1:
		yield [ collection ]
		return

	first = collection[0]
	for smaller in generate_possible_clusters(collection[1:],nr_ways):
		size_clustering=len(smaller)

		##Unfeasible solution
		if size_clustering>nr_ways:
			continue

		# insert `first` in each of the subpartition's subsets
		# This does not make the clustering any bigger
		for n, subset in enumerate(smaller):
			yield smaller[:n] + [[ first ] + subset]  + smaller[n+1:]
		
		# put `first` in its own subset if that does not make it grow
		if size_clustering<nr_ways:
			yield [[ first ]] + smaller

def determine_number_of_cluster_cache_partitioning(nr_ways,nr_items_in_set):
	collection=[i for i in range(nr_items_in_set)] ## generate numbers from 0 to N-1
	bound=nr_items_in_set if nr_ways>=nr_items_in_set else nr_ways
	sols_cp=number_of_solutions_partitioning_dp(nr_ways,bound)

	# For each solution clustering determine the number of possible cache partitionings
	possible_clusters=generate_possible_clusters(collection,nr_ways)
	nr_options=0

	for clustering in possible_clusters:
		nr_options+=sols_cp[nr_ways-1][len(clustering)-1]

	return nr_options	


# Private function that implements lookahead algorithm.
# Params: Apps data dictionary and total slots.
# Return: Partitioning solution for that workload.
def lookahead_algorithm(workload, nr_ways, use_bandwidth=False):
	solution = []
	available_ways = nr_ways
	there_is_maximal_marginal = True
	
	# Each app needs at least one 'slot'
	for benchmark in workload:
		solution.append(1)
		available_ways -= 1

	while available_ways > 0 and there_is_maximal_marginal:
		global_max_mu = 0
		global_mu_ways= -1

		# Find maximal marginal utility across applications 
		for idx,benchmark in enumerate(workload):
			(max_mu, mu_ways) = get_max_mu(benchmark, solution[idx], available_ways, nr_ways, use_bandwidth)
			if max_mu > global_max_mu:
				global_max_mu = max_mu
				selected_idx = idx
				global_mu_ways = mu_ways

		if global_max_mu == 0:
			there_is_maximal_marginal = False
		else:
			# Update partitions for selected benchmark and update available ways
			available_ways -= global_mu_ways - solution[selected_idx]
			solution[selected_idx] = global_mu_ways

	# If there are available ways yet, allocate remaining ways fairly
	if available_ways > 0:
		n_apps=len(workload)
		nr_fair_ways_per_benchmark = available_ways / n_apps
		remaining_ways = available_ways % n_apps
		i = 0

		while i < n_apps and available_ways > 0:
			nr_extra_ways = nr_fair_ways_per_benchmark
			
			if remaining_ways > 0:
				nr_extra_ways += 1
				remaining_ways -= 1
			
			solution[i] += nr_extra_ways
			available_ways -= nr_extra_ways		
			i += 1 
		
	return solution

# Private function that builds bandwidth gap between X nr_way and X+1 nr_way, for each benchmark.
# Params: App vector and total slots
# Return: BW matrix [NBENCHS][NR_WAYS].
def build_bandwidth_reduction_table(workload, nr_ways):
	bandwidth_reduction_table = []
	
	for benchmark in workload:
		bench_table=[]
		bw_tbl=benchmark.properties["bandwidth_mbps"]
		bandwidth_ref = bw_tbl[1]
		
		for i in range(1, nr_ways + 1):
			bench_table.append(bandwidth_ref - bw_tbl[i])

		bandwidth_reduction_table.append(bench_table)

	return bandwidth_reduction_table
			
# Private function that implements Yu-Petrof algorithm.
# Params: Bandwidth reduction matrix, total slots and optionally, a threshold.
# Return: Partitioning solution for that workload.
def yu_petrof_algorithm(workload, bandwidth_reduction_table, nr_ways, threshold = 0):
	solution = []
	nr_apps=len(bandwidth_reduction_table)
	available_ways = nr_ways
	there_is_maximal_value = True
		
	# Each app needs at least one 'slot'
	for i in range(nr_apps):
		solution.append(1)
		available_ways -= 1

	# Get bandwidth gap between 'slots' for each app
	inc_bandwidth_reduction_table = []
	for idx_benchmark in range(nr_apps):
		inc_bandwidth_reduction_values = []
		for i in range(nr_ways - 1):
			inc_bandwidth_reduction_values.append(bandwidth_reduction_table[idx_benchmark][i + 1] - bandwidth_reduction_table[idx_benchmark][i])
		inc_bandwidth_reduction_table.append( inc_bandwidth_reduction_values)

	while available_ways > 0 and there_is_maximal_value:
		maximal_value = 0
		selected_benchmark_idx = -1
		ways_for_selected_benchmark = -1
		
		for idx_benchmark in range(nr_apps):
			for i in range(available_ways):
				benefit = inc_bandwidth_reduction_table[idx_benchmark][i]
				if (maximal_value == 0 and benefit > 0) or (benefit > (maximal_value + threshold)):
					maximal_value = benefit
					selected_benchmark_idx = idx_benchmark
					ways_for_selected_benchmark = i + 1
		
		if maximal_value == 0:
			there_is_maximal_value = False
		else:
			# Remove values from bandwidth_reduction_table, update partitions for selected benchmark and update available ways
			for i in range(ways_for_selected_benchmark):
				inc_bandwidth_reduction_table[selected_benchmark_idx][i] = 0
			available_ways -= ways_for_selected_benchmark - solution[selected_benchmark_idx]
			solution[selected_benchmark_idx] = ways_for_selected_benchmark

	# If there are available ways yet, allocate remaining ways fairly
	if available_ways > 0:
		nr_fair_ways_per_benchmark = available_ways / nr_apps
		remaining_ways = available_ways % nr_apps
		i = 0

		while i < nr_apps and available_ways > 0:
			nr_extra_ways = nr_fair_ways_per_benchmark
			
			if remaining_ways > 0:
				nr_extra_ways += 1
				remaining_ways -= 1
			
			solution[i] += nr_extra_ways
			available_ways -= nr_extra_ways		
			i += 1 

	return solution
	
################################################################################
##                               PUBLIC FUNCTIONS                             ##
################################################################################

## Get array of partition masks (once this is done)
def get_partition_masks(partition_sizes):
	next_available_way=0
	masks=[]
	for size in partition_sizes:
		masks.append(hex(((1<<size)-1)<<next_available_way))
		next_available_way=next_available_way+size
	return masks

##
# Return equal LLC partitioning 
def get_equal_llc_schedule(nr_partitions, nr_ways):
	nr_fair_ways=nr_ways/nr_partitions
	nr_extra_ways=nr_ways % nr_partitions

	solution=[]

	for i in range(nr_partitions):
		if nr_extra_ways>0:
			ways_this_app=nr_fair_ways+1
			nr_extra_ways=nr_extra_ways-1		
		else:
			ways_this_app=nr_fair_ways

		solution.append(ways_this_app)
	return solution

##
# Return simple custom partitioning 
def get_simple_schedule(workload, nr_ways, req_threshold):
	solution = []
	nr_apps = len(workload)
	nr_harmless_apps = 0
	
	for llcrpki in map(lambda app:app.properties["llcrpki"][nr_ways], workload):
		if llcrpki < req_threshold:
			solution.append(1)
			nr_harmless_apps += 1
		else:
			solution.append(0)

	nr_available_ways = nr_ways - nr_harmless_apps
	nr_harm_apps = nr_apps - nr_harmless_apps
	nr_fair_ways = nr_available_ways / nr_harm_apps
	nr_extra_ways = nr_available_ways % nr_harm_apps

	j = 0
	
	for i in range(nr_harm_apps):
		if nr_extra_ways > 0:
			ways_this_app = nr_fair_ways+1
			nr_extra_ways = nr_extra_ways-1		
		else:
			ways_this_app = nr_fair_ways
		
		while solution[j] != 0:
			j += 1
		
		solution[j] = ways_this_app

	return solution

##
# Return simple 2 custom partitioning 
def get_simple2_schedule(workload, nr_ways, req_threshold, harmless_ratio):
	solution = []
	nr_apps = len(workload)
	nr_harmless_apps = 0
	nr_harmless_ways = int(( nr_ways / nr_apps ) * harmless_ratio)

	if nr_harmless_ways < 1:
		nr_harmless_ways = 1
	
	for llcrpki in map(lambda app:app.properties["llcrpki"][nr_ways], workload):
		if llcrpki < req_threshold:
			solution.append(nr_harmless_ways)
			nr_harmless_apps += 1
		else:
			solution.append(0)

	nr_available_ways = nr_ways - ( nr_harmless_apps * nr_harmless_ways )
	nr_harm_apps = nr_apps - nr_harmless_apps
	nr_fair_ways = nr_available_ways / nr_harm_apps
	nr_extra_ways = nr_available_ways % nr_harm_apps

	j = 0
	
	for i in range(nr_harm_apps):
		if nr_extra_ways > 0:
			ways_this_app = nr_fair_ways+1
			nr_extra_ways = nr_extra_ways-1		
		else:
			ways_this_app = nr_fair_ways
		
		while solution[j] != 0:
			j += 1
		
		solution[j] = ways_this_app

	return solution


##
# Return smartfake custom partitioning 
def get_smartfake_schedule(workload, nr_ways, req_threshold, misses_gap_ratio_threshold):
	solution = []
	nr_apps = len(workload)
	nr_harmless_apps = 0

	## HARMLESS APPS SECTION
	# Search harmless apps and assign 1 way for them
	for llcrpki in map(lambda app:app.properties["llcrpki"][nr_ways], workload):
		if llcrpki < req_threshold:
			solution.append(1)
			nr_harmless_apps += 1
		else:
			solution.append(0)

	## HARM APPS SECTION
	critical_working_set_list = []
	ways_harm_apps_list = []
	
	# Calculate critical working set for harm apps
	for i in range(nr_apps):
		if solution[i] == 0: # If app is a harm app, calculates his working set
			max_misses_gap = workload[i].properties["llcmpki"][1] - workload[i].properties["llcmpki"][2]
			ways_in_working_set = 2

			while ways_in_working_set < nr_ways:
				current_misses_gap = workload[i].properties["llcmpki"][ways_in_working_set] - workload[i].properties["llcmpki"][ways_in_working_set+1]
				misses_gap_ratio = current_misses_gap / max_misses_gap

				if misses_gap_ratio <= misses_gap_ratio_threshold:
					break

				ways_in_working_set += 1

			critical_working_set_list.append(ways_in_working_set)
				
	
	nr_available_ways = nr_ways - nr_harmless_apps
	nr_harm_apps = nr_apps - nr_harmless_apps

	# Calculate ways to assign for each harm app
	total_weigh = sum(critical_working_set_list)

	for i in range(nr_harm_apps):
		fair_ways_for_this_app = max(1, int((critical_working_set_list[i] * nr_available_ways) / total_weigh))
		ways_harm_apps_list.append(fair_ways_for_this_app)
		nr_available_ways -= fair_ways_for_this_app
		
	# Assign ways to harm apps, distributing the existing available ways
	j = 0
	
	for i in range(nr_harm_apps):
		if nr_available_ways > 0:
			if nr_available_ways > (nr_harm_apps - i):
				ways_this_app = ways_harm_apps_list[i] + 2
				nr_available_ways -= 2
			else:
				ways_this_app = ways_harm_apps_list[i] + 1
				nr_available_ways -= 1
		else:
			ways_this_app = ways_harm_apps_list[i]
		
		while solution[j] != 0:
			j += 1
		
		solution[j] = ways_this_app

	return solution


##
# Return smartfake custom partitioning 
def get_smartfake2_schedule(workload, nr_ways, req_threshold, harmless_ratio, misses_gap_ratio_threshold):
	solution = []
	nr_apps = len(workload)
	nr_harmless_apps = 0
	nr_harmless_ways = int(( nr_ways / nr_apps ) * harmless_ratio)

	if nr_harmless_ways < 1:
		nr_harmless_ways = 1

	## HARMLESS APPS SECTION
	# Search harmless apps and assign nr_harmless_ways for them
	for llcrpki in map(lambda app:app.properties["llcrpki"][nr_ways], workload):
		if llcrpki < req_threshold:
			solution.append(nr_harmless_ways)
			nr_harmless_apps += 1
		else:
			solution.append(0)

	## HARM APPS SECTION
	critical_working_set_list = []
	ways_harm_apps_list = []
	
	# Calculate critical working set for harm apps
	for i in range(nr_apps):
		if solution[i] == 0: # If app is a harm app, calculates his working set
			max_misses_gap = workload[i].properties["llcmpki"][1] - workload[i].properties["llcmpki"][2]
			ways_in_working_set = 2

			while ways_in_working_set < nr_ways:
				current_misses_gap = workload[i].properties["llcmpki"][ways_in_working_set] - workload[i].properties["llcmpki"][ways_in_working_set+1]
				misses_gap_ratio = current_misses_gap / max_misses_gap

				if misses_gap_ratio <= misses_gap_ratio_threshold:
					break

				ways_in_working_set += 1

			critical_working_set_list.append(ways_in_working_set)
				
	
	nr_available_ways = nr_ways - (nr_harmless_apps * nr_harmless_ways)
	nr_harm_apps = nr_apps - nr_harmless_apps

	# Calculate ways to assign for each harm app
	total_weigh = sum(critical_working_set_list)

	for i in range(nr_harm_apps):
		fair_ways_for_this_app = max(1, int((critical_working_set_list[i] * nr_available_ways) / total_weigh))
		ways_harm_apps_list.append(fair_ways_for_this_app)
		nr_available_ways -= fair_ways_for_this_app
		
	# Assign ways to harm apps, distributing the existing available ways
	j = 0
	
	for i in range(nr_harm_apps):
		if nr_available_ways > 0:
			if nr_available_ways > (nr_harm_apps - i):
				ways_this_app = ways_harm_apps_list[i] + 2
				nr_available_ways -= 2
			else:
				ways_this_app = ways_harm_apps_list[i] + 1
				nr_available_ways -= 1
		else:
			ways_this_app = ways_harm_apps_list[i]
		
		while solution[j] != 0:
			j += 1
		
		solution[j] = ways_this_app

	return solution


##
# Calculate LLC partitioning based on rate of demand (BW)
def on_demand_partitioning(workload,nr_ways,max_bandwidth=float('Inf')):
	bw_alone_vector=map(lambda app:app.properties["bandwidth_mbps"][nr_ways],workload)

	##Reserve one way for each benchmark
	nr_apps=len(workload)
	ways_to_share=nr_ways-nr_apps

	## Apply BW model if necessary
	if max_bandwidth != float('Inf'):	
		## Replace slowdown vector ...
		total_bw,bw_shared=determine_bw_shared(bw_alone_vector,max_bandwidth)
	else:
		total_bw=sum(bw_alone_vector)
		bw_shared=bw_alone_vector

	solution=map(lambda x: int(ways_to_share*(x/total_bw))+1,bw_shared) ## Plus one -> at least one each
	nr_assigned_ways=sum(solution)
	nr_remaining_ways=nr_ways-nr_assigned_ways

	## Assign remaining ways in a RR order
	while nr_remaining_ways>0:
		i=0
		while nr_remaining_ways>0 and i<nr_apps:
			solution[i]=solution[i]+1
			nr_remaining_ways=nr_remaining_ways-1
			i=i+1

	return solution


def generate_solutions_to_explore(nsols,nr_ways,nr_apps,threshold,partial_sol):
	initial_solution=[]
	for i in range(1, nr_ways - (nr_apps - 1) + 1):
		subsol=partial_sol+[i]
		remaining_ways=nr_ways-i
		remaining_apps=nr_apps-1
		nr_sols=nsols[remaining_ways-1][remaining_apps-1]
		if len(subsol)>=2 and nr_sols<=threshold:
			yield (subsol,nr_sols)
		else:
			for x in generate_solutions_to_explore(nsols,remaining_ways,remaining_apps,threshold,subsol):
				yield x


## Just one solution
def unroll_solution_cluster(nsols,nr_ways,nr_apps,threshold,sol):
	## Calculate number of solutions
	remaining_ways=nr_ways-sum(sol)
	remaining_apps=nr_apps-len(sol)
	nr_sols=nsols[remaining_ways-1][remaining_apps-1]
	if nr_sols<=threshold:
		return [sol]
	else:
		sols=[]
		iterator=generate_solutions_to_explore(nsols,remaining_ways,remaining_apps,threshold,sol)
		for x,siz in iterator:
			sols.append(x)

		return sols


def unroll_solutions(nsols,nr_ways,nr_apps,threshold,clusters):	
	rclusters=[]
	for cluster in clusters:
		for sol in cluster:
			sols=unroll_solution_cluster(nsols,nr_ways,nr_apps,threshold,sol)
			single_sols=[[x] for x in sols]
			rclusters.extend(single_sols)
	return rclusters


def fusion_tasks(partial_solutions,low_threshold,high_threshold,disable=False):
	partsols=[x for x in partial_solutions]
	## Copy for comparison purposes
	unmerged_sols=list(partsols)
	if disable:
		return ([ [x] for (x,y) in partsols],unmerged_sols)


	merged_solutions=[]

	## Traverse iterator
	
	remaining_solutions=len(partsols)

	while remaining_solutions>0:
		this_cluster=[]
		work=0
		
		## Merge as much as possible
		i=0
		while (i<remaining_solutions) and work<low_threshold:
			(psol,size)=partsols[i]

			##Try to add to this cluster without exceeding too much
			if work+size<=high_threshold:
				this_cluster.append(psol) # (psol,size))
				work=work+size
				partsols.pop(i)
				remaining_solutions=remaining_solutions-1
			elif size>=low_threshold:
				## If we find a big enough one, it goes to the cluster
				merged_solutions.append([psol]) #[(psol,size)])
				remaining_solutions=remaining_solutions-1
				partsols.pop(i)

			else: ## Continue
				i=i+1

		## Add cluster to merged_solutions
		merged_solutions.append(this_cluster)
		
	return (merged_solutions,unmerged_sols)



## Single-workload pruning
def subsol_is_promising(benchmark_set, cost_function, max_bandwidth, default_cost, op, nr_ways, partial_solution,**kwargs):

	if type(partial_solution) is tuple:
		partial_solution=partial_solution[0] ## Discard size...

	remaining_ways = nr_ways - sum(partial_solution)
	new_partial_solution = list(partial_solution)
	ideal_cost=cost_function(benchmark_set, new_partial_solution, max_bandwidth, remaining_ways)


	return not op(default_cost, ideal_cost)


## Wrappper for avoiding copy of parameters
def subsol_is_promising_parallel(default_cost,partial_solution):

	params=get_global_properties()
	params["default_cost"]=default_cost
	params["partial_solution"]=partial_solution

	return subsol_is_promising(**params)

## Multi-workload candidate filtering
def subsols_are_promising(benchmark_set, cost_function, max_bandwidth, default_cost, op, nr_ways, partial_solutions):
	return [ subsol_is_promising(benchmark_set, cost_function, max_bandwidth, default_cost, op, nr_ways, x) for  x in partial_solutions]


## function to set up global data
gbl_dict={}
#@dview.remote(block=True):
def set_global_properties(key_vals):
	import os
	global gbl_dict
	pid=os.getpid()
	gbl_dict.update(key_vals)
	## Establish BW model
	if "bw_model" in key_vals:
		select_bw_model(key_vals["bw_model"])

	#ttps://stackoverflow.com/questions/3061/calling-a-function-of-a-module-by-using-its-name-a-string
	#if "cost_function" in key_vals:
		#retrieve function symbol

	return pid	


def get_global_properties(foo=None):
	global gbl_dict
	return gbl_dict #[x.name for x in gbl_dict["benchmark_set"]]

def get_global_property(foo=None):
	return "BW model is %d" % glb_active_bw_model

## Generalized version of the branch and bound algorithm that processes a list of partial solutions and returns the best of all of them
def get_optimal_schedule_aux_multi(benchmark_set, cost_function, max_bandwidth, default_cost, op, nr_ways, partial_solutions, bound,**kwargs):
	total_branches=0
	best_solution=None
	best_cost=default_cost

	kwargs.setdefault("prepruning",True)
	prepruning=kwargs["prepruning"]

	## Get the best of each and every partial solution
	for it,partial_solution in enumerate(partial_solutions):
		## Initial pruning
		remaining_ways = nr_ways - sum(partial_solution)
		new_partial_solution = list(partial_solution)

		if prepruning:
			ideal_cost=cost_function(benchmark_set, partial_solution, max_bandwidth, remaining_ways)
		
			## Try to prune if bound enabled
			if bound and op(best_cost, ideal_cost):
				total_branches=total_branches+1 ## Pruning kind of counts
				#res.append("prune %d" % it)
				continue

		(cost, solution, branches) = get_optimal_schedule_aux(benchmark_set, cost_function, max_bandwidth, best_cost, op, remaining_ways, new_partial_solution, bound)

		#res.append("sol %s it %d with cost %d" % (str(solution),it,cost))

		total_branches += branches

		## Good solution so far is worse than default value	
		if not solution:
			#res.append("sol it %d is None" % it)
			continue

		## Replace best solution
		if op(cost, best_cost):
			best_cost = cost
			best_solution = list(solution)
			#res.append("sol it %d is better" % cost)
		else:
			del solution		

	return (best_cost, best_solution, total_branches) #,res)

def get_optimal_schedule_aux_multi_parallel(default_cost,partial_solutions,prepruning):
	params=get_global_properties()

	params["default_cost"]=default_cost
	params["partial_solutions"]=partial_solutions
	params["prepruning"]=prepruning

	## Just in case...
	if "partial_solution" in params:
		del params["partial_solution"] 

	return get_optimal_schedule_aux_multi(**params)


def get_start_end_micros(start,end,start_datet):
	# Get sequential time execution
	end_seqtime=pytz.utc.localize(end)
	end_seqtime=end_seqtime.astimezone(pytz.timezone("utc"))
	start_seqtime=pytz.utc.localize(start)
	start_seqtime=start_seqtime.astimezone(pytz.timezone("utc"))

	start_micros=int((start_seqtime-start_datet).total_seconds()*1000000)
	end_micros=int((end_seqtime-start_datet).total_seconds()*1000000)

	return (start_micros,end_micros)


##
# Public function that calculates best solution for workload passed as a parameter
# Params: benchmark set [Apps] or ([Apps,Clusters], cost function, maximize or not, slots available, max bandwidth, multiprocessing or sequential mode and activate bound or not.
# Return: Best cost and best solution for the workload.
def get_optimal_schedule(benchmark_set, cost_function, maximize, nr_ways, max_bandwidth=float('Inf'),**kwargs):
	## Set defaults
	kwargs.setdefault("bound",True)
	kwargs.setdefault("hint_ucp_slowdown",True)
	kwargs.setdefault("initial_bound",None)
	kwargs.setdefault("multiprocessing",False)
	kwargs.setdefault("user_options",None)
	kwargs.setdefault("low_threshold",100)
	kwargs.setdefault("high_threshold",125)
	kwargs.setdefault("fusion",False)
	kwargs.setdefault("print_times",False)
	kwargs.setdefault("async",False)
	kwargs.setdefault("prepruning",False)
	kwargs.setdefault("postpruning",False)
	kwargs.setdefault("chunk",1)
	kwargs.setdefault("parallel_debug",False)
	kwargs.setdefault("bw_model","simple")
	kwargs.setdefault("use_remote_params",True)
	kwargs.setdefault("paraver",False)

	metadatum = {}

	if "user_options" in kwargs:
		## Merge dict with user options... [Overriding defaults...]
		kwargs.update(kwargs["user_options"])

	## To make things simpler
	bound=kwargs["bound"]
	print_times=kwargs["print_times"]
	async=kwargs["async"]
	chunk=kwargs["chunk"]
	prepruning=kwargs["prepruning"]
	postpruning=kwargs["postpruning"]
	parallel_debug=kwargs["parallel_debug"]
	bw_model=kwargs["bw_model"]
	use_remote_params=kwargs["use_remote_params"]
	paraver =  kwargs["paraver"]

	## Disable options that do not make sense
	if parallel_debug:
		async=False

	if not bound:
		prepruning=False

	if not async:
		paraver=False

	if not prepruning:
		postpruning=False

	log_enabled=paraver or print_times

	## Algorithm main code 
	clustering=not(type(benchmark_set) is list)

	# Determine right operator
	op = operator.gt if maximize else operator.lt
	
	## Be careful, workload can be a clustering now
	if clustering:
		workload=benchmark_set[0]
		cluster_list=benchmark_set[1]
	else:
		workload=benchmark_set
		cluster_list=benchmark_set

	## Obtain a candidate solution to begin with
	if kwargs["initial_bound"]:
		default_cost=kwargs["initial_bound"]
	else:
		if kwargs["hint_ucp_slowdown"] and not clustering: ## UCP SLOWDOWN NOT SUPPORTED FOR NOW FOR CLUSTERING
			heuristic=get_schedule_UCP_gen(workload,nr_ways,metric="slowdown")
		else:
			heuristic=get_equal_llc_schedule(len(cluster_list), nr_ways)

		default_cost=cost_function(benchmark_set,heuristic, max_bandwidth)

	if kwargs["multiprocessing"]:
		
		## Important!!
		best_solution = None
		## Reduction ...
		best_cost = default_cost
		total_branches = 0

		## Start connection right away
		rc = ipp.Client(timeout=30000) #debug=True)
		nr_engines=len(rc.ids)


		dview = rc[:]
		dict_globals={"benchmark_set":benchmark_set,
				"cost_function": cost_function,
				"max_bandwidth": max_bandwidth,
				"op": op,
				"nr_ways": nr_ways,
				"bound": bound,
				"bw_model": bw_model}
		# Copy global data
		ret=dview.apply_sync(lambda x: set_global_properties(x),dict_globals)

		if parallel_debug:
			set_global_properties(dict_globals)


		lview = rc.load_balanced_view()
		lview.block = not async		

		# DEBUG code to check if variables are being set all right
		# ret = lview.map(lambda x: get_global_property(x),[3]*nr_engines, block=True)
		# print ret
		# rc.close()
		# exit(0)
		
		start = time.time()

		if log_enabled:
			if print_times and paraver:
				print >> sys.stderr, "Activated paraver trace generation"
			trace = []
			unlocalized_start = pytz.utc.localize(datetime.datetime.utcnow())
			start_datet = unlocalized_start.astimezone(pytz.timezone("utc"))
		
		nsols=number_of_solutions_partitioning_dp(nr_ways,len(workload))		
		## Unroll one level of the tree only for small scenarios
		if len(workload) <= 2:
			clusters=[]
			initial_solution=[]
			for i in range(1, nr_ways - (len(workload) - 1) + 1):
				initial_solution.append([i])
				clusters.append([[i]])
		else:
			initial_solution=generate_solutions_to_explore(nsols,nr_ways,len(workload),kwargs["high_threshold"],[])			
			
			# Unfold list (generator)
			subsols=[x for x in initial_solution]
			nr_initial_solutions=len(subsols)

			if prepruning:
				if prepruning==2:
					## Update total branches as we are considering them for pruning
					total_branches+=len(subsols)
					chunkval=1

					## Determine which are good
					if use_remote_params:
						if parallel_debug:
							promising= map(lambda x: subsol_is_promising_parallel(best_cost,x), subsols)
						else:	
							promising= lview.map(lambda x: subsol_is_promising_parallel(best_cost,x), subsols,chunksize=chunkval)
					else:
						promising= lview.map(lambda x: subsol_is_promising(benchmark_set, cost_function, max_bandwidth, best_cost, op, nr_ways , x), subsols,chunksize=chunkval)		

					if async:
						promising.get()

						if paraver:
							metadata=promising.metadata

							for i,stats in enumerate(metadata): 
								start_micros=int((stats["started"]-start_datet).total_seconds()*1000000)
								end_micros=int((stats["completed"]-start_datet).total_seconds()*1000000)
								if i==0:
									trace.append(["1","1","1","1","1","0",str(start_micros),"7"])
								trace.append(["1","%i"%(stats['engine_id']+1),"1","1","%i"%(stats['engine_id']+1),str(start_micros),str(end_micros),"6"])
					##Filter with array of bool
					initial_solution= [ (subsol,size)   for i,(subsol,size)   in enumerate(subsols)  if promising[i] ]		
				else:
					initial_solution= [ (subsol,size)   for (subsol,size)  in subsols  if subsol_is_promising(benchmark_set, cost_function, max_bandwidth, best_cost, op, nr_ways ,subsol)  ]
			else:
				#Overwrite iterator
				initial_solution=subsols
				if paraver:
					end_micros=0


			## Fusion tasks if kwargs["fusion"]
			(clusters,unmerged)=fusion_tasks(initial_solution,kwargs["low_threshold"],kwargs["high_threshold"],not kwargs["fusion"])
		
		end = time.time()
		sim_time=str(round(end - start,4))

		if print_times:
			if kwargs["fusion"]:
				str_merge=" (Merged from %d into %d clusters)" % (len(unmerged),len(clusters)) 
			else:
				str_merge=" "

			print >> sys.stderr, "Unroll/Preprune/Merge time: %s (Pruning rate: %s out of %d)%s" % (sim_time,len(clusters),nr_initial_solutions, str_merge)

		if postpruning:
			if log_enabled:
				start = datetime.datetime.utcnow()

			old_clusters=len(clusters)
			clusters=unroll_solutions(nsols,nr_ways,len(workload),kwargs["low_threshold"],clusters)


			if log_enabled:
				end = datetime.datetime.utcnow()
				if print_times:
					print >> sys.stderr, "Postpruning split: %s (Goes from %d to %d solutions)" % (str(round((end - start).total_seconds(),4)),old_clusters,len(clusters))


			if paraver:
				(start_micros,end_micros)=get_start_end_micros(start,end,start_datet)
				trace.append(["1","1","1","1","1",str(start_micros),str(end_micros),"7"])
				# Write header and trace to file


		### PURE PARALLEL STAGE
		nr_subtrees_remaining=len(clusters)
		if chunk==0:
			chunk=nr_subtrees_remaining ## process all
		max_subtrees=nr_engines*chunk
		max_index=len(clusters)-1
		low_idx=0

		while nr_subtrees_remaining> 0:
			if log_enabled:
				start = time.time()

			## Determine which range should be processed
			min_idx=low_idx
			max_idx=low_idx+max_subtrees

			if max_idx>max_index:
				max_idx=max_index+1

			if parallel_debug:
				local_node_solutions =[]
				for x in clusters[min_idx:max_idx]:
					sol=get_optimal_schedule_aux_multi(benchmark_set, cost_function, max_bandwidth, best_cost, op, nr_ways , x, bound)
					local_node_solutions.append(sol)
					print sol
			else:
			#https://ipyparallel.readthedocs.io/en/latest/api/ipyparallel.html#ipyparallel.LoadBalancedView
				if use_remote_params:
					local_node_solutions = lview.map(lambda x: get_optimal_schedule_aux_multi_parallel(best_cost,x,(min_idx!=0) or (not prepruning) or (postpruning)), clusters[min_idx:max_idx])	
				else:
					local_node_solutions = lview.map(lambda x: get_optimal_schedule_aux_multi(benchmark_set, cost_function, max_bandwidth, best_cost, op, nr_ways , x, bound), clusters[min_idx:max_idx])		
			if async: 
				## Now local_node_solutions is an AsyncMapResult object. But it can be iterated.
				## Wait for completion.... invoking get (later we'll use metadata!)
				d=local_node_solutions.get()

			## Debugging code
			if log_enabled:
				end = time.time()
				if print_times:
					print >> sys.stderr, "Parallel computation: %s" % str(round(end - start,4))

				if async:
					## https://ipyparallel.readthedocs.io/en/latest/details.html#extended-interface
					##Extract metadata
					metadata=local_node_solutions.metadata

					for i,stats in enumerate(metadata):
						computation_time=(stats['completed']-stats['started']).total_seconds()
						#transfer_time=(stats['started']-stats['submitted']).total_seconds()+(stats['received']-stats['completed']).total_seconds()
						total_time=(stats['received']-stats['submitted']).total_seconds()
						if print_times:
							print  >> sys.stderr, clusters[low_idx+i],"%i,%i,%.4f,%.4f" % (i,stats['engine_id'],computation_time,total_time),local_node_solutions[i],stats["stdout"]
						if paraver:							
							start_micros=int((stats["started"]-start_datet).total_seconds()*1000000)
							if i == 0:
								trace.append(["1","1","1","1","1",str(end_micros),str(start_micros),"7"])
							end_micros=int((stats["completed"]-start_datet).total_seconds()*1000000)
							# 1(=state record),cpu,app,task,thread,begin,end,state
							trace.append(["1","%i"%(stats['engine_id']+1),"1","1","%i"%(stats['engine_id']+1),str(start_micros),str(end_micros),"1"])


			if log_enabled:
				start = datetime.datetime.utcnow()

			## Sequential reduction
			for local_node_solution in local_node_solutions: 
				(local_best_cost, local_best_solution, local_total_branches)=local_node_solution

				total_branches += local_total_branches
		
				if local_best_solution and op(local_best_cost, best_cost):
					best_cost = local_best_cost
					best_solution = local_best_solution


			if log_enabled:
				end = datetime.datetime.utcnow()
				if print_times:
					print >> sys.stderr, "Sequential reduction: %s" % str(round((end - start).total_seconds(),4))

				if paraver:
					(start_micros,end_micros)=get_start_end_micros(start,end,start_datet)
					trace.append(["1","1","1","1","1",str(start_micros),str(end_micros),"7"])
					# Write header and trace to file
					header_str="#Paraver ({}/{}/{} at {}:{}):{}:1({}):1:1({}:1)\n".format(start_datet.day,start_datet.month,start_datet.year,start_datet.hour,start_datet.minute,end_micros,nr_engines,nr_engines)
					metadatum["header"]=header_str
					metadatum["trace"]=trace

			## Update low_idx for next iteration
			nr_subtrees_remaining=nr_subtrees_remaining-(max_idx-low_idx+1)
			low_idx=max_idx+1

			## END WHILE

		rc.close()
	else:
		## Sequential version
		(best_cost, best_solution, total_branches) = get_optimal_schedule_aux(benchmark_set, cost_function, max_bandwidth, default_cost, op, nr_ways, [], bound)
	
	## Case when best solution matches that of heuristic algorithm
	if not kwargs["initial_bound"] and not best_solution:
		best_solution=heuristic

	return (best_cost, best_solution, total_branches, metadatum)


##
# Public function that calculates a solution for workload passed as a parameter using the lookahead algorithm.
# Params: workload (array of Apps) and slots available.
# Return: cache partitioning for the workload
def get_schedule_UCP(workload, nr_ways):
	return lookahead_algorithm(workload,nr_ways,False)

##
# Public function that calculates a solution for workload passed as a parameter using the lookahead algorithm.
# Params: workload (array of Apps) and slots available.
# Return: cache partitioning for the workload
def get_schedule_UCP_gen(workload, nr_ways, metric="llcmpki"):
	curves=[]
	nr_ways_equal=int(round(nr_ways/len(workload)))

	for app in workload:
		if metric=="inv-speedup":
			curves.append(1/app.get_speedup_table(nr_ways_equal))
		else:
			curves.append(app.get_metric_table(metric))

	return lookahead_algorithm_gen(curves,nr_ways)


# Public function that calculates a solution for each workload using the Yu-Petrof algorithm.
# Params: workload (array of Apps) and optionally, a threshold.
# Return: cache partitioning for the workload
def get_schedule_yu_petrof(workload, nr_ways, threshold = 0):
	bandwidth_reduction_table = build_bandwidth_reduction_table(workload, nr_ways)
	return yu_petrof_algorithm(workload, bandwidth_reduction_table, nr_ways, threshold)


# Weighted access of n mrcs with n float indexes
def access_indexes(mrcs, indexes):
	v_ns = []
	for i in range(len(indexes)):
		r = indexes[i] % 1
		int_i = int(indexes[i])
		if int_i >= len(mrcs[i]):
			v_ns.append(mrcs[i][int_i-1])
		elif (r == 0 or int_i+1 >= len(mrcs[i])):
			v_ns.append(mrcs[i][int_i])
		else:
			v_ns.append((1-r)*mrcs[i][int_i] + r*mrcs[i][int_i+1])
	return v_ns

# Original function from kpart code, combines curves in pairs
def whirlpool_combined_curve(mrc1, mrc2):
	n = len(mrc1)
	assert(n==len(mrc2))
	index1 = 0.0
	index2 = 0.0
	mrc = []
	buckets = []
	for cur_ways in range(n):
		v1,v2 = access_indexes([mrc1,mrc2], [index1,index2])
		v_sum = v1+v2
		index1 += v1/v_sum
		index2 += v2/v_sum
		mrc.append(sum(access_indexes([mrc1,mrc2], [index1,index2])))
		b1 = int(round(index1))
		buckets.append((b1,cur_ways+1-b1))
	return mrc,buckets

def whirlpool_combine_ncurves_f(mrcs,nr_ways):
	nr_apps = len(mrcs)
	indexes = [0.0 for i in range(nr_apps)]
	mrc = []
	buckets = []
	for cur_ways in range(nr_ways):
		v_ns = access_indexes(mrcs,indexes)
		v_sum = sum(v_ns)
		for i in range(nr_apps):
			indexes[i] += v_ns[i]/v_sum
		mrc.append(sum(access_indexes(mrcs,indexes)))
		# Save index to return the buckets in the same order as the mrcs
		total_ways = cur_ways +1
		cur_buckets = []
		total_weight = sum(indexes)
		for i in range(nr_apps):
			app_bucket=total_ways*indexes[i]/total_weight
			cur_buckets.append(app_bucket)
		buckets.append(cur_buckets)
	return (mrc,buckets)


def whirlpool_combine_ncurves(mrcs,nr_ways):
	nr_apps = len(mrcs)
	indexes = [0.0 for i in range(nr_apps)]
	mrc = []
	buckets = []
	for cur_ways in range(nr_ways):
		v_ns = access_indexes(mrcs,indexes)
		v_sum = sum(v_ns)
		for i in range(nr_apps):
			indexes[i] += v_ns[i]/v_sum
		mrc.append(sum(access_indexes(mrcs,indexes)))
		# Save index to return the buckets in the same order as the mrcs
		sorted_ixs = [(i,indexes[i]) for i in range(nr_apps)]
		sorted_ixs = sorted(sorted_ixs,key=lambda x:x[1],reverse=True)
		total_ways = cur_ways +1
		# Assign ways at least 1 way per app
		cur_buckets = [0 for i in range(nr_apps)]
		for i in range(nr_apps):
			if total_ways > 0:
				cur_buckets[sorted_ixs[i][0]] = 1
				total_ways -= 1
		if(total_ways > 0):
			total_weight = sum(indexes)
			for i in range(nr_apps):
				app_bucket = int(round((total_ways*sorted_ixs[i][1])/total_weight))
				cur_buckets[sorted_ixs[i][0]]+=int(app_bucket)
				total_ways -= app_bucket
				total_weight -= sorted_ixs[i][1]

		buckets.append(cur_buckets)
	return (mrc,buckets)

def get_schedule_whirlpool(workload, nr_ways, metric="llcmpki"):
	curves=[]
	for app in workload:
		curves.append(app.get_metric_table(metric).values)

	(mrc,buckets)=whirlpool_combine_ncurves(curves,nr_ways)
	return buckets[nr_ways-1]	

def get_schedule_whirlpool_float(workload, nr_ways, metric="llcmpki"):
	patched_workload=get_scaled_properties_cluster(workload, nr_ways, metric=metric)
	# The table entry to access is max_ways
	per_app_ways = [nr_ways for app in workload]
	# They are all in the same cluster
	full_mask = get_partition_masks([nr_ways])
	per_app_masks = [full_mask[0] for app in workload]
	cluster_id = [1 for app in workload]
	return (patched_workload,per_app_ways,per_app_masks,cluster_id)


def start_ipcluster():
	call("cd /home/bench/Jorge/TFM_JorgeC/simulator && . ./shrc && ipcluster start -n 20 &", shell=True)
	time.sleep(10)

def stop_ipcluster():
	call("pkill ipcluster --signal SIGINT", shell=True)
	time.sleep(30)

# Input, list of way assignments...
def minimize_overlapping_partitioning(partitions,nr_ways,fix_intel_bug=True):
	left=True
	nr_partitions=len(partitions)
	overlap_counter=[]
	masks=[]
	clusters=[]

	for i in range(nr_ways):
		overlap_counter.append(0)

	for i in range(nr_partitions):
		nr_ways_to_assign=partitions[i]
		##find the way with the lowest degree of overlapping
		if left:
			j=1
			cur_min=overlap_counter[0]
			min_idx=0
			
			while j<nr_ways and cur_min>0:
				if cur_min>overlap_counter[j]:
					cur_min=overlap_counter[j]
					min_idx=j
				j=j+1
			## Fill to the left
			start=min_idx

		else:
			j=nr_ways-2
			cur_min=overlap_counter[nr_ways-1]
			min_idx=nr_ways-1		

			while j>=0 and cur_min>0:
				if cur_min>overlap_counter[j]:
					cur_min=overlap_counter[j]
					min_idx=j
				j=j-1
			## Fill to the left
			start=min_idx-nr_ways_to_assign+1
	
		#Intel fix (Not alone...)	
		if fix_intel_bug and nr_ways==11 and start==10 and nr_ways_to_assign==1:
			start-=1 # Move it one to the left..


		clusters.append((start,nr_ways_to_assign))
		masks.append(hex((1<<nr_ways_to_assign)-1<<start))

		## Increment overlapping counters
		for k in range(start,start+nr_ways_to_assign):
			overlap_counter[k]+=1

		## REverse order for next iteration
		left=(not left)

	return (masks,clusters)

def get_dunn_index(n_apps,k,normalized_stalls,partitions):
	# Precompute [min,max] of each cluster
	clusters = [[1,0] for i in range(k)]
	for i in range(n_apps):
		if normalized_stalls[i] < clusters[partitions[i]][0]:
			clusters[partitions[i]][0] = normalized_stalls[i]
		if normalized_stalls[i] > clusters[partitions[i]][1]:
			clusters[partitions[i]][1] = normalized_stalls[i]

	# Order to have adjacent clusters
	clusters = sorted(clusters,key=lambda x:x[0])

	# In the selfa17 paper: Dunn(k) = d_min/d_max
	d_min = 1 # min distance between points of diff clusters
	d_max = 0 # max within-cluster distance
	for i in range(k):
		intra_cluster_distance = clusters[i][1] - clusters[i][0]
		if (intra_cluster_distance > d_max):
			d_max = intra_cluster_distance
		if (i<k-1):
			inter_cluster_distance = clusters[i+1][0] - clusters[i][1]
			if inter_cluster_distance < d_min:
				d_min = inter_cluster_distance
	return d_min / d_max


def estimate_ways_exponential(x, max_ways):
	if max_ways == 20:
		nways = int(round(19.26833327*np.exp(0.67431637*x) - 17.72167744))
	elif max_ways == 11:
		# parameters for 10 ways scaling down the paper ways
		nways = int(round(4.09659516*np.exp(1.19540837*x) - 2.48157709))
	return nways

def get_schedule_dunn(workload, nr_ways):
	n_apps = len(workload)
	curves=[]
	for app in workload:
		curves.append(app.get_metric_table("llcmpki").values)

	(mrc,buckets)=whirlpool_combine_ncurves(curves,nr_ways)

	app_ways=buckets[nr_ways-1]
	stalls = []
	for i,app in enumerate(workload):
		stalls.append(app.properties["stalls_l2_miss"][app_ways[i]])

	stalls_min = min(stalls)
	stalls_max = max(stalls)
	normalized_stalls = []
	for stalls_app in stalls:
		if stalls_max != stalls_min:
			n_stalls = (stalls_app-stalls_min)/(stalls_max-stalls_min)
		else:
			n_stalls = 1/n_apps
		normalized_stalls.append(n_stalls)


	dunn_min = 100000
	best_k = 2
	best_centroids = []
	best_partitions = []
	for k in range(2,n_apps):
		kmeans = KMeans(n_clusters=k)
		kmeans = kmeans.fit(np.reshape(normalized_stalls,(-1,1)))
		# Array with the centroid of each partition (odered each time different)
		centroids = kmeans.cluster_centers_
		# Position of the centroid array to which each app belongs
		partitions = kmeans.labels_
		dunn=get_dunn_index(n_apps,k,normalized_stalls,partitions)
		if dunn < dunn_min:
			best_centroids = centroids
			best_partitions = partitions
			best_k = k
			dunn_min = dunn
	# Each iteration best_partition can reference to different centroids(ordered != each time)
	o_centroids = [(old_i,c) for old_i,c in enumerate(best_centroids)]
	o_centroids = sorted(o_centroids,key=lambda x:x[1],reverse=True)
	cluster_new_i = {old_i:ordered_i for ordered_i,(old_i,c) in enumerate(o_centroids)}
	# reassign apps cluster index in partitions so it references to the reordered cluster
	for i in range(len(best_partitions)):
		best_partitions[i] = cluster_new_i[best_partitions[i]]
	centroids_ways = [estimate_ways_exponential(x[1],nr_ways) for x in o_centroids]
	
	# Check if there is left ways and assign them to the smallest clusters
	left_ways = nr_ways - sum(centroids_ways)
	i_smallest = len(centroids_ways)-1
	while (left_ways > 0):
		centroids_ways[i_smallest] += 1
		left_ways -= 1
		i_smallest = (i_smallest - 1) % len(centroids_ways)

	masks, part_info =  minimize_overlapping_partitioning(centroids_ways,nr_ways)

	per_app_masks = []
	per_app_ways = []
	cluster_ids = []
	per_app_centroids = []
	for i_app in range(len(best_partitions)):
		# Assign CLOS to PID
		per_app_masks.append(masks[best_partitions[i_app]])
		per_app_ways.append(centroids_ways[best_partitions[i_app]])
		per_app_centroids.append("%0.3f"%o_centroids[best_partitions[i_app]][1][0])
		cluster_ids.append(i_app)

	return (per_app_ways,per_app_masks,per_app_centroids,cluster_ids)

def merge_clusters(cluster1,cluster2,curve,buckets):
	apps=cluster1.apps + cluster2.apps
	scaled_apps = []
	for i,app in enumerate(apps):
		assigned_ways = [frac_ways[i] for frac_ways in buckets]
		scaled_apps.append(App(app.name,app.build_scaled_properties(assigned_ways),parent=app))
	return Cluster(scaled_apps,curve,buckets)

class Cluster:
	def __init__(self,apps,curve,buckets):
		self.apps = apps
		self.curve = curve
		self.buckets = buckets

	def __repr__(self):
		names=[]
		for idx,app in enumerate(self.apps):
			names.append(app.name)
		return repr(names)

	def distance(self,cluster2):
		curves=[]
		for app in self.apps:
			curves.append(app.original_app.get_metric_table("llcmpki"))
		for app in cluster2.apps:
			curves.append(app.original_app.get_metric_table("llcmpki"))

		combined_curve,buckets = whirlpool_combine_ncurves_f([c.values for c in curves],len(self.curve))

		partitioned_curve=[]
		# Calculate partitioned curve
		for nr_ways in range(1,len(self.curve)+1):
			ways_distr = lookahead_algorithm_gen(curves,nr_ways)
			partitioned_sum = 0
			for idxapp,assigned_ways in enumerate(ways_distr):
				partitioned_sum += curves[idxapp][assigned_ways]
			partitioned_curve.append(partitioned_sum)

		distance_sum = 0
		for i in range(len(self.curve)):
			distance_sum += abs(combined_curve[i] - partitioned_curve[i])

		return (distance_sum,combined_curve,buckets,partitioned_curve)	


	def get_cluster_speedup_curve(self,nr_ways,total_apps,max_bandwidth):
		speedup_curves=[]
		nr_ways_equal=int(round(nr_ways/total_apps))
		# Now each app has stored the proportional metrics values to the float ways
		sp=[]
		for i_ways in range(1,nr_ways+1):
			aggregate_val=0
			for app in self.apps:
				aggregate_val+=app.get_speedup_table(nr_ways_equal)[i_ways]
			sp.append(aggregate_val)
		return pd.Series(sp,index=range(1,nr_ways+1))


def determine_aggregate_speedup(partitioning,nr_ways,max_bandwidth):
	bw_exceeded=False
	nr_ways_equal=int(round(nr_ways/len(partitioning)))

	## Calculate total bandwidth for that partitioning
	bw_aggregate=0
	sp_aggregate=0
	bw_app=[]

	for app,ways in partitioning:
		bw_cur=app.get_metric_table("bandwidth_mbps")[ways]
		bw_aggregate+=bw_cur

	if max_bandwidth != float('Inf') and bw_aggregate>max_bandwidth:
		bw_exceeded=True
		scale_factor=bw_aggregate/max_bandwidth
		for bw,app in bw_app:
			app.scaling_factor=1.0-((bw/bw_aggregate)*scale_factor)	

	for app,ways in partitioning:
		sp_aggregate+=app.get_speedup_table(nr_ways_equal)[ways]

	if bw_exceeded:
		for bw,app in bw_app:
			app.scaling_factor=1.0

	return sp_aggregate


def determine_best_partitioning(clusters,nr_ways,total_apps,max_bandwidth=float('Inf')):
	speedup_curves=map(lambda x: x.get_cluster_speedup_curve(nr_ways,total_apps,max_bandwidth), clusters)
	inverse_sp_curves=[1/curve for curve in speedup_curves]
	nr_ways_equal=int(round(nr_ways/total_apps))

	# Lookahead
	if len(clusters)==1:
		partitioning=[nr_ways]
	else:
		partitioning=lookahead_algorithm_gen(inverse_sp_curves,nr_ways)

	## Determine Per-application way-assignment
	per_app_partitioning=[]
	bw_aggregate=0
	sp_aggregate=0
	bw_app=[]

	for idxc,nr_ways_cluster in enumerate(partitioning):
		## Traverse applications in a cluster to determine 
		cl=clusters[idxc]
		buckets=cl.buckets
		for idxa,app in enumerate(cl.apps):
			# Now stores the cluster ways so it can later retrieve the scaled value from properties
			per_app_partitioning.append((app,nr_ways_cluster,idxc))
			sp_aggregate+=app.get_speedup_table(nr_ways_equal)[nr_ways_cluster]
			bw_cur=app.get_metric_table("bandwidth_mbps")[nr_ways_cluster]
			bw_aggregate+=bw_cur
			bw_app.append((bw_cur,app))

	## See if we exceed the Total system BW
	if max_bandwidth != float('Inf') and bw_aggregate>max_bandwidth:
		scale_factor=bw_aggregate/max_bandwidth
		## Per-app scaling factor
		for bw,app in bw_app:
			app.scaling_factor=1.0-((bw/bw_aggregate)*scale_factor)

		(partitioning,per_app_partitioning,bw_app,sp_aggregate)=determine_best_partitioning(clusters,nr_ways,total_apps,float('Inf'))

		## Restore scaling factors .... 
		for bw,app in bw_app:
			app.scaling_factor=1.0
	return (partitioning,per_app_partitioning,bw_app,sp_aggregate)



def get_kpart_schedule(workload,nr_ways,max_bandwidth=float('Inf')):
	curClusters = []
	total_apps = len(workload)	
	for i,app in enumerate(workload):
		buckets = [[nw] for nw in range(1,nr_ways + 1)]
		curve = app.get_metric_table("llcmpki").values
		app.bench_id = i
		cl=Cluster([app],curve,buckets)
		curClusters.append(cl)
	
	solutions = [(list(curClusters),determine_best_partitioning(curClusters,nr_ways,total_apps,max_bandwidth))]

	while len(curClusters) > 1:

		cluster_data=None
		min_found=30000000
		min_idx=(-1,-1)
		
		## Compute distance and keep track of min
		for idx,clusterRef in enumerate(curClusters):
			for i in range(idx+1,len(curClusters)):
				cluster=curClusters[i]
				(distance,combined_curve,buckets,partitioned_curve)=clusterRef.distance(cluster)

				if distance < min_found:
					min_found=distance
					min_idx=(idx,i)
					cluster_data=(distance,combined_curve,buckets,partitioned_curve)

		## Merge 2 closest clusters
		new_cluster=merge_clusters(curClusters[min_idx[0]],curClusters[min_idx[1]],cluster_data[1],cluster_data[2])

		del curClusters[min_idx[0]]
		del curClusters[min_idx[1]-1]

		#Add merged cluster ...
		curClusters.insert(0,new_cluster)
		solutions.append((list(curClusters),determine_best_partitioning(curClusters,nr_ways,total_apps,max_bandwidth)))	

	## Determine BEST K
	max_sp=-1
	best_k=0

	for (idx,sol) in enumerate(solutions):
		(clusters,data)=sol
		(partitioning,per_app_partitioning,bw_app,sp_aggregate)=data

		if sp_aggregate > max_sp:
			max_sp=sp_aggregate
			best_k=idx

	return (solutions,best_k)


def get_kpart_best_gen(workload,nr_ways,max_bandwidth=float('Inf'),debugging=False):
	## Invoke kpart.....
	(solutions,k)=get_kpart_schedule(workload,nr_ways,max_bandwidth)

	if debugging:
		for sol in solutions:
			clusters,(partitioning,per_app_partitioning,bw_app,sp_aggregate) = sol
			bw_app = [(bw,app.name) for (bw,app) in bw_app]
			(workload_alt,solution_alt)=zip(*per_app_partitioning)
			print sp_aggregate, unfairness(list(workload_alt),list(solution_alt),max_bandwidth), clusters, partitioning, list(solution_alt)
		print "Best k =",k

	best_sol=solutions[k]
	clusters=best_sol[0]
	per_cluster_part=best_sol[1][0]
	per_cluster_masks=get_partition_masks(per_cluster_part)
	per_app_partitioning=best_sol[1][1]


	per_app_masks=[None]*len(per_app_partitioning)
	per_app_ways=[None]*len(per_app_partitioning)
	cluster_ids=[None]*len(per_app_partitioning)
	clusters = [[] for _ in range(len(per_cluster_masks))]
	## Traverse apps in the same order as the per_app_partitioning vector
	for (app,way_clust,idxc) in per_app_partitioning:
		orig_i = app.original_app.bench_id
		per_app_ways[orig_i] = way_clust
		per_app_masks[orig_i] = per_cluster_masks[idxc]
		cluster_ids[orig_i] = idxc
		clusters[idxc].append(app.original_app)

	patched_workload=[None]*len(per_app_partitioning)
	# Return patched workload with apps scaled with llcmpkc whirlpool 
	for cluster in clusters:
		scaled_apps = get_scaled_properties_cluster(cluster,nr_ways)
		for app in scaled_apps:
			patched_workload[app.original_app.bench_id] = app

	return (patched_workload,per_app_ways,per_app_masks,cluster_ids)

def get_kpart_best(workload,nr_ways,max_bandwidth=float('Inf'),debugging=False):
	return get_kpart_best(workload,nr_ways,max_bandwidth,debugging)[1]

			
## Function that builds a hacked set of apps (whirlpool properties) that belong to a cluster
## Cluster is just a plain list of apps
## Return list of "Shadow" apps with the hacked properties...
def get_scaled_properties_cluster(cluster, max_ways, metric="llcmpkc"):
	## Caso trivial...
	if len(cluster)==1:
		return cluster
	## Get MRC for each APP
	curves=map(lambda x: x.get_metric_table(metric).values,cluster)
	## Apply whirlpool
	comb_mrc,buckets = whirlpool_combine_ncurves_f(curves,max_ways)
	hacked_apps=[]
	for i,app in enumerate(cluster):
		props=app.build_scaled_properties([b_ways[i] for b_ways in buckets])
		hacked_apps.append(App(app.name,props,app))


	return hacked_apps

if __name__ == "__main__":
#	l=[c for c in generate_possible_clusters(["A","B","C","D","E"],3)]
#	print l
#	print len(l)
#	exit(0)
	write_dot_tree(6,4)
	#print determine_number_of_cluster_cache_partitioning(20,20)	
	#nsols=number_of_solutions_partitioning_dp(11,4)		
	#solutions=[[[1]]]
	#sols=unroll_solutions(nsols,11,4,10000,solutions)	
	#print sols
	#sol=number_of_solutions_partitioning_dp(11,11)
	#print "***All solutions:***"
	#print sol
	#print "***Solutions recursive:***"
	#print map(lambda x: number_of_solutions_partitioning(20,x),range(1,10))
	#print "***Solutions Dynamic Programming:***"
	#print sol[10]	
	#print "***Solutions Clustering 11 11 y 20 20:***"
	#print map(lambda x: determine_number_of_cluster_cache_partitioning(11,x),range(1,13))
	#print map(lambda x: determine_number_of_cluster_cache_partitioning(20,x),range(1,13))