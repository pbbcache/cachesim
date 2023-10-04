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
from __future__ import print_function
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
#from sklearn.cluster import KMeans
import time
import functools
import datetime
import pytz
import os
from bisect import bisect_left
import threading
import zmq
import socket
import uuid
from ipyparallel.util import disambiguate_url
import ctypes 
from functools import reduce

# Function to parse a harness file
def parse_harness_file(filename):
	# Best 3, 6 (split / and get last)
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
	


glb_algorithms_supported=["opt-stp","yu-petrov","opt-unf","ucp","ucp-slowdown","kpart","whirlpool","whirlpool-c","equal-part","optc-stp","optc-unf","opt-unf-np","opt-stp-np","yu-petrov","user","kpart-opt","opt-stp-bf","opt-unf-bf","opt-bf","optc-stp-bf","optc-unf-bf","optc","opt-map-stp","opt-map-unf","opt-map-antt","opt-map-fair","opt-map-cov","opt-map"]

# Function to add extra algorithms to supported global list
def add_extra_algorithms(*algs):
	global glb_algorithms_supported
	for alg in algs:
		glb_algorithms_supported.append(alg)

# Function to obtain the global supported algorithms
def get_algorithm_names():
	global glb_algorithms_supported
	return glb_algorithms_supported


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
	#bw_alone_vector=list(bw_alone_vector)
	equations=[]
	U=sp.Symbol('U')
	variables=[U]
	guess=[sum(bw_alone_vector)]
	## To build equations dynamically
	eq=-U

	#old loop
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

## Topology is another global simulator parameter
glb_active_topology=0
glb_available_topologies=["uma","numa"]

# Function to select which topology to use
def select_topology(str_topo):
	global glb_active_topology
	if str_topo in glb_available_topologies:
		glb_active_topology=glb_available_topologies.index(str_topo)
	else:
		print("No such topology: ", str_topo, file=sys.stderr)
		exit(1)

## Global simulator parameter
glb_active_bw_model=2
glb_available_bw_models=["morad","stalls","simple"]

# Function to select which bandwidth model to use (traditional, improved, etc.)
def select_bw_model(str_model):
	global glb_active_bw_model
	if str_model in glb_available_bw_models:
		glb_active_bw_model=glb_available_bw_models.index(str_model)
	else:
		print("No such BW model: ", str_model, file=sys.stderr)
		exit(1)

# Function that applies the bandwidth model for a specific application list, the assigned ways and the maximum observed bandwidth
def apply_bw_model_gen(benchmark_list,ways_assigned,max_bandwidth):
	ways_assigned = list(map(int, ways_assigned))
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


	if glb_active_bw_model<=1:
		bw_alone_vnorm=list(map(lambda x: x/max_bandwidth,bw_alone_vector))

		## Invoke solver to determine approximate solution
		solver_solution=estimate_bw_shared(bw_alone_vnorm)
	
		## Get max predicted BW from solution
		predicted_total_bw=solver_solution[0]

		# Get BW shared for each app
		predicted_bw_shared=solver_solution[1:]
	else:
		aggregate_alone_bw=sum(bw_alone_vector)
		if aggregate_alone_bw<=max_bandwidth:
			## No saturation, no problem 
			predicted_bw_shared=bw_alone_vector
		else:
			## KPART SIMPLE MODEL...
			max_b_alone=max(bw_alone_vector)
			sum_square=(1.0/max_b_alone)*functools.reduce(operator.add,map(lambda bw: bw**2,bw_alone_vector))
			xm=(max_bandwidth+sum_square-aggregate_alone_bw)/sum_square
			predicted_bw_shared=[bi*(1-(bi/max_b_alone)*(1-xm)) for bi in bw_alone_vector]

		## In this case, predicted_bw_shared and bw_alone_vnorm are not really normalized
		## But maintaining a uniform format is better
		bw_alone_vnorm=bw_alone_vector


	overall_slowdown=[]
	
	## Determine actual slowdown by considering BW-related slowdown 
	for i,bw_alone in enumerate(bw_alone_vnorm): 
		## Alter bw_slowdown_factor according to stall-based model
		if glb_active_bw_model==1:
			s_alone=total_stalls_alone_v[i]
			m_alone=mem_access_stalls_alone_v[i]
			## bw_slowdown_factor is Old value provided by Morad's model
			## Determine new bw_slowdown factor
			old_morad=float(bw_alone/predicted_bw_shared[i])
			bw_slowdown_factor=(s_alone+m_alone*(old_morad-1))/s_alone
		else:
			## Simple models (Slowdown is proportional to BW reduction)
			bw_slowdown_factor=bw_alone/predicted_bw_shared[i]			

		## Calculate performance reduction 
		overall_slowdown.append(slowdown_vector[i]*bw_slowdown_factor)

	return overall_slowdown

# Function to estimate the bandwidth shared of each application form an array of bandwidth alone and the maximum observed bandwidth
def determine_bw_shared(bw_alone_vector,max_bandwidth):
	bw_alone_vnorm=list(map(lambda x: x/max_bandwidth,bw_alone_vector))

	## Invoke solver to determine approximate solution
	solver_solution=estimate_bw_shared(bw_alone_vnorm)

	# Return Agg BW and BW shared for each app
	return (solver_solution[0]*max_bandwidth,map(lambda x: x* max_bandwidth,solver_solution[1:]))		



def apply_bw_model_topology(benchmark_list,ways_assigned,max_bandwidth):
	if glb_active_topology==0: ## UMA
		return apply_bw_model_gen(benchmark_list,ways_assigned,max_bandwidth)
	elif glb_active_topology==1: ## NUMA
		## Traverse applications to determine which one goes with
		## and treat each thing separately
		llc_id_count=0
		app_llc_id=[]
		dict_groups={}
		for i,app in enumerate(benchmark_list):
			llc_id=app.llc_id
			if llc_id==-1:
				llc_id=0 ## Normalize
			if llc_id in dict_groups:
				(list_apps,idx_apps,ways_apps)=dict_groups[llc_id]
			else:
				list_apps=[]
				idx_apps=[]
				ways_apps=[]
				dict_groups[llc_id]=(list_apps,idx_apps,ways_apps)
				## New llc ID detected
				llc_id_count=llc_id_count+1 
			## Keep track of normalized LLC_ID
			app_llc_id.append(llc_id)
			## Store per LLCID info
			list_apps.append(app)
			idx_apps.append(i)
			ways_apps.append(ways_assigned[i])

		if llc_id_count==1: ## Single group?
			return apply_bw_model_gen(benchmark_list,ways_assigned,max_bandwidth)
		else:

			## Populate and combine slowdowns
			slowdowns=[1 for i in range(len(benchmark_list))]

			for llc_id in range(llc_id_count):
				(list_apps,idx_apps,ways_apps)=dict_groups[llc_id]
				## Calculate BW and cache contention assuming local memory allocation
				sv=apply_bw_model_gen(list_apps,ways_apps,max_bandwidth)
				for j,idx in enumerate(idx_apps):
					slowdowns[idx]=sv[j]

			return slowdowns


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
		self.bench_id=-1
		self.cluster_in_llc=-1
		self.llc_id=-1
		## For kpart
		if parent:
			self.original_app=parent.original_app
			self.bench_id=parent.bench_id
			self.llc_id=parent.llc_id
		else:
			self.original_app=self		
		## Added for efficiency reasons
		self.nr_ways=len(self.original_app.properties["ipc"])


	def __repr__(self):
		return repr((self.name))

	def show(self):
		print("==== %s ====" % self.name) 
		print(self.properties)

	def get_metric_table(self,metric):
		return self.properties.sort_index()[metric]

	def get_speedup_table(self,nr_ways_fair):
		if nr_ways_fair==0:
			nr_ways_fair=1
		ipc_ref=self.original_app.properties["ipc"][nr_ways_fair]
		return 	(self.properties["ipc"]*(self.scaling_factor/ipc_ref)).sort_index()

	def get_ipc_alone(self):
		ipc_table=self.original_app.properties["ipc"]
		return ipc_table[self.nr_ways]
     
	## Scaled ways is a python list of floats
	## (proportional ways in cluster the application belongs to )
	# metrics for Kpart: "bandwidth_mbps","ipc","slowdown","llcmpki"
	def build_scaled_properties(self,scaled_ways,metrics=[]):
		## Work with original app!!!
		#start= time.time()

		if metrics:
			bdf = self.original_app.properties.loc[:, metrics]
		else:
			bdf = self.original_app.properties.drop("BENCH",axis=1)

		# Initiallize dict to create df
		res_dict = {c:[] for c in bdf.columns}
		res_dict["float_ways"] = []
		# Build zero interpolation
		projected_values = bdf.loc[1]+(bdf.loc[1]-bdf.loc[2])
		projected_values = projected_values.apply(lambda x: 0.000001 if x<=0 else x)
		float_ways = []
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
			for c in new_row.index:
				res_dict[c].append(new_row[c])

		float_df = pd.DataFrame(res_dict,index=[i+1 for i in range(len(res_dict["float_ways"]))])

		#end = time.time()
		#sim_time = str(round(end - start, 4))
		#print "Build_props time: %s" % sim_time
		#print "fl_df: %s" % float_df.columns
		#print "fl_df: %s" % scaled_ways


		return float_df

def get_slowdown_vector(benchmark_list, solution, max_bandwidth=float('Inf')):
	solution = list(map(int, solution))
	## Apply BW model if necessary
	if max_bandwidth != float('Inf'):	
		## Replace slowdown vector ...
		return apply_bw_model_topology(benchmark_list[0:len(solution)],solution,max_bandwidth)
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
		order_apps=[]
		for i,cluster in enumerate(cluster_list):
			cluster_ways=final_solution[i]
			for app in cluster:
				plain_solution.append(cluster_ways)
				order_apps.append(app)

		## Overwrite variables...
		final_solution=plain_solution
		apps=order_apps
		
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

	throughput_value=round(throughput_value,3)

	return throughput_value

# Throughput function.
# Params: Benchmarks data, solution dictionary, max bandwidth and remaining ways (to know if it is a partial solution, used for bounding)
# Return: Throughput value calculated.
def stp(benchmark_set, solution, max_bandwidth=float('Inf'), remaining_ways=0):
	return throughput(benchmark_set, solution, max_bandwidth=max_bandwidth, remaining_ways=remaining_ways)

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

# Geometric mean function.
# Params: Benchmarks data, solution dictionary, max bandwidth and remaining ways (to know if it is a partial solution, used for bounding)
# Return: Geometric mean value calculated.
def gmean_speedup(benchmark_set, solution, max_bandwidth=float('Inf'), remaining_ways=0):
	##Keep track of 
	(benchmark_list,final_solution,apps_in_partial_solution)=normalize_solution(benchmark_set, solution,remaining_ways)
	slowdown_vector=get_slowdown_vector(benchmark_list,final_solution,max_bandwidth)

	## Calculate gmean
	gmean_value=1
	for slowdown in slowdown_vector:
		gmean_value *= 1/slowdown

	gmean_value=gmean_value**(1/len(slowdown_vector))

	del slowdown_vector
	del final_solution

	gmean_value=round(gmean_value,3)

	return gmean_value


# Throughput function.
# Params: Benchmarks data, solution dictionary, max bandwidth and remaining ways (to know if it is a partial solution, used for bounding)
# Return: Throughput value calculated.
def gmean_ipc(benchmark_set, solution, max_bandwidth=float('Inf'), remaining_ways=0):
	##Keep track of 
	(benchmark_list,final_solution,apps_in_partial_solution)=normalize_solution(benchmark_set, solution,remaining_ways)
	slowdown_vector=get_slowdown_vector(benchmark_list,final_solution,max_bandwidth)		
	predicted_ipc=[ app.get_ipc_alone()/slowdown_vector[i] for i,app in enumerate(benchmark_list)]

	## Calculate throughput
	throughput_value=1.0
	for ipc in predicted_ipc:
		throughput_value *= ipc

	throughput_value=round(throughput_value**(1.0/len(predicted_ipc)),3)

	del slowdown_vector
	del final_solution
	del predicted_ipc

	return throughput_value

# Harmonic mean function.
# Params: Benchmarks data, solution dictionary, max bandwidth and remaining ways (to know if it is a partial solution, used for bounding)
# Return: Harmonic mean value calculated.
def hmean_speedup(benchmark_set, solution, max_bandwidth=float('Inf'), remaining_ways=0):
    ##Keep track of 
    (benchmark_list,final_solution,apps_in_partial_solution)=normalize_solution(benchmark_set, solution,remaining_ways)
    slowdown_vector=get_slowdown_vector(benchmark_list,final_solution,max_bandwidth)

    hmean_value=len(slowdown_vector)/sum([s for s in slowdown_vector])  ## slowdown is the inverse of speedup

    hmean_value=round(hmean_value,3)

    del slowdown_vector
    del final_solution

    return hmean_value


# Aggregate IPCs function.
# Params: Benchmarks data, solution dictionary, max bandwidth and remaining ways (to know if it is a partial solution, used for bounding)
# Return: Aggregate IPC value calculated.
def aggregate_ipc(benchmark_set, solution, max_bandwidth=float('Inf'), remaining_ways=0):
	## Keep track of
	(benchmark_list, final_solution, apps_in_partial_solution) = normalize_solution(benchmark_set, solution, remaining_ways)
	slowdown_vector = get_slowdown_vector(benchmark_list, final_solution, max_bandwidth)
	ipc_values=[ app.get_ipc_alone()/slowdown_vector[i] for i,app in enumerate(benchmark_list)]

	## Calculate aggregate IPC
	agg_ipc_value=sum(ipc_values)

	del final_solution
	del ipc_values

	return agg_ipc_value

#  M1 fairness metric function.
#  Params: benchmark data, solution dictionary, max bandwidth and remaining ways (to know if it is a partial solution, used for bounding).
#  Return: M1 metric calculated.
def m1_metric(benchmark_set, solution, max_bandwidth=float('Inf'), remaining_ways=0):
    ## Keep track of
    (benchmark_list, final_solution, apps_in_partial_solution) = normalize_solution(benchmark_set, solution, remaining_ways)
    slowdown_vector = get_slowdown_vector(benchmark_list, final_solution, max_bandwidth)

    ## Calculate M1 metric
    m1_metric = 0
    ## Calculate correct lower bound... (assuming the rest slowdowns can be made as close as possible)
    if remaining_ways>0:
        num_apps = apps_in_partial_solution + 1 ##len(slowdown_vector)
        seq=slowdown_vector[apps_in_partial_solution:]
        ## Patch slowdown average
        slowdown_vector[apps_in_partial_solution]=reduce(lambda x, y: x + y, seq ) / float(len(seq))
    else:
        num_apps = len(slowdown_vector)		
  
    for i in range(num_apps):
        for j in range(i+1, num_apps):
            m1_metric += abs(slowdown_vector[i] - slowdown_vector[j])

    del slowdown_vector
    del final_solution

    m1_metric = round(m1_metric, 3)

    return m1_metric

# ANTT METRIC function.
# Params: Benchmarks data, solution dictionary, max bandwidth and remaining ways (to know if it is a partial solution, used for bounding)
# Return: Unfairness value calculated.
def antt_metric(benchmark_set, solution, max_bandwidth=float('Inf'), remaining_ways=0):
	(benchmark_list,final_solution,apps_in_partial_solution)=normalize_solution(benchmark_set, solution,remaining_ways) 
	slowdown_vector=get_slowdown_vector(benchmark_list,final_solution,max_bandwidth)

	antt_value = sum(slowdown_vector) / len(slowdown_vector)
	del slowdown_vector
	return antt_value

def antt(benchmark_set, solution, max_bandwidth=float('Inf'), remaining_ways=0):
	return antt_metric(benchmark_set, solution, max_bandwidth, remaining_ways)

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

	unfairness_value=round(unfairness_value,3)
	throughput_value=round(throughput_value,3)
	del slowdown_vector
	return (unfairness_value,-throughput_value)


def unf_stp(benchmark_set, solution, max_bandwidth=float('Inf'), remaining_ways=0):
	return unfairness_max_throughput(benchmark_set, solution, max_bandwidth, remaining_ways)


def fairness_metric(benchmark_set, solution, max_bandwidth=float('Inf'), remaining_ways=0):
	(benchmark_list,final_solution,apps_in_partial_solution)=normalize_solution(benchmark_set, solution,remaining_ways) 
	slowdown_vector=get_slowdown_vector(benchmark_list,final_solution,max_bandwidth)
	slowdown_vector = np.array(slowdown_vector, dtype='float')
	fairness=1-np.std(slowdown_vector)/np.mean(slowdown_vector)
	del slowdown_vector
	return fairness


def cov_unfairness_metric(benchmark_set, solution, max_bandwidth=float('Inf'), remaining_ways=0):
	(benchmark_list,final_solution,apps_in_partial_solution)=normalize_solution(benchmark_set, solution,remaining_ways) 
	slowdown_vector=get_slowdown_vector(benchmark_list,final_solution,max_bandwidth)
	slowdown_vector = np.array(slowdown_vector, dtype='float')
	unfairness=np.std(slowdown_vector)/np.mean(slowdown_vector)
	del slowdown_vector
	return unfairness

def jain_fairness(benchmark_set, solution, max_bandwidth=float('Inf'), remaining_ways=0):
	cov=cov_unfairness_metric(benchmark_set, solution, max_bandwidth, remaining_ways)
	return 1.0/(1.0+cov**2)


# Params: Benchmarks data, solution dictionary, max bandwidth and remaining ways (to know if it is a partial solution, used for bounding)
# Return: Unfairness and -throughput tuple calculated.
def max_slowdown_unfairness(benchmark_set, solution, max_bandwidth=float('Inf'), remaining_ways=0):
	(benchmark_list,final_solution,apps_in_partial_solution)=normalize_solution(benchmark_set, solution,remaining_ways) 
	throughput_value = 0	
	
	slowdown_vector=get_slowdown_vector(benchmark_list,final_solution,max_bandwidth)
	

	max_slowdown=max(slowdown_vector)
	## Pick min that cannot grow with new slowdown values...
	if remaining_ways>0:
		## Let's assume that min slowdown unfairness will not get any smaller
		unfairness_value = max_slowdown / min(slowdown_vector[0:len(solution)])
	else:
		unfairness_value = max_slowdown / min(slowdown_vector)

	unfairness_value=round(unfairness_value,3)
	max_slowdown=round(max_slowdown,3)

	return (max_slowdown,unfairness_value)

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
		benchmark_dict[benchmark_name]=benchmark_dframe


	## Traverse workloads first
	i=0
	for workload in workloads:
		for benchmark_name in workload:
			if not benchmark_name in benchmarks_list:
				print("Can't find benchmark name (%s) in global table" % benchmark_name)
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
		benchmark_dict[benchmark_name]=benchmark_dframe


	## Traverse workloads first
	i=0
	for workload in workloads:
		for benchmark_name in workload:
			if not benchmark_name in benchmarks_list:
				print("Can't find benchmark name (%s) in global table" % benchmark_name)
				return None
			## Get right Dframe
			benchmark_dframe=benchmark_dict[benchmark_name]
			## Safe to continue	
			workloads_table[i].append(App(benchmark_name,benchmark_dframe))
		i += 1		

	return (workloads_table,nr_ways)


def get_workloads_table_from_list(summary_csv_path, workloads, separator=None):
		return get_workloads_table_from_list_gen(summary_csv_path, workloads)[0]


# To do analysis manually with the simulator
def get_application_info_from_file(summary_csv_path):
	summary_dframe=pd.read_csv(summary_csv_path,sep=',')

	# Get list of benchmarks
	benchmarks_list = summary_dframe["BENCH"].unique()
	apps=[]

	# Get info for each benchmark in the global table
	for benchmark_name in benchmarks_list:
		# Get benchmark_name's DataFrame using NR_WAYS as index
		benchmark_dframe = summary_dframe.loc[summary_dframe["BENCH"] == benchmark_name].set_index("NR_WAYS")
		
		## Calculate max IPS and slowdown and save on benchmark_dframe column
		ipc_vals=benchmark_dframe[["ipc"]].values
		ipc_max = max(ipc_vals)[0]
		benchmark_dframe["slowdown"] = ipc_max/benchmark_dframe["ipc"]
		nr_ways=len(ipc_vals)

		apps.append(App(benchmark_name,benchmark_dframe))

	return apps

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
	if isinstance(curves, map):
		curves = list(curves)

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
		nr_fair_ways_per_benchmark = available_ways // n_apps
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


def number_of_nodes_partitioning(nr_ways,nr_partitions):
	assert nr_partitions<=nr_ways

	if nr_partitions==1 or nr_ways==nr_partitions:
		return 1		
	elif nr_ways==nr_partitions+1: ## Just one app can take the extra way...
		return nr_partitions
	else:
		## Recursive case
		total_nodes= 0

		for i in range(1,nr_ways-(nr_partitions-1)+1):
			total_nodes+=number_of_nodes_partitioning(nr_ways-i,nr_partitions-1)+1

		return total_nodes 


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


def generate_possible_clusters(collection,nr_ways,max_cluster=0):
	if max_cluster==0:
		max_cluster=len(collection)

	if len(collection) == 1:
		yield [ collection ]
		return

	first = collection[0]
	for smaller in generate_possible_clusters(collection[1:],nr_ways,max_cluster):
		size_clustering=len(smaller)

		##Unfeasible solution
		if size_clustering>nr_ways:
			continue

		# insert `first` in each of the subpartition's subsets
		# This does not make the clustering any bigger
		for n, subset in enumerate(smaller):
			if len(subset)< max_cluster:
				yield smaller[:n] + [[ first ] + list(subset)]  + smaller[n+1:]
		
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

def generate_possible_mappings(collection,min_items,max_items,max_groups,level=1):
	if len(collection) == 1:
		yield [ collection ]
		return

	first = collection[0]
	for smaller in generate_possible_mappings(collection[1:],min_items,max_items,max_groups,level+1):  
		# insert `first` in each of the subpartition's subsets
		#if level==1:
		#    set_trace()
		if len(smaller)>max_groups:
			continue
		for n, subset in enumerate(smaller):
			if len(subset)< max_items:
				left=smaller[:n]
				right=smaller[n+1:]
				unfeasible=False
				if level==1:
					i=0
					m=len(left)
					while not unfeasible and i<m:
						if len(left[i])<min_items:
							unfeasible=True
						i=i+1
					i=0
					m=len(right)
					while not unfeasible and i<m:
						if len(right[i])<min_items:
							unfeasible=True
						i=i+1
				if not unfeasible:                    
					yield left[:n] + [[ first ] + list(subset)]  + right
		
				# put `first` in its own subset if that it is allowed
		if  (len(smaller)<max_groups) and ((level==1 and min_items==1) or (level>1)):
			yield [[ first ]] + smaller

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
def yu_petrov_algorithm(workload, bandwidth_reduction_table, nr_ways, threshold = 0):
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
	partition_sizes = list(map(int, partition_sizes))
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
			nr_extra_ways-=1		
		else:
			ways_this_app=nr_fair_ways

		# Round down the number of ways assigned
		ways_this_app = int(ways_this_app)
		solution.append(ways_this_app)
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


class EngineCommunicator(object):
	def __init__(self, interface='tcp://*', identity=None):
		self._ctx = zmq.Context()
		self.pub = self._ctx.socket(zmq.PUB)
		self.sub = self._ctx.socket(zmq.SUB)
		#self.sub.setsockopt( zmq.LINGER, 0 )

		# configure sockets
		self.identity = identity or bytes(uuid.uuid4())
		self.sub.SUBSCRIBE = b''

		# bind to ports
		pub_port = self.pub.bind_to_random_port(interface)
		self.pub_url = interface+":%i"%pub_port
		# guess first public IP from socket
		self.location = socket.gethostbyname_ex(socket.gethostname())[-1][0]
		self.peers = {}
   
	def __del__(self):
		self.pub.close()
		self.sub.close()
		self._ctx.term()
    
	def info(self):
		"""return the connection info for this object's sockets."""
		return (self.identity, self.pub_url, self.location)
    
	def connect(self, peers, master=False):
		"""connect to peers.  `peers` will be a dict of 4-tuples, keyed by name.
		{peer : (ident, addr, pub_addr, location)}
		where peer is the name, ident is the XREP identity, addr,pub_addr are the
		"""
		for peer, (ident, pub_url, location) in peers.items():
			self.peers[peer] = ident
			if ident != self.identity or master:
				self.sub.connect(disambiguate_url(pub_url, location))

	def publish(self, msg, flags=0):
		self.pub.send_pyobj(msg,flags=flags)

	def consume(self, flags=0):
		return self.sub.recv_pyobj(flags=flags)

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

	if "topology" in key_vals:
		select_topology(key_vals["topology"])

	if "broadcast" in key_vals and key_vals["broadcast"]:
		com=EngineCommunicator()
		gbl_dict["com"]=com
		return com.info()

	#ttps://stackoverflow.com/questions/3061/calling-a-function-of-a-module-by-using-its-name-a-string
	#if "cost_function" in key_vals:
		#retrieve function symbol

	return pid	


def start_communicator(peers):
	global gbl_dict

	## Retrieve communicator object
	com=gbl_dict["com"]
	## Connect to peers
	com.connect(peers)

	## Create and start background thread
	cthread=ComThread(com)
	cthread.start()

	##Store cthread ref in dict
	gbl_dict["cthread"]=cthread

def stop_communicator(engine_id):
	global gbl_dict

	if "com" in gbl_dict:
		## Send empty list to close communication
		com=gbl_dict["com"]
		cthread=gbl_dict["cthread"]

		## This is no longer used as the master notifies itself
		#if engine_id==0:
		#	com.publish((None,None))
			
		## Wait for cthread to finish
		cthread.stop()
		cthread.join()

		del gbl_dict["com"]
		del gbl_dict["cthread"]
		del com

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
	kwargs.setdefault("paraver",False)

	metadatum = {}

	if "user_options" in kwargs:
		## Merge dict with user options... [Overriding defaults...]
		kwargs.update(kwargs["user_options"])

	## To make things simpler
	bound=kwargs["bound"]
	print_times=kwargs["print_times"]
	asyncr=kwargs["async"]
	chunk=kwargs["chunk"]
	prepruning=kwargs["prepruning"]
	postpruning=kwargs["postpruning"]
	parallel_debug=kwargs["parallel_debug"]
	bw_model=kwargs["bw_model"]
	paraver =  kwargs["paraver"]

	## Disable options that do not make sense
	if parallel_debug:
		asyncr=False

	if not bound:
		prepruning=False

	if not asyncr:
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
		if kwargs["hint_ucp_slowdown"] and not clustering:
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
		lview.block = not asyncr		

		# DEBUG code to check if variables are being set all right
		# ret = lview.map(lambda x: get_global_property(x),[3]*nr_engines, block=True)
		# print ret
		# rc.close()
		# exit(0)
		
		start = time.time()

		if log_enabled:
			if print_times and paraver:
				print("Activated paraver trace generation", file=sys.stderr)
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
					if parallel_debug:
						promising= map(lambda x: subsol_is_promising_parallel(best_cost,x), subsols)
					else:	
						promising= lview.map(lambda x: subsol_is_promising_parallel(best_cost,x), subsols,chunksize=chunkval)

					if asyncr:
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

			print("Unroll/Preprune/Merge time: %s (Pruning rate: %s out of %d)%s" % (sim_time,len(clusters),nr_initial_solutions, str_merge), file=sys.stderr)

		if postpruning:
			if log_enabled:
				start = datetime.datetime.utcnow()

			old_clusters=len(clusters)
			clusters=unroll_solutions(nsols,nr_ways,len(workload),kwargs["low_threshold"],clusters)


			if log_enabled:
				end = datetime.datetime.utcnow()
				if print_times:
					print("Postpruning split: %s (Goes from %d to %d solutions)" % (str(round((end - start).total_seconds(),4)),old_clusters,len(clusters)), file=sys.stderr)


			if paraver:
				(start_micros,end_micros)=get_start_end_micros(start,end,start_datet)
				trace.append(["1","1","1","1","1",str(start_micros),str(end_micros),"7"])
				# Write header and trace to file


		### Purely parallel stage
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
					print(sol)
			else:			
				local_node_solutions = lview.map(lambda x: get_optimal_schedule_aux_multi_parallel(best_cost,x,(min_idx!=0) or (not prepruning) or (postpruning)), clusters[min_idx:max_idx])	

			if asyncr: 
				## Now local_node_solutions is an AsyncMapResult object. But it can be iterated.
				## Wait for completion.... invoking get (later we'll use metadata!)
				d=local_node_solutions.get()

			## Debugging code
			if log_enabled:
				end = time.time()
				if print_times:
					print("Parallel computation: %s" % str(round(end - start,4)), file=sys.stderr)

				if asyncr:
					## https://ipyparallel.readthedocs.io/en/latest/details.html#extended-interface
					##Extract metadata
					metadata=local_node_solutions.metadata

					for i,stats in enumerate(metadata):
						computation_time=(stats['completed']-stats['started']).total_seconds()
						#transfer_time=(stats['started']-stats['submitted']).total_seconds()+(stats['received']-stats['completed']).total_seconds()
						total_time=(stats['received']-stats['submitted']).total_seconds()
						if print_times:
							print(clusters[low_idx+i],"%i,%i,%.4f,%.4f" % (i,stats['engine_id'],computation_time,total_time),local_node_solutions[i],stats["stdout"], file=sys.stderr)
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
					print("Sequential reduction: %s" % str(round((end - start).total_seconds(),4)), file=sys.stderr)

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


################################################################################
##                           FUNCTIONS FOR PARAVER-TRACING INSTRUMENTATION                          ##
################################################################################


### Code to refactor paraver instrumentation code
def add_trace_item(trace,stats,start_datet,event_type,include_master=False):
	offset=2 if include_master else 1
	start_micros=int((stats["started"]-start_datet).total_seconds()*1000000)
	end_micros=int((stats["completed"]-start_datet).total_seconds()*1000000)
	trace.append(["1","%i"% (stats['engine_id']+offset),"1","1","%i"%(stats['engine_id']+offset),str(start_micros),str(end_micros),event_type])

def add_trace_metadata(trace,metadata,start_datet,event_type,include_master=False):
	for stats in metadata: 
		add_trace_item(trace,stats,start_datet,event_type,include_master)

def get_trace_timestamp():
	return datetime.datetime.utcnow()

def add_trace_seq_item(trace,start,end,start_datet,event_type):
	(start_micros,end_micros)=get_start_end_micros(start,end,start_datet)
	trace.append(["1","%i"%1,"1","1","%i"%(1),str(start_micros),str(end_micros),event_type])

def generate_start_datet():
	unlocalized_start = pytz.utc.localize(datetime.datetime.utcnow())
	start_datet = unlocalized_start.astimezone(pytz.timezone("utc"))
	return start_datet

def generate_paraver_header(start_datet,nr_engines):
	end = datetime.datetime.utcnow()
	# Get sequential time execution
	end_seqtime=pytz.utc.localize(end)
	end_seqtime=end_seqtime.astimezone(pytz.timezone("utc"))
	end_micros=int((end_seqtime-start_datet).total_seconds()*1000000)

	return "#Paraver ({}/{}/{} at {}:{}):{}:1({}):1:1({}:1)\n".format(start_datet.day,start_datet.month,start_datet.year,start_datet.hour,start_datet.minute,end_micros,nr_engines,nr_engines)



################################################################################
##                 SEQUENTIAL BEST-FIRST B&B SEARCH APPROACH                  ##
################################################################################


## Entail separate set for sequence and queues...
def bb_insert_node(seq, keys, item, keyfunc=lambda v: v,key=None):
	"""Insert an item into a sorted list using a separate corresponding
		sorted keys list and a keyfunc() to extract the key from each item.

	Based on insert() method in SortedCollection recipe:
	http://code.activestate.com/recipes/577197-sortedcollection/
	"""
	if not key is None:
		k=key
	else:
		k = keyfunc(item)  # Get key.
	i = bisect_left(keys, k)  # Determine where to insert item.
	keys.insert(i, k)  # Insert key of item to keys list.
	seq.insert(i, item)  # Insert the item itself in the corresponding place.

## Extra parameter lower bounds 
def get_optimal_schedule_bf_seq(benchmark_set, cost_function, max_bandwidth, default_cost, op, total_nr_ways, partial_solutions, lower_bounds, bound):

	# To determine n_clusters
	if type(benchmark_set) is list:
		cluster_list=benchmark_set
	else:
		cluster_list=benchmark_set[1] ## Use cluster rather than that

	n_clusters=len(cluster_list)

	pq=[] ## Priority queue
	lb=[] ## Lower bound associated with prio queue
	best_cost=default_cost
	best_solution=None
	total_nodes_expanded=0

	## Need to insert root of the tree
	if len(partial_solutions)==0:
		pq.append([])
		lb.append(1) ## It does not really matter
		total_nodes_expanded=1
	else:
		## Transfer solutions thus far to pool
		for i in range(len(partial_solutions)):
			bb_insert_node(pq,lb,partial_solutions[i],key=lower_bounds[i])

	## While list is not empty
	while (pq):
		partial_solution=pq.pop(0)
		## Do nothing with ps_lower_bound here
		ps_lower_bound=lb.pop(0)
		nr_ways=total_nr_ways-sum(partial_solution)
		## We visit a new node 
		level=len(partial_solution)
		leaf_node=False

		if level == n_clusters - 1: # Trivial Base case
			partial_solution.append(nr_ways) # Being the last app, we give you all the available 'slots'
			## Potential optimization (lower bound should match cost in this case, right!)
			this_cost = cost_function(benchmark_set, partial_solution, max_bandwidth)
			leaf_node=True
			total_nodes_expanded+=1
		elif nr_ways == n_clusters - level: ## Optimization, one way for each app....
		##Easy ...
			for i in range(nr_ways):
				partial_solution.append(1)
			this_cost = cost_function(benchmark_set, partial_solution, max_bandwidth)
			leaf_node=True
			total_nodes_expanded+=1
		else: 
		# Non leaf node
			for nr_ways_assigned in range(1, nr_ways - (n_clusters - level - 1) + 1):
				new_partial_solution = list(partial_solution)
				new_partial_solution.append(nr_ways_assigned)
				remaining_ways = nr_ways - nr_ways_assigned

				this_lower_bound=cost_function(benchmark_set, new_partial_solution, max_bandwidth, remaining_ways)

				total_nodes_expanded+=1

				## Prune
				if bound and op(best_cost, this_lower_bound):
					continue
				else:
					## Enqueue
					bb_insert_node(pq,lb,new_partial_solution,key=this_lower_bound)

			depth_first = True
		if leaf_node and op(this_cost, best_cost):
			best_cost = this_cost
			best_solution = partial_solution

			## Purge unpromising nodes in the priority queue
			for lower_bound in reversed(lb):
				if  not op(lower_bound, best_cost):
					lb.pop()
					pq.pop()
				else:
					break ## Sorted list!

	return (best_cost, best_solution, total_nodes_expanded)


## Maintain separate lists for lbs and nodes; depths to control indexing
def bb_insert_node_hyb(pq, lbs, d_pq, d_lb, d_depths, item, lb):
	"""
	key is now the depth of the node
	lb is the extra argument to insert in lbs
	item is the node (sol/partial_sol)
	"""
	i = bisect_left(lbs, lb)  # Determine where to insert item.
	pq.insert(i, item)  # Insert the item itself in the corresponding place.
	lbs.insert(i, lb)

	depth = -(len(item))
	j = bisect_left(d_depths, depth)
	d_depths.insert(j, depth)  # Insert key of item to keys list.
	d_pq.insert(j, item)  # Insert the item itself in the corresponding place.
	d_lb.insert(j, lb)

def bb_pop_first_node_hyb(pq, lbs, d_pq, d_lb, d_depths, depth_first):
	if depth_first:
		# Obtain first element from the depth ordered
		partial_solution = d_pq.pop(0)
		d_lb.pop(0)
		d_depths.pop(0)
		ix = pq.index(partial_solution)
		pq.pop(ix)
		lbs.pop(ix)
	else:
		partial_solution = pq.pop(0)
		lbs.pop(0)
		jx = d_pq.index(partial_solution)
		d_pq.pop(jx)
		d_lb.pop(jx)
		d_depths.pop(jx)
	return partial_solution

def purge_unpromising(pq, lbs, d_pq, d_lb, d_depths, best_cost, op):
	## Purge unpromising nodes in the depth queue
	i = len(d_lb) - 1
	for lower_bound in reversed(d_lb):
		if not op(lower_bound, best_cost):
			sol = d_pq.pop(i)
			d_depths.pop(i)
			d_lb.pop(i)
			ix = pq.index(sol)
			pq.pop(ix)
			lbs.pop(ix)
		i -= 1


## Extra parameter lower bounds
def get_optimal_schedule_bf_seq_hybrid(benchmark_set, cost_function, max_bandwidth, default_cost, op, total_nr_ways,
								partial_solutions, lower_bounds, bound):
	# To determine n_clusters
	if type(benchmark_set) is list:
		cluster_list = benchmark_set
	else:
		cluster_list = benchmark_set[1]  ## Use cluster rather than that

	n_clusters = len(cluster_list)

	depth_first = False
	first_division = False

	# Indexed by depth
	d_pq = []  ## Priority queue
	d_lb = []  ## Lower bound associated with prio queue
	d_depths = [] ## Depth list so the depth of each item is tracked so the other lists are indexed
	# Indexed by lb
	pq = []
	lbs = []

	best_cost = default_cost
	best_solution = None
	total_nodes_expanded = 0

	## Need to insert root of the tree
	if len(partial_solutions) == 0:
		d_pq.append([])
		d_lb.append(1)  ## It does not really matter
		d_depths.append(0)
		pq.append([])
		lbs.append(1)
		total_nodes_expanded = 1
	else:
		## Transfer solutions thus far to pool
		for i in range(len(partial_solutions)):
			bb_insert_node_hyb(pq, lbs, d_pq, d_lb, d_depths, partial_solutions[i], lower_bounds[i])

	## While list is not empty
	while (pq):
		partial_solution = bb_pop_first_node_hyb(pq, lbs, d_pq, d_lb, d_depths, depth_first)
		nr_ways = total_nr_ways - sum(partial_solution)
		## We visit a new node

		level = len(partial_solution)
		leaf_node = False

		# Wait until we visit a node in the second level with best-first to start DF to arrive to the leaves
		if level > 1:
			first_division = True
		# After chosing the next node with BF we enable depth_first to arrive faster to a leaf (update upper bound more freq.)
		if first_division and not depth_first:
			depth_first = True

		if level == n_clusters - 1:  # Trivial Base case
			partial_solution.append(nr_ways)  # Being the last app, we give you all the available 'slots'
			this_cost = cost_function(benchmark_set, partial_solution, max_bandwidth)
			leaf_node = True
			total_nodes_expanded += 1
		elif nr_ways == n_clusters - level:  ## Optimization, one way for each app....
			##Easy ...
			for i in range(nr_ways):
				partial_solution.append(1)
			this_cost = cost_function(benchmark_set, partial_solution, max_bandwidth)
			leaf_node = True
			total_nodes_expanded += 1
		else:
			# Non leaf node
			for nr_ways_assigned in range(1, nr_ways - (n_clusters - level - 1) + 1):
				new_partial_solution = list(partial_solution)
				new_partial_solution.append(nr_ways_assigned)
				remaining_ways = nr_ways - nr_ways_assigned

				this_lower_bound = cost_function(benchmark_set, new_partial_solution, max_bandwidth, remaining_ways)

				total_nodes_expanded += 1
				depth = len(new_partial_solution)
				## Prune
				if bound and op(best_cost, this_lower_bound):
					continue
				else:
					bb_insert_node_hyb(pq, lbs, d_pq, d_lb, d_depths, new_partial_solution, this_lower_bound)

		if leaf_node:
			# Disable the depth_first (next node will be picked with BF)
			depth_first = False
			if op(this_cost, best_cost):
				# We have reached a leaf node that improves the upper-bound
				best_cost = this_cost
				best_solution = partial_solution
				purge_unpromising(pq, lbs, d_pq, d_lb, d_depths, best_cost, op)

	return (best_cost, best_solution, total_nodes_expanded)

################################################################################
##                 PARALLEL BEST-FIRST B&B SEARCH APPROACH                    ##
################################################################################

## Background thread to process data from socket
class ComThread(threading.Thread):
	def __init__(self,communicator):
		## Queue of items processed from the socket
		threading.Thread.__init__(self)
		self._stop_event = threading.Event() 
		self.items=[]
		self.lock=threading.Lock()
		self.communicator=communicator

	def stop(self):
		self._stop_event.set()
		#self.communicator.sub.close()
	
	def stopped(self):
		return self._stop_event.is_set()

	## Thread's main()
	def run(self):
		while not self._stop_event.is_set( ):
			item=None
			try:
				item = self.communicator.consume()
				#print "Item received", item
				if item[0] is None:
					break
			except Exception as e: 
				print(e)
				break 
			
			if item is None:
				return 
			self.lock.acquire()
			self.items.append(item)
			self.lock.release()

	## Master uses lower_bound==None
	def try_to_update_solution(self,local_solution,lower_bound):
		## Dirty read first: optimization for the most common case
		if len(self.items)==0:
			return (False,False)
		
		remote_sols=self.consume_items()

		found_better=False
		master=lower_bound is None ## To be considered later

		## Process most recent solutions first
		for remote_sol in reversed(remote_sols):
			(sol,cost)=remote_sol
			## I am going to update the cost only
			if (local_solution.is_worse(cost)):
				found_better=True
				if master:
					local_solution.update(sol,cost,False) ## Update sol but do not broadcast sol...
				else:
					local_solution.update_cost(cost)
		
		if master:
			return (found_better,False)
		else:
			return (found_better,local_solution.is_better(lower_bound))


	## Safe method to consume items processed by the thread
	def consume_items(self):
		ret_items=[]
		self.lock.acquire()
		for item in self.items:
			ret_items.append(item)

		del self.items[:] ## clear() method is python 3
	
		self.lock.release()

		return ret_items


## Abstraction to maintain solution neatly
class SolutionManager:
	def __init__(self,sol=None,cost=None,op=operator.lt,communicator=None):
		self.sol=sol
		self.cost=cost
		self.op=op
		self.communicator=communicator

	def get_solution(self):
		return self.solution

	def get_cost(self):
		return self.cost

	def update_cost(self,new_cost):
		self.cost=new_cost
		self.sol=None ## Not worried about cost

	def update(self,solution,cost,broadcast=True):
		self.sol=solution
		self.cost=cost

		if broadcast and self.communicator:
			self.communicator.publish((self.sol,self.cost))
	
	def is_worse(self,this_cost):
		## this cost cannot be none
		return self.op(this_cost, self.cost)

	def is_better(self,this_cost):
		## this cost cannot be none
		return self.op(self.cost,this_cost)


## Function to update global properties after the node processing
def bb_process_node_result(task,props):
	if task.cancelled():
		#TODO consider what to do with cancelled tasks regarding metadata
		#meta=task.metadata
		#print meta.keys()
		#print "Task was cancelled in process"
		return False

	#node_queue
	#lb_queue
	pq=props["node_queue"]
	lbs=props["lb_queue"]
	best_solution=props["best_solution"]
	total_nodes=props["total_nodes"]
	better_solution=True

	## Get result
	(nodes_lbs,this_sol,this_cost,node_count)=task.get()
	#print(task.stdout)

	## Get metadata
	if "trace" in props:
		add_trace_item(props["trace"],task.metadata,props["start_datet"],"1",include_master=props["trace_master"])

	## Update global node count 
	props["total_nodes"]=total_nodes+node_count

	## Found a better solution
	if (not this_sol is None) and best_solution.is_worse(this_cost):
			## Update dict as well
			best_solution.update(this_sol, this_cost, broadcast=False) 

			#print "New solution found:", this_sol,this_cost
			better_solution=True
			## Purge unpromising nodes in the priority queue
			for lower_bound in reversed(lbs):
				if  best_solution.is_better(lower_bound):
					lbs.pop()
					pq.pop()
				else:
					break ## Sorted list!

	## Add to common pool 
	if better_solution:
		for node,lb in nodes_lbs:
			bb_insert_node(pq,lbs,node,key=lb)
	else:
		for node,lb in nodes_lbs:
			##Add if potentially better
			if  best_solution.is_worse(lb): 
				bb_insert_node(pq,lbs,node,key=lb)			

	return better_solution


def bb_is_leaf_node(partial_solution, nr_clusters, nr_ways):
	level=len(partial_solution)
	remain_nr_ways=nr_ways-sum(partial_solution)

	return (level == nr_clusters - 1) or (remain_nr_ways == (nr_clusters - level))


## Returns almost the same thing that bb_node_processing_aux
## (best_solution,best_cost,total_nodes_expanded)	
def bb_process_leaf_node(benchmark_set, cost_function, max_bandwidth, default_cost, op, nr_ways, partial_solution,nr_clusters,**kwargs):
	level=len(partial_solution)
	remain_nr_ways=nr_ways-sum(partial_solution)

	best_solution=None
	best_cost=default_cost
	remain_nr_ways=nr_ways-sum(partial_solution)
	total_nodes_expanded=0
	leaf_node=False

	if level == nr_clusters - 1: # Trivial Base case
		partial_solution.append(remain_nr_ways) # Being the last app, we give you all the available 'slots'
		## Potential optimization (lower bound should match cost in this case, right!)
		this_cost = cost_function(benchmark_set, partial_solution, max_bandwidth)
		leaf_node=True
		total_nodes_expanded+=1
	elif remain_nr_ways == nr_clusters - level: ## Optimization, one way for each app....
	##Easy ...
		for i in range(remain_nr_ways):
			partial_solution.append(1)
		this_cost = cost_function(benchmark_set, partial_solution, max_bandwidth)
		leaf_node=True
		total_nodes_expanded+=1
	

	if leaf_node and op(this_cost, best_cost):
		best_cost = this_cost
		best_solution = partial_solution

	return (best_solution,best_cost,total_nodes_expanded)	


def bb_node_processing_aux(benchmark_set, cost_function, max_bandwidth, engine_solution, nr_ways, partial_solution,nr_clusters,lower_bound,**kwargs):
	best_solution=None
	remain_nr_ways=nr_ways-sum(partial_solution)
	total_nodes_expanded=0

	if "cthread" in kwargs:
		cthread=kwargs["cthread"]
	else:
		cthread=None
	
	## Local queue for this node: tuples -> (subnode,lb)
	local_queue=[]

	## We visit a new node 
	level=len(partial_solution)
	new_solution=None


	if bb_is_leaf_node(partial_solution,nr_clusters,nr_ways):
		## Process leaf node
		(new_solution,this_cost,total_nodes_expanded)=bb_process_leaf_node(benchmark_set, cost_function, max_bandwidth, engine_solution.cost, engine_solution.op, nr_ways, partial_solution,nr_clusters)
		## Update better solution
		if (not new_solution is None):
			engine_solution.update(new_solution,this_cost)
			best_solution = new_solution		
	else: 
	# Non leaf node
		for nr_ways_assigned in range(1, remain_nr_ways - (nr_clusters - level - 1) + 1):
			new_partial_solution = list(partial_solution)
			new_partial_solution.append(nr_ways_assigned)
			remaining_ways = remain_nr_ways - nr_ways_assigned

			## Update local solution
			if cthread is not None:
				(sol_updated,prune_node)=cthread.try_to_update_solution(engine_solution,lower_bound)
				## Exit right away
				if prune_node:
					return ([],best_solution,engine_solution.cost,total_nodes_expanded)

			## Optimization: Detect terminal nodes automatically and record solution if better
			if bb_is_leaf_node(new_partial_solution,nr_clusters,nr_ways):
				## Process leaf node
				(new_solution,this_cost,count)=bb_process_leaf_node(benchmark_set, cost_function, max_bandwidth, engine_solution.cost, engine_solution.op, nr_ways, new_partial_solution, nr_clusters)

				if count==0:
					raise Exception('pete')

				## Update best cost if necessary
				total_nodes_expanded=+1
				if not new_solution is None:
					## Update local incumbent 
					engine_solution.update(new_solution,this_cost)
					best_solution = new_solution

			else:
				this_lower_bound=cost_function(benchmark_set, new_partial_solution, max_bandwidth, remaining_ways)

				total_nodes_expanded+=1

				## Prune
				if engine_solution.is_better(this_lower_bound):
					continue
				else:
					## Enqueue
					local_queue.append((new_partial_solution,this_lower_bound))


	## Better solution found when processing leaf node
	if (not best_solution is None):
		local_queue=list(filter(lambda item: engine_solution.is_worse(item[1]),local_queue))

	return (local_queue,best_solution,engine_solution.cost,total_nodes_expanded)


## Lower bound is not part of the subnode specification
def bb_break_into_subnodes(node,nr_clusters, total_nr_ways, max_children):
	if bb_is_leaf_node(node,nr_clusters,total_nr_ways):
		return [(node,0,0)]
	else:
		level=len(node)
		nr_remaining_ways=total_nr_ways-sum(node)
		num_children=nr_remaining_ways- (nr_clusters - level - 1) 

		if max_children<=0 or num_children<=max_children:
		##Subnode is the full node
			return [(node,1,num_children+1)]
		else:
			subnodes=[]
			## Loop
			low_limit=1
			## Round up -------> in Python2 it was coded as (express.)/max_children, which gives only the int part of division, so it always rounds down 
			nr_partitions=(num_children+(max_children-1))//max_children
			for i in range(nr_partitions):
				high_limit=min(low_limit+max_children,num_children+1)
				subnodes.append((node,low_limit,high_limit))
				low_limit=high_limit

			return subnodes

def bb_subnode_processing(benchmark_set, cost_function, max_bandwidth, engine_solution, nr_ways, subnode ,nr_clusters, max_children,lower_bound,**kwargs):
	best_solution=None
	total_nodes_expanded=0
	(partial_solution, low_limit, high_limit)=subnode
	remain_nr_ways=nr_ways-sum(partial_solution)

	if "cthread" in kwargs:
		cthread=kwargs["cthread"]
	else:
		cthread=None

	## Local queue for this node: tuples -> (subnode,lb)
	local_queue=[]

	## We visit a new node 
	level=len(partial_solution)
	new_solution=None


	if bb_is_leaf_node(partial_solution,nr_clusters,nr_ways):
		## Process leaf node
		(new_solution,this_cost,total_nodes_expanded)=bb_process_leaf_node(benchmark_set, cost_function, max_bandwidth, engine_solution.cost, engine_solution.op, nr_ways, partial_solution,nr_clusters)

		## Update better solution
		if (not new_solution is None):
			engine_solution.update(new_solution,this_cost)
			best_solution = new_solution	
	else: 
	# Non leaf node
		for nr_ways_assigned in range(low_limit, high_limit):
			new_partial_solution = list(partial_solution)
			new_partial_solution.append(nr_ways_assigned)
			remaining_ways = remain_nr_ways - nr_ways_assigned

			## Update local solution
			if cthread is not None:
				(sol_updated,prune_node)=cthread.try_to_update_solution(engine_solution,lower_bound)
#				if sol_updated:
#					print "Solution updated"
				## Exit right away
				if prune_node:
					return ([],best_solution,engine_solution.cost,total_nodes_expanded)

			## Optimization: Detect terminal nodes automatically and record solution if better
			if bb_is_leaf_node(new_partial_solution,nr_clusters,nr_ways):
					## Process leaf node
					(new_solution,this_cost,count)=bb_process_leaf_node(benchmark_set, cost_function, max_bandwidth, engine_solution.cost, engine_solution.op, nr_ways, new_partial_solution, nr_clusters)

					if count==0:
						raise Exception('pete')

					## Update best cost if necessary
					total_nodes_expanded=+1
					if not new_solution is None:
						## Update local incumbent 
						engine_solution.update(new_solution,this_cost)
						best_solution = new_solution	

			else:
				this_lower_bound=cost_function(benchmark_set, new_partial_solution, max_bandwidth, remaining_ways)

				total_nodes_expanded+=1

				## Prune
				if engine_solution.is_better(this_lower_bound):
					continue
				else:
					## Break down into subnodes and enqueue....
					subnodes=bb_break_into_subnodes(new_partial_solution,nr_clusters,nr_ways,max_children)

					## Enqueue
					for subnode in subnodes:
						local_queue.append((subnode,this_lower_bound))

	## Better solution found when processing leaf node
	if (not best_solution is None):
		## Filter out unpromising subnodes
		local_queue=list(filter(lambda item: engine_solution.is_worse(item[1]),local_queue))

	return (local_queue,best_solution,engine_solution.cost,total_nodes_expanded)


def  bb_node_processing(arg):
	(node, best_cost,lower_bound)=arg
	## Get globals and invoke function 
	params=get_global_properties()
	engine_solution=params["engine_solution"]


	## Create best_solution object the very first time
	if engine_solution is None:
		## Connect with communicator in broadcast mode
		if "broadcast" in params and params["broadcast"]:
			com=params["com"]
		else:
			com=None
		engine_solution=SolutionManager(None,best_cost,params["op"],com)
		params["engine_solution"]=engine_solution

	else:
		## See if we can do it better
		if engine_solution.is_worse(best_cost):
			engine_solution.update_cost(best_cost)
	
	## Set up lower_bound
	params["lower_bound"]=lower_bound
	if params["max_children"]<0:
		params["partial_solution"]=node
		return bb_node_processing_aux(**params)
	else:
		params["subnode"]=node
		return bb_subnode_processing(**params)


def bb_prune_pending(pending_tasks, pending_lbs, best_solution):
	cancelled_tasks=0

	i=len(pending_lbs)-1
	for lower_bound in reversed(pending_lbs):
		if best_solution.is_better(lower_bound):
			task=pending_tasks[i]
			if not task.cancel():
				#print("Failed to cancel task")
				pass
			else:
				#meta=task.metadata
				#print meta.keys()
				#print meta
				#print "Task was cancelled in prune"
				pass
				##print "Cancelled task"
			cancelled_tasks=cancelled_tasks+1	
			## Removed from structure 
			pending_tasks.pop()
			pending_lbs.pop()
			i=i-1
		else:
			break ## Sorted list!

	return cancelled_tasks


nr_completed_bb_tasks=0
cond_var=threading.Condition()
tasks_completed=[]

def reset_completion_variables():
	global nr_completed_bb_tasks
	global cond_var
	global tasks_completed
	nr_completed_bb_tasks=0
	del tasks_completed[:]


def wait_until_task_completed():
	global nr_completed_bb_tasks
	global cond_var
	global tasks_completed
	completed=[]
	cond_var.acquire()

	while nr_completed_bb_tasks==0:
		cond_var.wait()

	## Reset bb task count
	nr_completed_bb_tasks=0
	for task in tasks_completed:
		completed.append(task)

	del tasks_completed[:] ## clear() method is python 3


	cond_var.release()

	return completed


def notify_task_completion(task):
	global nr_completed_bb_tasks
	global cond_var
	cond_var.acquire()

	## Reset bb task count
	nr_completed_bb_tasks=nr_completed_bb_tasks+1

	tasks_completed.append(task)

	## In this case there is only one thread, so wake it up if necessary
	cond_var.notify()

	cond_var.release()


## Note that the callback is executed in the current thread (launcher)
## if the task was cancelled. Otherwise the callback is executed by a different
## thread. That's why we need to verify whether the task has completed or not
def bb_task_completed(task):
	if task.cancelled():
		#print "Cancelled", threading.current_thread().ident 
		pass
	else:
		notify_task_completion(task)
		#print "Completed", threading.current_thread().ident 

def bb_process_pending(pending_tasks, pending_lbs,bb_props,op, completed=None):
	## Traverse pending tasks
	i=0
	deleted=0
	better_solution_found=False
	best_solution=bb_props["best_solution"]
	
	if completed is None:
		## First stage (check pending tasks that completed)
		for task in pending_tasks[:]:
			cur_cost=pending_lbs[i]
			if task.ready():
				if bb_process_node_result(task,bb_props):
					better_solution_found=True
				pending_tasks.remove(task)	
				del pending_lbs[i]
				deleted=deleted+1
				#print "task", i, "completed. Removing from list"	
			else:
				i=i+1
	else:
		for task in completed:
			i=pending_tasks.index(task)
			if i==-1 or not task.ready():
				print("Error")
				exit(1)
			if bb_process_node_result(task,bb_props):
				better_solution_found=True
			del pending_tasks[i]
			del pending_lbs[i]
			deleted=deleted+1

	if better_solution_found:
		## Cancel task if pending is worse than incumbent
		deleted+=bb_prune_pending(pending_tasks, pending_lbs, best_solution)

	return deleted


## Busqueda en anchura del arbol
def bb_generate_solutions_to_explore(nr_ways,nr_apps,target_count):
	
	queue=[[]]
	leaves=[]
	nodes_found=0

	while queue and nodes_found<target_count:
		node=queue.pop(0)
		level=len(node)
		nr_remaining_ways=nr_ways-sum(node)

		if nr_remaining_ways>0:
			## We will expand one node (so it is no longer in queue)
			nodes_found-=1
			for nr_ways_assigned in range(1, nr_remaining_ways- (nr_apps - level - 1) + 1):
				subsol=node+[nr_ways_assigned]
				queue.append(subsol)
				nodes_found+=1
		else:
			## Push in leaves
			leaves.append(node)
			nodes_found+=1

	return queue+leaves

## Single-workload pruning
def determine_lower_bound(benchmark_set, cost_function, max_bandwidth, default_cost, op, nr_ways, partial_solution,**kwargs):

	if type(partial_solution) is tuple:
		partial_solution=partial_solution[0] ## Discard size...

	remaining_ways = nr_ways - sum(partial_solution)
	new_partial_solution = list(partial_solution)
	lower_bound=cost_function(benchmark_set, new_partial_solution, max_bandwidth, remaining_ways)


	return (lower_bound,not op(default_cost, lower_bound))


## Wrappper for avoiding copy of parameters
def determine_lower_bound_parallel(default_cost,partial_solution):

	params=get_global_properties()
	params["default_cost"]=default_cost
	params["partial_solution"]=partial_solution

	return determine_lower_bound(**params)


## Extra parameter lower bounds 
def get_optimal_schedule_bf_par(benchmark_set, cost_function, max_bandwidth, default_cost, op, total_nr_ways, **kwargs):
	debug=False
	kwargs.setdefault("bound",True)
	kwargs.setdefault("print_times",False)
	kwargs.setdefault("bw_model","simple")
	kwargs.setdefault("paraver",False)
	kwargs.setdefault("trace_master",False)
	kwargs.setdefault("max_children",3)
	kwargs.setdefault("initial_load_factor",2)
	kwargs.setdefault("dyn_load_factor",2)
	kwargs.setdefault("broadcast",False)

	metadatum={}

	if "user_options" in kwargs:
		## Merge dict with user options... [Overriding defaults...]
		kwargs.update(kwargs["user_options"])

	## To make things simpler
	bound=kwargs["bound"]
	print_times=kwargs["print_times"]
	bw_model=kwargs["bw_model"]
	paraver =  kwargs["paraver"]
	max_children = kwargs["max_children"]
	initial_load_factor= kwargs["initial_load_factor"]
	dyn_load_factor= kwargs["dyn_load_factor"]
	log_enabled=paraver or print_times
	broadcast=kwargs["broadcast"]
	trace_master=kwargs["trace_master"]

	reset_completion_variables()

	# To determine n_clusters
	if type(benchmark_set) is list:
		cluster_list=benchmark_set
	else:
		cluster_list=benchmark_set[1] ## Use cluster rather than that

	nr_clusters=len(cluster_list)


	pq=[] ## Priority queue
	lb=[] ## Lower bound associated with prio queue
	pending_tasks=[]
	pending_lbs=[]

	#best_cost=default_cost
	best_solution=SolutionManager(None,default_cost,op)

	bb_props= { "node_queue": pq,
						"lb_queue": lb,
						"best_solution": best_solution,
						"total_nodes":0}

	rc = ipp.Client(timeout=30000) #debug=True)
	dview = rc[:]
	nr_engines=len(rc.ids)
	dict_globals={"benchmark_set":benchmark_set,
			"cost_function": cost_function,
			"max_bandwidth": max_bandwidth,
			"op": op,
			"nr_ways": total_nr_ways,
			"bound": bound,
			"bw_model": bw_model,
			"nr_clusters": nr_clusters,
			"max_children": max_children,
			"engine_solution": None, ## Critical not to reuse previous crap...
			"broadcast": broadcast
			}

	if log_enabled:
		trace=[]
		start_datet=generate_start_datet()
		start=get_trace_timestamp()			
		bb_props["trace"]=trace
		bb_props["start_datet"]=start_datet
		bb_props["trace_master"]=trace_master

	## Setup global variables for engines 
	# Copy global data
	ret=dview.apply_async(lambda x: set_global_properties(x),dict_globals)

	if broadcast:
		peers = ret.get_dict()
		## Create a local receiver of the solution
		com=EngineCommunicator()

		## Add this com to peers
		peers[nr_engines]=com.info()

		if debug:
			print(peers)

		## Start communicator remotely
		dview.apply_sync(lambda x: start_communicator(x),peers)

		## Connect to all publishers, including myself !
		com.connect(peers,master=True)

		## Initiate local socket thread
		cthread=ComThread(com)
		cthread.start()
	else:
		ret.get()
		cthread=None

	## Use Load Balanced View from now on in Async mode
	lview = rc.load_balanced_view()
	lview.block = False	

	## Generate initial nodes (to keep processors busy)
	initial_subsols=bb_generate_solutions_to_explore(total_nr_ways,nr_clusters,initial_load_factor*nr_engines)

	## Determine initial bounds remotely
	if trace_master:
		add_trace_seq_item(trace,start,get_trace_timestamp(),start_datet,"7")
		
	result=lview.map(lambda x: determine_lower_bound_parallel(best_solution.cost,x), initial_subsols)

	## Wait until everything is done
	seeds=result.get()

	if log_enabled:
		if trace_master:
			start=get_trace_timestamp()	
		
		add_trace_metadata(trace,result.metadata,start_datet,"6",include_master=trace_master)

	## We discard those intermediate nodes automatically expanded
	## Just account for those whose bound function was evaluated
	bb_props["total_nodes"]=len(initial_subsols)

	for subsol,seed in zip(initial_subsols,seeds):
		(lbound,promising)=seed
		if promising:
			if max_children<0:
				## Use nodes 
				bb_insert_node(pq,lb,subsol,key=lbound)
			else:
				## Break down into subnodes and enqueue....
				subnodes=bb_break_into_subnodes(subsol,nr_clusters,total_nr_ways,max_children)
				## Enqueue subnode
				for subnode in subnodes:
					bb_insert_node(pq,lb,subnode,key=lbound)

	## Populate 
	limit_queue=dyn_load_factor*nr_engines

	## While list is not empty
	while pq or pending_tasks:

		tasks_submitted=0

		## Retrieve remote solutions and prune if necessary
		if cthread is not None:
			(sol_updated,_)=cthread.try_to_update_solution(best_solution,None)
			
			## Cancel stuff remotely
			if sol_updated:
				#print "Trying to prune from master"
				bb_prune_pending(pending_tasks, pending_lbs,best_solution)

		## Send everything that we have until a certain limit
		while pq and len(pending_tasks)<=limit_queue:
			node=pq.pop(0)
			ps_lower_bound=lb.pop(0)
			
			## Using subnodes really
			if max_children>=0: 
				(subsol,ll,hl)=node
			else:
				subsol=node

			## Avoid processing leaf nodes remotely
			if bb_is_leaf_node(subsol,nr_clusters, total_nr_ways):
				(new_solution,this_cost,total_nodes_expanded)=bb_process_leaf_node(benchmark_set, cost_function, max_bandwidth, best_solution.cost, op, total_nr_ways, subsol,nr_clusters)

				bb_props["total_nodes"]=bb_props["total_nodes"]+1

				## Cancel tasks that can be prunned
				if not new_solution is None:
					## Update incumbent and cost
					best_solution.update(new_solution,this_cost, broadcast=False) 

					if debug:
						print("New solution found:", new_solution, this_cost)
					bb_prune_pending(pending_tasks, pending_lbs,best_solution)

			elif (ps_lower_bound is None) or best_solution.is_worse(ps_lower_bound):
				## No callbacks for now
				#print node,ps_lower_bound,bb_props["upper_bound"]
				task = lview.apply_async(bb_node_processing, (node,best_solution.cost,ps_lower_bound))

				task.add_done_callback(bb_task_completed)	

				## Add pending task
				bb_insert_node(pending_tasks,pending_lbs,task,key=ps_lower_bound)
				tasks_submitted=tasks_submitted+1
			else:
				if debug:
					print("Task pruned in origin")


		if tasks_submitted==0 and pending_tasks:
			if trace_master:
				add_trace_seq_item(trace,start,get_trace_timestamp(),start_datet,"7")
			
			completed=wait_until_task_completed()
			
			if trace_master:
				start=get_trace_timestamp()			
			
			if len(completed)>0:
				bb_process_pending(pending_tasks, pending_lbs, bb_props, op) # completed)

	if broadcast:
		## Send  terminating message
		com.publish((None,None))
		## Use non-blocking here to try to speed things up
		ret=lview.map(lambda x: stop_communicator(x),[i for i in range(nr_engines)])
		cthread.stop() 
		cthread.join()
		#ret.get()
		del com
	rc.close()

	if log_enabled:
		if trace_master:
			add_trace_seq_item(trace,start,get_trace_timestamp(),start_datet,"7")
		metadatum["header"]=generate_paraver_header(start_datet,nr_engines+trace_master)
		metadatum["trace"]=trace

	return (best_solution.cost, best_solution.sol, bb_props["total_nodes"],metadatum)



def get_optimal_schedule_bf(benchmark_set, cost_function, maximize, nr_ways, max_bandwidth=float('Inf'),**kwargs):

	kwargs.setdefault("bound",True)
	kwargs.setdefault("hint_ucp_slowdown",True)
	kwargs.setdefault("initial_bound",None)
	kwargs.setdefault("multiprocessing",False)
	kwargs.setdefault("user_options",None)
	kwargs.setdefault("bw_model","simple")
	kwargs.setdefault("hybrid",False)

	if "user_options" in kwargs:
		## Merge dict with user options... [Overriding defaults...]
		kwargs.update(kwargs["user_options"])

	bound=kwargs["bound"]
	bw_model=kwargs["bw_model"]

	## To maintain compatibility with the rest of the algorithms
	metadatum={}

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
		if kwargs["hint_ucp_slowdown"] and not clustering:
			heuristic=get_schedule_UCP_gen(workload,nr_ways,metric="slowdown")
		else:
			heuristic=get_equal_llc_schedule(len(cluster_list), nr_ways)

		default_cost=cost_function(benchmark_set,heuristic, max_bandwidth)

	## Sequential version
	if kwargs["multiprocessing"]:
		(best_cost, best_solution, total_branches, metadatum) = get_optimal_schedule_bf_par(benchmark_set, cost_function, max_bandwidth, default_cost, op, nr_ways, user_options=kwargs)
	elif kwargs["hybrid"]:
		(best_cost, best_solution, total_branches) = get_optimal_schedule_bf_seq_hybrid(benchmark_set, cost_function,
																				 max_bandwidth, default_cost, op,
																				 nr_ways, [], [], bound)
	else:
		(best_cost, best_solution, total_branches) = get_optimal_schedule_bf_seq(benchmark_set, cost_function, max_bandwidth, default_cost, op, nr_ways, [], [], bound)

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
def get_schedule_yu_petrov(workload, nr_ways, threshold = 0):
	bandwidth_reduction_table = build_bandwidth_reduction_table(workload, nr_ways)
	return yu_petrov_algorithm(workload, bandwidth_reduction_table, nr_ways, threshold)


# Weighted access of n mrcs with n float indexes
def access_indexes(mrcs, indexes):
	v_ns = []
	for i in range(len(indexes)):
		r = indexes[i] % 1
		int_i = int(indexes[i])
		#if r==0.0:
		#	v_ns.append(mrcs[i][0])
		#el
		if int_i >= len(mrcs[i]):
			v_ns.append(mrcs[i][int_i-1])
		elif (r == 0 or int_i+1 >= len(mrcs[i])):
			v_ns.append(mrcs[i][int_i])
		else:
			v_ns.append((1-r)*mrcs[i][int_i] + r*mrcs[i][int_i+1])
	return v_ns

def access_indexes_acc(mrcs, indexes):
    v_ns = []
    for i in range(len(indexes)):
        curve=mrcs[i]
        nr_ways=len(curve)
        r = indexes[i] % 1
        q = int(indexes[i])
            
        if q >= nr_ways:
            v_ns.append(curve[nr_ways-1])
        elif (r == 0) and (q>=1):
            v_ns.append(curve[q-1])
        elif q==0:
            value=(1-r)*curve[0] + r*curve[1]
            delta=value-curve[1]
            value=curve[0]+delta;
            
            if value<0:
                value=1
            v_ns.append(value)
        else:
            v_ns.append((1-r)*curve[q-1] + r*curve[q])
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
		mrc.append(sum(access_indexes_acc([mrc1,mrc2], [index1,index2])))
		b1 = int(round(index1))
		buckets.append((b1,cur_ways+1-b1))
	return mrc,buckets


def whirlpool_combine_ncurves_f(mrcs,nr_ways):
	mrcs=list(mrcs)
	nr_apps = len(mrcs)
	indexes = [0.0 for i in range(nr_apps)]
	mrc = []
	buckets = []
	for cur_ways in range(nr_ways):
		v_ns = access_indexes(mrcs,indexes)
		v_sum = sum(v_ns)
		for i in range(nr_apps):
			indexes[i] += v_ns[i]/v_sum
		mrc.append(sum(access_indexes_acc(mrcs,indexes)))
		# Save index to return the buckets in the same order as the mrcs
		total_ways = cur_ways +1
		cur_buckets = []
		total_weight = sum(indexes)
		for i in range(nr_apps):
			app_bucket=total_ways*indexes[i]/total_weight
			cur_buckets.append(app_bucket)
		buckets.append(cur_buckets)
	return (mrc,buckets)

def whirlpoolc_combine_ncurves_f(mrcs,llcmpkcs,nr_ways,agg_function=sum):
	nr_apps = len(mrcs)
	indexes = [0.0 for i in range(nr_apps)]
	mrc = []
	buckets = []
	for cur_ways in range(nr_ways):
		s_ns = access_indexes(llcmpkcs,indexes)
		s_sum = sum(s_ns)
		for i in range(nr_apps):
			indexes[i] += s_ns[i]/s_sum
		mrc.append(agg_function(access_indexes_acc(mrcs,indexes)))
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


def merge_clusters(cluster1,cluster2,curve,buckets):
	apps=cluster1.apps + cluster2.apps
	scaled_apps = []
	for i,app in enumerate(apps):
		assigned_ways = [frac_ways[i] for frac_ways in buckets]
		scaled_apps.append(App(app.name,app.build_scaled_properties(assigned_ways),parent=app))
	return Cluster(scaled_apps,curve,buckets)

def raw_unfairness(vals):
	if len(vals)==1:
		return vals[0]
	else:
		return max(vals)/min(vals)

# Class that stores a cluster of applications, the combined curve that represents it, and the buckets defining how the space is divided between apps
class Cluster:
	def __init__(self,apps,curve,buckets):
		self.apps = apps
		self.curve = curve
		self.buckets = buckets
		self.idx_cluster=min([app.bench_id for app in apps])

	def __repr__(self):
		names=[]
		for idx,app in enumerate(self.apps):
			names.append(app.name)
		return repr(names)

	def distance(self,cluster2,fix_partitioned=True):
		curves=[]
		for app in self.apps:
			curves.append(app.original_app.get_metric_table("llcmpki"))
		for app in cluster2.apps:
			curves.append(app.original_app.get_metric_table("llcmpki"))

		combined_curve,buckets = whirlpool_combine_ncurves_f([c.values for c in curves],len(self.curve))

		if fix_partitioned:
		
			# Critical fix: Normalize 1...n
			# It turns out that the cluster's curves are mere arrays
			# They need to be accessed in the range 1..n inside 
			# lookahead_algorithm_gen(), so we need to prepend
			# a dummy value as item 0.

			curve_a=list(self.curve)
			curve_b=list(cluster2.curve)
			curve_a.insert(0,len(self.curve))
			curve_b.insert(0,len(self.curve))
			pcurves=[curve_a,curve_b]
		else:
			pcurves=curves

		partitioned_curve=[]
		# Calculate partitioned curve
		for nr_ways in range(1,len(self.curve)+1):
			if nr_ways>1:
				ways_distr = lookahead_algorithm_gen(pcurves,nr_ways)
				partitioned_sum = 0
				for idxapp,assigned_ways in enumerate(ways_distr):
					partitioned_sum += pcurves[idxapp][assigned_ways]
			else:
					partitioned_sum=combined_curve[nr_ways-1]

			partitioned_curve.append(partitioned_sum)

		distance_sum = 0
		for i in range(len(self.curve)):
			distance_sum += abs(combined_curve[i] - partitioned_curve[i])

		return (distance_sum,combined_curve,buckets,partitioned_curve)	


	def distance_gen(self,cluster2,metric="llcmpki",metric_space="llcmpki",agg_function=sum,use_abs=True,fix_partitioned=True):
		curves=[]
		space_curves=[]
		for app in self.apps:
			curves.append(app.original_app.get_metric_table(metric))
			space_curves.append(app.original_app.get_metric_table(metric_space))

		for app in cluster2.apps:
			curves.append(app.original_app.get_metric_table(metric))
			space_curves.append(app.original_app.get_metric_table(metric_space))

		combined_curve,buckets = whirlpoolc_combine_ncurves_f([c.values for c in curves],[s.values for s in space_curves],len(self.curve),agg_function=agg_function)

		if fix_partitioned:
			
			# Critical fix: Normalize 1...n
			# It turns out that the cluster's curves are mere arrays
			# They need to be accessed in the range 1..n inside 
			# lookahead_algorithm_gen(), so we need to prepend
			# a dummy value as item 0.

			curve_a=list(self.curve)
			curve_b=list(cluster2.curve)
			curve_a.insert(0,len(self.curve))
			curve_b.insert(0,len(self.curve))
			pcurves=[curve_a,curve_b]
		else:
			pcurves=curves

		partitioned_curve=[]
		# Calculate partitioned curve
		for nr_ways in range(1,len(self.curve)+1):
			if nr_ways>1:
				ways_distr = lookahead_algorithm_gen(pcurves,nr_ways)
				ways_distr = list(map(int, ways_distr))
				values_partitioned=[]
				for idxapp,assigned_ways in enumerate(ways_distr):
					values_partitioned.append(pcurves[idxapp][assigned_ways])
				new_value=agg_function(values_partitioned)
			else:
				new_value=combined_curve[nr_ways-1]

			partitioned_curve.append(new_value)

		distance_curve=[]
		distance_sum = 0
		for i in range(len(self.curve)):
			if use_abs:
				distance_sum += abs(combined_curve[i] - partitioned_curve[i])
			else:
				distance_sum += combined_curve[i] - partitioned_curve[i]
			distance_curve.append(distance_sum)

		return (distance_sum,combined_curve,buckets,partitioned_curve,distance_curve)		


	def slowdown_distance(self,cluster2):
		return self.distance_gen(cluster2,metric="slowdown",metric_space="llcmpkc",agg_function=max,use_abs=False)

	def get_cluster_speedup_curve(self,nr_ways,total_apps,max_bandwidth):
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


def determine_best_partitioning(clusters,nr_ways,total_apps,max_bandwidth=float('Inf'),uoptions={}):
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

	if len(clusters) >  nr_ways:
		sp_aggregate = -1

	return (partitioning,per_app_partitioning,bw_app,sp_aggregate)

def determine_best_partitioning_slowdown(clusters,nr_ways,total_apps,max_bandwidth=float('Inf')):
	slowdown_curves=map(lambda x: x.curve, clusters)
	# Lookahead
	if len(clusters)==1:
		partitioning=[nr_ways]
	else:
		partitioning=lookahead_algorithm_gen(slowdown_curves,nr_ways)

	## Determine Per-application way-assignment
	per_app_partitioning=[]
	bw_aggregate=0
	slowdown_aggregate=0
	bw_app=[]
	slowdownv=[]

	for idxc,nr_ways_cluster in enumerate(partitioning):
		## Traverse applications in a cluster to determine 
		cl=clusters[idxc]
		buckets=cl.buckets
		for idxa,app in enumerate(cl.apps):
			# Now stores the cluster ways so it can later retrieve the scaled value from properties
			per_app_partitioning.append((app,nr_ways_cluster,idxc))
			#slowdown_aggregate+=
			slowdownv.append(app.get_metric_table("slowdown")[nr_ways_cluster])
			bw_cur=app.get_metric_table("bandwidth_mbps")[nr_ways_cluster]
			bw_aggregate+=bw_cur
			bw_app.append((bw_cur,app))

	if len(clusters) >  nr_ways:
		unfairness_value = 1000
	else:
		unfairness_value=max(slowdownv)/min(slowdownv)

	return (partitioning,per_app_partitioning,bw_app,unfairness_value)


def get_kpart_schedule(workload,nr_ways,max_bandwidth=float('Inf'),uoptions={}):
	curClusters = []
	total_apps = len(workload)

	for i,app in enumerate(workload):
		buckets = [[nw] for nw in range(1,nr_ways + 1)]
		curve = app.get_metric_table("llcmpki").values
		app.bench_id = i
		cl=Cluster([app],curve,buckets)
		curClusters.append(cl)

	if uoptions["best_metric"] == "speedup":
		solutions = [(list(curClusters),determine_best_partitioning(curClusters,nr_ways,total_apps,max_bandwidth,uoptions=uoptions))]
		# (partitioning,per_app_partitioning,bw_app,sp_aggregate) <- determine_best_part()
	else:
		solutions = [(list(curClusters), determine_best_partitioning_slowdown(curClusters, nr_ways, total_apps, max_bandwidth))]

	while len(curClusters) > 1:

		cluster_data=None
		min_found=30000000
		min_idx=(-1,-1)

		## Compute distance and keep track of min
		for idx,clusterRef in enumerate(curClusters):
			for i in range(idx+1,len(curClusters)):
				cluster=curClusters[i]
				#distance,combined_curve,buckets,partitioned_curve)=clusterRef.distance_gen(cluster,metric=uoptions["dist_metric"],metric_space=uoptions["buckets_metric"]
				#																			,agg_function=uoptions["agg_function"],use_abs=uoptions["use_abs"])
				(distance, combined_curve, buckets, partitioned_curve) = clusterRef.distance(cluster)
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
		if uoptions["best_metric"] == "speedup":
			solutions.append((list(curClusters), determine_best_partitioning(curClusters,nr_ways,total_apps,max_bandwidth)))
		else:
			solutions.append((list(curClusters), determine_best_partitioning_slowdown(curClusters, nr_ways, total_apps, max_bandwidth)))


	max_sp=-1
	best_k=0

	for (idx,sol) in enumerate(solutions):
		(clusters,data)=sol
		(partitioning,per_app_partitioning,bw_app,sp_aggregate)=data

		if sp_aggregate > max_sp:
			max_sp=sp_aggregate
			best_k=idx

	#assert sum(solutions[best_k][1][0])==nr_ways, "The algorithm is assigning all available ways to clusters"

	return (solutions,best_k)

def get_kpart_schedule_optimized(workload,nr_ways,max_bandwidth=float('Inf'),uoptions={}):
	curClusters = []
	total_apps = len(workload)	
	
	##Create distance matrix
	distance_matrix=np.full((total_apps,total_apps),float('Inf'))

	## Create lists for curves - triplets (combined_curve,buckets,partitioned_curve)
	kpart_curves=[]
	for i in range(total_apps):
		row=[]
		for j in range(total_apps):
			row.append(None)
		kpart_curves.append(row)


	for i,app in enumerate(workload):
		buckets = [[nw] for nw in range(1,nr_ways + 1)]
		curve = app.get_metric_table("slowdown").values
		## Supercritical
		app.bench_id = i
		cl=Cluster([app],curve,buckets)
		curClusters.append(cl)

	##Calculate initial distances and curves
	for i,clusterRef in enumerate(curClusters):
		idx_i=clusterRef.idx_cluster
		for j in range(i+1,len(curClusters)):
			cluster=curClusters[j]
			idx_j=cluster.idx_cluster
			(distance,combined_curve,buckets,partitioned_curve,distance_curve)=clusterRef.distance_gen(cluster,metric="slowdown",metric_space="llcmpkc",agg_function=raw_unfairness,use_abs=False)
			##Store in matrices
			distance_matrix[idx_i,idx_j]=distance
			kpart_curves[idx_i][idx_j]=(combined_curve,buckets,partitioned_curve,clusterRef,cluster)

	solutions = [(list(curClusters),determine_best_partitioning_slowdown(curClusters,nr_ways,total_apps,max_bandwidth))]
	# (partitioning,per_app_partitioning,bw_app,sp_aggregate) <- determine_best_part()

	it=1
	while len(curClusters) > 1:

		## Find min distance
		min_i,min_j=np.unravel_index(distance_matrix.argmin(), distance_matrix.shape)
		distance=distance_matrix[min_i,min_j]
		(combined_curve,buckets,partitioned_curve,clusterA,clusterB)=kpart_curves[min_i][min_j]

		##Determine location of prev_clusters (traverse all, it could be more efficient...)
		min_idx=[0,0]
		for idx,cluster in enumerate(curClusters):
			if cluster.idx_cluster == clusterA.idx_cluster:
				min_idx[0]=idx
			elif cluster.idx_cluster == clusterB.idx_cluster:
				min_idx[1]=idx

		print("***Iteration %d***" % it)
		print(distance_matrix)
		print(distance)
		print("Selection:",clusterA,clusterB)
		print("***********************************")
		it=it+1

		## Merge 2 closest clusters
		new_cluster=merge_clusters(clusterA,clusterB,combined_curve,buckets)

		## Sorted is critical!
		min_idx=sorted(min_idx)

		## Remove prev clusters...
		del curClusters[min_idx[0]]
		del curClusters[min_idx[1]-1]

		##One of the indexes is going away (the minimum one). Update the distance matrix to discard values in the future
		removed_index=clusterA.idx_cluster if clusterA.idx_cluster>clusterB.idx_cluster else clusterB.idx_cluster
		distance_matrix[removed_index,:]=float('Inf')
		distance_matrix[:,removed_index]=float('Inf')

		#Add merged cluster ...
		curClusters.insert(0,new_cluster)
		solutions.append((list(curClusters),determine_best_partitioning_slowdown(curClusters,nr_ways,total_apps,max_bandwidth)))	

		## Calculate new distance if this is not the last iteration
		if len(curClusters) > 1:
			## Update distance info
			##Calculate initial distances and curves
			idx_new=new_cluster.idx_cluster

			for cluster in curClusters:
				if cluster.idx_cluster!=idx_new:
					idx_other=cluster.idx_cluster
					if (idx_new<idx_other):
						(idx_i,idx_j)=(idx_new,idx_other)
						(cluster_i,cluster_j)=(new_cluster,cluster)
					else:
						(idx_i,idx_j)=(idx_other,idx_new)
						(cluster_i,cluster_j)=(cluster,new_cluster)
		
					(distance,combined_curve,buckets,partitioned_curve,distance_curve)=cluster_i.distance_gen(cluster_j,metric="slowdown",metric_space="llcmpkc",agg_function=raw_unfairness,use_abs=False)
		
					##Store in matrices
					distance_matrix[idx_i,idx_j]=distance
					kpart_curves[idx_i][idx_j]=(combined_curve,buckets,partitioned_curve,cluster_i,cluster_j)


	## Determine BEST K
	use_min=True
	if use_min:
		best_score=-float('Inf')
	else:
		best_score=float('Inf')
	best_k=0

	for (idx,sol) in enumerate(solutions):
		(clusters,data)=sol
		(partitioning,per_app_partitioning,bw_app,score)=data

		if use_min and score < best_score:
			best_score=score
			best_k=idx
		if not use_min and score > best_score:
			best_score=score
			best_k=idx

	#assert sum(solutions[best_k][1][0])==nr_ways, "The algorithm is assigning all available ways to clusters"

	return (solutions,best_k)



def get_kpart_best_gen(workload,nr_ways,max_bandwidth=float('Inf'),debugging=False,variant=0,**kwargs):
	kwargs.setdefault("dist_metric", "llcmpki")  # Curve distance /slowdown
	kwargs.setdefault("buckets_metric", "llcmpki")  # Bucket's space distribution /llcmpkc
	#kwargs.setdefault("ucp_metric", "llcmpki") # Metric to decide partitions' sizes /llcmpkc
	kwargs.setdefault("use_abs", True) # Whether to use absolute value in distance calcs
	kwargs.setdefault("agg_function", sum) # sum/max (function to use when the curves
	kwargs.setdefault("best_metric", "speedup") # sp/stp/unf

	if "user_options" in kwargs:
		## Merge dict with user options... [Overriding defaults...]
		kwargs.update(kwargs["user_options"])

	## Invoke kpart.....
	if variant==0:
		(solutions,k)=get_kpart_schedule(workload,nr_ways,max_bandwidth,uoptions=kwargs)
	else:
		(solutions,k)=get_kpart_schedule_optimized(workload,nr_ways,max_bandwidth,uoptions=kwargs)

	if debugging:
		for sol in solutions:
			clusters,(partitioning,per_app_partitioning,bw_app,sp_aggregate) = sol
			bw_app = [(bw,app.name) for (bw,app) in bw_app]
			print (clusters, sp_aggregate)
			#(workload_alt,solution_alt)=zip(*per_app_partitioning)
			#print sp_aggregate, unfairness(list(workload_alt),list(solution_alt),max_bandwidth), clusters, partitioning, list(solution_alt)
		print("Best k =",k)

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
	cluster = list(cluster)
	if len(cluster)==1:
		return cluster
	## Get MRC for each APP
	curves=map(lambda x: x.get_metric_table(metric).values,cluster)
	## Apply whirlpool
	comb_mrc,buckets = whirlpool_combine_ncurves_f(curves,max_ways)
	hacked_apps=[]
	for i,app in enumerate(cluster):
		## Potential optimization in reducing the number of metrics
		props=app.build_scaled_properties([b_ways[i] for b_ways in buckets], []) #["ipc","llcmpki","llcmpkc","slowdown"])
		hacked_apps.append(App(app.name,props,app))

	return hacked_apps


# Function that scales down the metric tables of each cluster in an array to a specific number of ways
# Params: Array of clusters and number of ways
# Return: Scaled clusters and apps
def get_patched_clustering_and_workload(clustering,nr_ways):
	patched_apps=[]
	patched_clustering=[]	
	
	for cluster in clustering:
		patched_cluster=get_scaled_properties_cluster(cluster,nr_ways)
		patched_apps.extend(patched_cluster)
		patched_clustering.append(patched_cluster)

	return (patched_clustering,patched_apps)


# Function that normalizes the output for a given clustering solution
def normalize_output_for_clustering_solution(workload,clusters,clustering_sol,nr_ways,llc_mapping=False):
	for i,app in enumerate(workload):
		app.bench_id = i
		
	## Calculate patched apps and so on
	(patched_clustering,patched_apps)=get_patched_clustering_and_workload(clusters,nr_ways)

	if 	llc_mapping:
		##Set cache ids... based on clusters
		cache_id=0
		per_cluster_masks=[]
		for clust in patched_clustering:
			per_cluster_masks.append(hex((1<<nr_ways)-1))
			for app in clust:
				app.llc_id=cache_id
			cache_id+=1

			## Each llc id can be considered its own mask...

	else:
		## First determine the LLC ID associated with each cluster
		## and treat each thing separately
		llc_cluster_ids=[]
		for clust in patched_clustering:
			cluster_llcid=clust[0].llc_id
			if cluster_llcid==-1:
				cluster_llcid=0 ## Normalize
			llc_cluster_ids.append(cluster_llcid)

		## Build cool stuff
		clustering_partial_solutions=[[] for i in range(len(llc_cluster_ids))]
		for i,clust in enumerate(patched_clustering):
			cluster_llcid=llc_cluster_ids[i]
			clustering_partial_solutions[cluster_llcid].append(clustering_sol[i])

		per_cluster_masks=[]
		for psol in clustering_partial_solutions:
			per_cluster_masks.extend(get_partition_masks(psol))

	per_app_masks = [None] * len(patched_apps)
	per_app_ways = [None] * len(patched_apps)
	cluster_ids = [None] * len(patched_apps)
	patched_workload=[None]*len(patched_apps)

	## Build per app ways and per app masks
	for i,cluster in enumerate(patched_clustering):
		cluster_ways=clustering_sol[i]
		cluster_mask=per_cluster_masks[i]
		for app in cluster:
			orig_i = app.original_app.bench_id
			per_app_ways[orig_i] = cluster_ways
			per_app_masks[orig_i] = cluster_mask
			cluster_ids[orig_i] = i

	for app in patched_apps:
		orig_i = app.original_app.bench_id
		#print orig_i,app
		patched_workload[orig_i] = app

	return (patched_workload,per_app_ways,per_app_masks,cluster_ids)


def write_rewind(fd, val):
	os.write(fd,val)
	os.lseek(fd,0,os.SEEK_SET)

## Read pmc config to find nr_cos and nr_ways
def read_pmc_config():
	nr_cos = 0
	nr_ways = 0
	fd = open("/proc/pmc/config")
	config = fd.read()
	for line in config.split("\n"):
		if "cat_nr_cos_available" in line:
			nr_cos = int(line.split("cat_nr_cos_available")[1][1:])
		if "llc_cbm0" in line:
			all_ways_mask = line.split("=")[1]
			nr_ways = len(bin(int(all_ways_mask,16)))-2
	fd.close()
	return (nr_cos,nr_ways)

def get_user_assignment(workload,nr_ways,workload_name,**kwargs):
	if "user_options" in kwargs:
		## Merge dict with user options... [Overriding defaults...]
		kwargs.update(kwargs["user_options"])
	workloads_assignment = kwargs["user_assignment"]
	w_id = int(workload_name.strip("W"))-1
	if w_id > len(workloads_assignment):
		print("Your user-sched file has not enough workloads specified w.r.t. your workload.csv file")
		exit(1)
	(per_app_ways,per_app_masks,cluster_id) = workloads_assignment[w_id]

	patched_workload=[None]*len(per_app_ways)
	for cl_ix in range(len(cluster_id)):
		# Get the corresponding indexes of this cluster to get the specific apps in the workload
		cl_apps_ixs = [i for i,clust_i in enumerate(cluster_id) if clust_i==cl_ix]
		# Obtain the apps in the workload of this cluster
		clust_apps = [workload[app_ix] for app_ix in cl_apps_ixs]
		# Re build scaled apps to return the result and be able to obtain the performance metrics
		scaled_apps = get_scaled_properties_cluster(clust_apps, nr_ways)
		# Rebuild patched workload
		for i in range(len(cl_apps_ixs)):
			patched_workload[cl_apps_ixs[i]] = scaled_apps[i]

	return (patched_workload,per_app_ways,per_app_masks,cluster_id)


if __name__ == "__main__":
#	l=[c for c in generate_possible_clusters(["A","B","C","D","E"],3)]
#	print l
#	print len(l)
#	exit(0)
	#bb_generate_solutions_to_explore(20,6,28)
	subnodes=bb_break_into_subnodes([2],4,20,5)
	print(subnodes)
	#write_dot_tree(6,4)
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
