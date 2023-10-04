#!/usr/bin/env python
# -*- coding: utf-8 -*-
from simulator_core import *
import csv

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

		## Reverse order for next iteration
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
	d_min = 1.0 # min distance between points of diff clusters
	d_max = 0.0 # max within-cluster distance
	for i in range(k):
		intra_cluster_distance = abs(clusters[i][1] - clusters[i][0])
		if (intra_cluster_distance > d_max):
			d_max = intra_cluster_distance
		if (i<k-1):
			inter_cluster_distance = abs(clusters[i+1][0] - clusters[i][1])
			if inter_cluster_distance < d_min:
				d_min = inter_cluster_distance
	#assert d_max!=0.0,"d_max is 0, it will causing division by 0"
	if d_max==0.0:
		return 300000.0
	return d_min / d_max


def estimate_ways_exponential(x, max_ways):
	if max_ways == 20:
		nways = int(np.round(19.26833327*np.exp(0.67431637*x) - 17.72167744))
	elif max_ways == 11:
		# parameters for 10 ways scaling down the paper ways
		nways = int(np.round(4.09659516*np.exp(1.19540837*x) - 2.48157709))
	return nways

def get_schedule_dunn(workload, nr_ways):
	n_apps = len(workload)
	hacked_apps=get_scaled_properties_cluster(workload,nr_ways)
	stalls = []

	if nr_ways>=len(workload): 
		for i,app in enumerate(workload):
			stalls.append(app.properties["stalls_l2_miss"][1])
	else:
		for i,app in enumerate(hacked_apps):
			stalls.append(app.properties["stalls_l2_miss"][nr_ways])

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
	apps_ix = []
	cluster_id = []
	for i_app in range(len(best_partitions)):
		# Assign CLOS to PID
		per_app_masks.append(masks[best_partitions[i_app]])
		per_app_ways.append(centroids_ways[best_partitions[i_app]])
		cluster_id.append(best_partitions[i_app])

	# (workload,per_app_ways,per_app_masks,cluster_id)
	return (per_app_ways,per_app_masks,cluster_id)



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
# Calculate LLC partitioning based on rate of demand (BW)
def on_demand_partitioning(workload,nr_ways,max_bandwidth=float('Inf')):
	bw_alone_vector=list(map(lambda app:app.properties["bandwidth_mbps"][nr_ways],workload))

	## Reserve one way for each benchmark
	nr_apps=len(workload)
	ways_to_share=nr_ways-nr_apps

	## Apply BW model if necessary
	if max_bandwidth != float('Inf'):	
		## Replace slowdown vector ...
		total_bw,bw_shared=determine_bw_shared(bw_alone_vector,max_bandwidth)
	else:
		total_bw=sum(bw_alone_vector)
		bw_shared=bw_alone_vector

	solution=list(map(lambda x: int(ways_to_share*(x/total_bw))+1,bw_shared)) ## Plus one -> at least one each
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


## Return category from traces
## 0->light_sharing
## 1->streaming
## 2-> cache-sensitive
def get_app_category(app,nr_ways):
	df=app.properties
	memory_intensive_misses=df[(df.slowdown<=1.03) & (df.llcmpkc>=14)]
	benchmark_slowdown_alto=df[df.slowdown>=1.06]
	sensible=df[(df.slowdown>=1.05) & (df.index>=2)]
	is_mi=len(memory_intensive_misses)>0
	slowdown_alto=len(benchmark_slowdown_alto)>0

	## Streaming
	if is_mi and not slowdown_alto:
		return 1
	elif len(sensible)>0:
		return 2
	else:
		return 0

## Return category from traces
## 0->light_sharing
## 1->streaming
## 2-> sensitive
## 3-> highly_sensitive
def get_app_category2(app,nr_ways):
	df=app.properties
	memory_intensive_misses=df[(df.slowdown<=1.03) &  (df.llcmpkc>=14)]
	benchmark_slowdown_alto=df[df.slowdown>=1.06]
	sensible=df[(df.slowdown>=1.05) & (df.index>=2)]
	highly_sensitive=df[(df.slowdown>=1.15) & (df.index>=2)]
	is_mi=len(memory_intensive_misses)>0
	slowdown_alto=len(benchmark_slowdown_alto)>0

	## Streaming
	if is_mi and not slowdown_alto:
		return 1
	elif len(highly_sensitive)>0:
		return 2
	elif len(sensible)>0:
		return 3
	else:
		return 0


## Variants....
## 0 -> balance among all
## 1 -> balance among streaming only
## 2 -> balance among cache sensitive only
def classification_clustering(workload,nr_ways,max_bandwidth,variant=0):
	classification=list(map(lambda bench: get_app_category(bench,nr_ways),workload))
	nr_clusters=0
	light_sharing_apps=[]
	nr_light_sharing=0
	clusters=[]
	pos_streaming=[]
	pos_sensitive=[]
	apps_assigned=[]

	## As many clusters as the number of apps that are not light sharing
	for i,app in enumerate(workload):
		app.bench_id=i ## for sorting
		app_class=classification[i]
		if app_class==0:
			light_sharing_apps.append(app)
			nr_light_sharing=nr_light_sharing+1
		else:
			clusters.append([app])
			apps_assigned.append(app)
			if app_class==1:
				pos_streaming.append(nr_clusters)
			else:
				pos_sensitive.append(nr_clusters)
			nr_clusters=nr_clusters+1


	## Si excece numero vias assert
	## TODO for now
	assert(nr_clusters<=nr_ways)

	if variant==0:
		target_cluster_ids=pos_streaming+pos_sensitive
	elif variant==1:
		target_cluster_ids=pos_streaming
	else:
		target_cluster_ids=pos_sensitive

	nr_target_clusters=len(target_cluster_ids)

	if nr_target_clusters==0:
		nr_reserved_ways=1
	else:
		nr_reserved_ways=0		

	clustering_sol=get_schedule_UCP_gen(apps_assigned,nr_ways-nr_reserved_ways,metric="slowdown")

	## Separate partition for light
	if nr_reserved_ways>0:
		clusters.append(light_sharing_apps)
		clustering_sol.append(nr_reserved_ways)
	else:
		## Even out streaming clusters in rotatory assignments
		idx=0
		nr_remaining=nr_light_sharing

		while nr_remaining>0:
			app=light_sharing_apps.pop(0)
			id_cluster=target_cluster_ids[idx%nr_target_clusters]
			clusters[id_cluster].append(app)
			nr_remaining=nr_remaining-1
			idx=idx+1

	return normalize_output_for_clustering_solution(workload,clusters,clustering_sol,nr_ways)


def benefit_way_stealing(workload,merged,idx_app,ucp_partitioning):
	slowdown_red=np.full(len(workload),float('Inf'))
	this_app=workload[idx_app]
	app_assigned_ways=ucp_partitioning[idx_app]
	my_slowdown_curve=this_app.properties.slowdown
	cur_slowdown=this_app.properties.slowdown[app_assigned_ways]

	for i,app in enumerate(workload):
		if merged[i] or i==idx_app:
			slowdown_red[i]=float('Inf')
		else: ## Calculate local distance
			## What if it was possible to transfer one way
			that_app=workload[i]
			those_ways=ucp_partitioning[i]
			slowdown_curve=that_app.properties.slowdown

			if those_ways==1:
				slowdown_red[i]=float('Inf')
			else:
				slowdown_that=slowdown_curve[those_ways-1]

				if slowdown_that>=cur_slowdown:
					slowdown_red[i]=float('Inf')
				else:
					slowdown_scenario_a=my_slowdown_curve[app_assigned_ways+1]-my_slowdown_curve[app_assigned_ways]
					slowdown_scenario_b=slowdown_curve[those_ways]-slowdown_curve[those_ways-1]
					slowdown_red[i]=slowdown_scenario_a-slowdown_scenario_b

	return slowdown_red

##
# Public function that calculates a solution for workload passed as a parameter using the lookahead algorithm.
# Params: workload (array of Apps) and slots available.
# Return: cache partitioning for the workload
def  ucp_unfairness(workload, nr_ways, uoptions={}):
	uoptions.setdefault("min_slowdown",True)
	curves=[]
	nr_ways_equal=int(round(nr_ways/len(workload)))
	for app in workload:
		curves.append(app.properties.slowdown)
	ucp_solution=lookahead_algorithm_gen(curves,nr_ways)

	if uoptions["min_slowdown"]:
		slowdown_curves=[app.properties.slowdown for app in workload]
		slowdown_vector=np.array([app.properties.slowdown[ucp_solution[i]] for i,app in enumerate(workload)])
		nr_apps=len(workload)
		merged=[False for _ in workload]

		for i in range(nr_apps):
			idx_max=slowdown_vector.argmax()

			if merged[idx_max]:
				break

			slowdown_reductions=benefit_way_stealing(workload,merged,idx_max,ucp_solution)


			## Find min ex
			idx_min=slowdown_reductions.argmin()
			benefit=slowdown_reductions[idx_min]

			if benefit>0.07:
				break

			## Transfer way
			ucp_solution[idx_max]=ucp_solution[idx_max]+1
			ucp_solution[idx_min]=ucp_solution[idx_min]-1
			merged[idx_max]=True

			## Update slowdown 
			slowdown_vector[idx_max]=slowdown_curves[idx_max][ucp_solution[idx_max]]
			slowdown_vector[idx_min]=slowdown_curves[idx_min][ucp_solution[idx_min]]



	return ucp_solution


def get_benchmark_categories_from_file(catfile, label_format=False):
	csv_file = open(catfile)
	csv_entry = csv.reader(csv_file)

	classes={}

	for entry in csv_entry:
		classes[entry[0]]=int(entry[1])

	csv_file.close()

	if label_format:
		categories_labels = {0:"light-sharing", 1:"streaming", 2:"cache-sensitive", 3:"unknown"}
		for key in classes.keys():
			if key not in classes:
				print("Warning key not found (%d). Defaulting to light-sharing", key)
				classes[key]=0
			else:
				classes[key] = categories_labels[classes[key]]

	return classes


def lfoc(workload,nr_ways,max_bandwidth,uoptions={}):
	from opt_clustering import get_optimal_clustering_bf
	uoptions.setdefault("max_streaming",5)
	uoptions.setdefault("light_per_streaming",2) ## Number of light apps for each missing streaming in streaming reserved partitions
	uoptions.setdefault("max_ways_streaming",2)
	uoptions.setdefault("streaming_part_size",1)	
	uoptions.setdefault("use_pair_clustering",False)
	uoptions.setdefault("collide_streaming_partitions",True)
	uoptions.setdefault("opt",False)
	uoptions.setdefault("simple_output",False)

	ucp_metric="slowdown"
	if "ucp_metric" in uoptions:
		ucp_metric=uoptions["ucp_metric"]

	if "benchmark_categories" in uoptions:
		class_dict=get_benchmark_categories_from_file(uoptions["benchmark_categories"])
		classification=[class_dict[bench.name] for bench in workload]
	else:
		classification=list(map(lambda bench: get_app_category(bench,nr_ways),workload))

	streaming_part_size=uoptions["streaming_part_size"]

	nr_clusters=0
	light_sharing_apps=[]
	nr_light_sharing=0
	nr_streaming=0
	nr_sensitive=0
	streaming_apps=[]
	clusters=[]
	pos_sensitive=[]
	pos_streaming=[]
	apps_assigned=[]

	## As many clusters as the number of sensitive benchmarks
	for i,app in enumerate(workload):
		if not uoptions["simple_output"]:
			app.bench_id=i ## for sorting
		app_class=classification[i]
		if app_class>=2:
			#clusters.append([app])
			apps_assigned.append(app)
			#pos_sensitive.append(nr_clusters)
			#nr_clusters=nr_clusters+1
			#nr_sensitive=nr_sensitive+1
		else:
			if app_class==0:
				light_sharing_apps.append(app)
				nr_light_sharing=nr_light_sharing+1
			else:
				streaming_apps.append(app)
				nr_streaming=nr_streaming+1


	## Calculate the number of ways for streaming apps
	if nr_streaming>0:
		streaming_per_part=uoptions["max_streaming"]
		nr_reserved_ways=((nr_streaming+streaming_per_part-1)*streaming_part_size)//streaming_per_part  ## Round up

		if nr_reserved_ways>uoptions["max_ways_streaming"]:
			nr_reserved_ways=uoptions["max_ways_streaming"]
			# Raise cap
			streaming_per_part=(nr_streaming+nr_reserved_ways-1)//(nr_reserved_ways/streaming_part_size)
			## Load factor depends on how many streamings 
	else:
		nr_reserved_ways=0

	## Corner case, no sensitive applications
	if len(apps_assigned)==0: 
		clustering_sol=[nr_ways]
		one_cluster=[app for app in workload]
		clusters=[one_cluster]
		if uoptions["simple_output"]:
			return (clusters,clustering_sol)
		else:
			return normalize_output_for_clustering_solution(workload,clusters,clustering_sol,nr_ways)


	if uoptions["use_pair_clustering"]:
		## Invoke pair_clustering for sensitive
		(solutions,k)=pair_clustering_core2(apps_assigned,nr_ways-nr_reserved_ways,max_bandwidth,False,uoptions=uoptions)

		best_sol=solutions[k]
		sol_clusters=best_sol[0]
		clustering_sol=best_sol[1][0]

		## Update LFOC variables
		## nr_sensitive now indicates the number of clusters 
		## with sensitive apps
		nr_clusters=nr_sensitive=len(sol_clusters)
		
		## Turn each app into "original app"
		clusters=[list(map(lambda app:app.original_app,cluster.apps)) for cluster in sol_clusters]

		pos_sensitive=[i for i in range(nr_clusters)]

	elif uoptions["opt"]:
		(patched_workload, per_app_ways, per_app_masks, cluster_id, total_branches, metadata) = get_optimal_clustering_bf(apps_assigned, max_slowdown_unfairness, False, nr_ways-nr_reserved_ways, max_bandwidth, multiprocessing=True,	user_options=uoptions)
	
		nr_clusters=nr_sensitive=max(cluster_id)+1
		pos_sensitive=[i for i in range(nr_clusters)]
		clusters = [[] for _ in range(nr_clusters)]
		clustering_sol = [[] for _ in range(nr_clusters)]


		for idxc in range(nr_clusters):
			# Search for first occurrence of idxc
			cluster_i = cluster_id.index(idxc)
			clustering_sol[idxc] = per_app_ways[cluster_i]
		
		for i,app in enumerate(patched_workload):
			clusters[cluster_id[i]].append(app.original_app)


		
	else:
		## Add sensitive clusters right here
		for app in apps_assigned:
			clusters.append([app])
			pos_sensitive.append(nr_clusters)
			nr_clusters=nr_clusters+1
			nr_sensitive=nr_sensitive+1

		assert(nr_clusters<=nr_ways)

		## Run UCP only with sensitive benchmarks
		clustering_sol=get_schedule_UCP_gen(apps_assigned,nr_ways-nr_reserved_ways,metric=ucp_metric)

	## Create partitions for streaming and assign streaming
	if nr_reserved_ways:
		if uoptions["collide_streaming_partitions"]:
			## Create single cluster with all streaming apps
			streaming_clust=[]
			for st in  streaming_apps:
				streaming_clust.append(st)
			clustering_sol.append(nr_reserved_ways)
			clusters.append(streaming_clust)
			pos_streaming.append(nr_clusters)
		else:
			stream_idx=0
			remaining_streaming=nr_streaming

			## For each free range
			for i in range(nr_reserved_ways):

				streaming_clust=[]
				stream_to_assign=min(streaming_per_part,remaining_streaming)

				for j in range(0,stream_to_assign):
					streaming_clust.append(streaming_apps[stream_idx])
					stream_idx=stream_idx+1
					remaining_streaming=remaining_streaming-1

				# Create streaming cluster
				clustering_sol.append(streaming_part_size)
				clusters.append(streaming_clust)
				pos_streaming.append(nr_clusters+i)



	## Assign lights  (load_balance())
	## Traverse streaming first ... [Require sorted stuff]

	if nr_light_sharing>0:

		if len(pos_streaming)>0: 

			idx_stream=0
			idx_light=0
			while  nr_light_sharing>0 and idx_stream<len(pos_streaming):
				streaming_clust=clusters[pos_streaming[idx_stream]]
				## Determine how many we can actually fit in here
				nr_stream_in_part=len(streaming_clust)
				nr_ways_streaming_part=clustering_sol[-1]
				room=uoptions["max_streaming"]*nr_ways_streaming_part-nr_stream_in_part

				if room>0:
					light_to_assign=min(room*uoptions["light_per_streaming"],nr_light_sharing)
					
					for i in range(0,light_to_assign):
						app=light_sharing_apps.pop(0)
						streaming_clust.append(app)
						nr_light_sharing=nr_light_sharing-1

				idx_stream=idx_stream+1

		## Go ahead and assign in RR the LIGHT SHARING apps among rest
		idx=0
		nr_remaining=nr_light_sharing

		while nr_light_sharing>0:
			app=light_sharing_apps.pop(0)
			id_cluster=pos_sensitive[idx%nr_sensitive]
			clusters[id_cluster].append(app)
			nr_light_sharing=nr_light_sharing-1
			idx=idx+1

	if uoptions["simple_output"]:
		return (clusters,clustering_sol)
	else:
		return normalize_output_for_clustering_solution(workload,clusters,clustering_sol,nr_ways)

##
# Removes light applications before applying clustering optimal algorithm and then assigns them
# to the resulting clusters giving priority to clusters with streaming apps as focc
def clustering_pseudo_optimal(algorithm,original_workload,nr_ways,max_bandwidth=float('Inf'), multiprocessing=False,user_options={},opt_spec=(unfairness,"False")):
	from opt_clustering import get_optimal_clustering_bf
	user_options.setdefault("max_streaming", 5)
	user_options.setdefault("light_per_streaming",
						2)  ## Number of light apps for each missing streaming in streaming reserved partitions

	app_class = {}
	light_sharing_apps = []
	workload = []
	# classes: 0=light, 1=streaming, 2=sensitive
	for i,app in enumerate(original_workload):
		# Be careful and DO NOT use bench_id, it is changed inside get_optimal_clustering
		app.original_app.ix = i
		app_class[app.name] = get_app_category(app, nr_ways)
		if app_class[app.name] == 0: # light-sharing benchmarks
			light_sharing_apps.append(app)
		else:
			workload.append(app)

	## Too many MI-apps lead to single part all ...
	if len(workload)>11:
		## all apps in the same cluster
		clusters=[original_workload]
		cluster_way={0: nr_ways}
		cluster_mask={0: hex((1<<nr_ways) -1)}
		total_branches=0
		metadata={}
	else:
		if algorithm=="poptc-unf":
			(patched_workload, per_app_ways, per_app_masks, cluster_id, total_branches, metadata) = get_optimal_clustering_bf(workload, unfairness_max_throughput, False, nr_ways, max_bandwidth, multiprocessing=multiprocessing,	user_options=user_options)
		elif algorithm=="poptc-stp":
			(patched_workload, per_app_ways, per_app_masks, cluster_id, total_branches, metadata) = get_optimal_clustering_bf(workload, throughput, True, nr_ways, max_bandwidth, multiprocessing, user_options=user_options)
		elif algorithm=="poptc-cov":
			(patched_workload, per_app_ways, per_app_masks, cluster_id, total_branches, metadata) = get_optimal_clustering_bf(workload, max_slowdown_unfairness, False, nr_ways, max_bandwidth, multiprocessing=multiprocessing,	user_options=user_options)
		else:
			(patched_workload, per_app_ways, per_app_masks, cluster_id, total_branches, metadata) = get_optimal_clustering_bf(workload, opt_spec[0], opt_spec[1], nr_ways, max_bandwidth, multiprocessing=multiprocessing,	user_options=user_options)			 

		nr_clusters = max(cluster_id)+1
		pos_streaming = []
		pos_sensitive = []
		clusters = [[] for _ in range(nr_clusters)]
		nr_light_sharing = len(light_sharing_apps)
		cluster_way = {}
		cluster_mask = {}


		for idxc in range(nr_clusters):
			cluster_i = cluster_id.index(idxc)
			cluster_way[idxc] = per_app_ways[cluster_i]
			cluster_mask[idxc] = per_app_masks[cluster_i]

		for i,app in enumerate(patched_workload):
			clusters[cluster_id[i]].append(app.original_app)
			if app_class[app.name] == 1:
				pos_streaming.append(cluster_id[i])
			elif app_class[app.name] == 2:
				pos_sensitive.append(cluster_id[i])

		nr_sensitive = len(pos_sensitive)

		## Assign lights  (load_balance())
		## Traverse streaming first ... [Require sorted stuff]

		if nr_light_sharing > 0:

			if len(pos_streaming) > 0:
				#print pos_streaming
				#print clusters
				idx_stream = 0
				while nr_light_sharing > 0 and idx_stream < len(pos_streaming):
					#print idx_stream, pos_streaming[idx_stream]
					streaming_clust = clusters[pos_streaming[idx_stream]]
					## Determine how many we can actually fit in here
					nr_stream_in_part = len(streaming_clust)
					room = user_options["max_streaming"] - nr_stream_in_part

					if room > 0:
						light_to_assign = min(room * user_options["light_per_streaming"], nr_light_sharing)

						for i in range(0, light_to_assign):
							app = light_sharing_apps.pop(0)
							streaming_clust.append(app)
							nr_light_sharing = nr_light_sharing - 1

					idx_stream = idx_stream + 1

			## Go ahead and assign in RR the LIGHT SHARING apps among rest
			idx = 0
			while nr_light_sharing > 0:
				app = light_sharing_apps.pop(0)
				id_cluster = pos_sensitive[idx % nr_sensitive]
				clusters[id_cluster].append(app)
				nr_light_sharing = nr_light_sharing - 1
				idx = idx + 1

	# Re-build solution with additional light apps
	per_app_masks = [None] * len(original_workload)
	per_app_ways = [None] * len(original_workload)
	cluster_id = [None] * len(original_workload)
	patched_workload = [None] * len(original_workload)

	for idxc,cluster in enumerate(clusters):
		patched_cluster = get_scaled_properties_cluster(cluster,nr_ways)
		for app in patched_cluster:
			orig_i = app.original_app.ix
			per_app_ways[orig_i] = cluster_way[idxc]
			per_app_masks[orig_i] = cluster_mask[idxc]
			patched_workload[orig_i] = app
			cluster_id[orig_i] = idxc

	return (patched_workload, per_app_ways, per_app_masks, cluster_id, total_branches, metadata)


##
# Removes light applications before applying clustering optimal algorithm and then assigns them
# to the resulting clusters giving priority to clusters with streaming apps as focc
def mapping_pseudo_optimal(algorithm,original_workload,nr_ways,max_bandwidth=float('Inf'), multiprocessing=False,user_options={},opt_spec=(unfairness,"False")):
	from opt_clustering import get_optimal_clustering_bf
	user_options.setdefault("max_streaming", 5)
	user_options.setdefault("light_per_streaming",
						2)  ## Number of light apps for each missing streaming in streaming reserved partitions
	user_options.setdefault("cores_per_llc",4)

	if "benchmark_categories" in user_options:
		class_dict=get_benchmark_categories_from_file(user_options["benchmark_categories"])
		classification=[class_dict[bench.name] for bench in original_workload]
	else:
		classification=list(map(lambda bench: get_app_category(bench,nr_ways),original_workload))


	app_class = {}
	light_sharing_apps = []
	workload = []
	# classes: 0=light, 1=streaming, 2=sensitive
	for i,app in enumerate(original_workload):
		# Be careful and DO NOT use bench_id, it is changed inside get_optimal_clustering
		app.original_app.ix = i
		app_class[app.name] = classification[i]
		if app_class[app.name] == 0: # light-sharing benchmarks
			light_sharing_apps.append(app)
		else:
			workload.append(app)

	##Add a ls app to get the min...
	if len(light_sharing_apps)==1:
		app = light_sharing_apps.pop(0)
		workload.append(app)

	# Update number of cores per LLC to reserve space for light sharing (balanced approach)
	cores_per_llc=user_options["cores_per_llc"]
	nr_light_per_set=len(light_sharing_apps)//cores_per_llc

	uopt_patched=dict(user_options)
	uopt_patched["opt_mapping"]=True
	max_core_groups=(len(original_workload)+cores_per_llc-1)//cores_per_llc
	uopt_patched["max_core_groups"]=max_core_groups
	## Dejando máximo....
	#if nr_light_per_set!=0:
	#	uopt_patched["cores_per_llc"]=cores_per_llc-nr_light_per_set

	if algorithm=="popt-map-unf":	
		(patched_workload,per_app_ways,per_app_masks,cluster_id,total_branches,metadata)=get_optimal_clustering_bf(workload,unfairness_max_throughput,False,nr_ways,max_bandwidth,multiprocessing=multiprocessing,user_options=uopt_patched)
	elif algorithm=="popt-map-stp":
		(patched_workload,per_app_ways,per_app_masks,cluster_id,total_branches,metadata)=get_optimal_clustering_bf(workload,throughput,True,nr_ways,max_bandwidth,multiprocessing=multiprocessing,user_options=uopt_patched)
	elif algorithm=="popt-map-cov":
		(patched_workload,per_app_ways,per_app_masks,cluster_id,total_branches,metadata)=get_optimal_clustering_bf(workload,cov_unfairness_metric,False,nr_ways,max_bandwidth,multiprocessing=multiprocessing,user_options=uopt_patched)
	else:
		(patched_workload,per_app_ways,per_app_masks,cluster_id,total_branches,metadata)=get_optimal_clustering_bf(workload,opt_spec[0],opt_spec[1],nr_ways,max_bandwidth,multiprocessing=multiprocessing,user_options=uopt_patched)	

	nr_clusters = max(cluster_id)+1
	clusters = [[] for _ in range(nr_clusters)]
	nr_light_sharing = len(light_sharing_apps)
	cluster_way = {}
	cluster_mask = {}


	for idxc in range(nr_clusters):
		cluster_i = cluster_id.index(idxc)
		cluster_way[idxc] = per_app_ways[cluster_i]
		cluster_mask[idxc] = per_app_masks[cluster_i]

	for i,app in enumerate(patched_workload):
		app.llc_id=cluster_id[i]
		clusters[cluster_id[i]].append(app.original_app)

	## Assign lights  (load_balance())
	## Traverse streaming first ... [Require sorted stuff]

	## Go ahead and assign in RR the LIGHT SHARING apps among rest
	while nr_light_sharing > 0:
		app = light_sharing_apps.pop(0)
		## Find free cluster
		clx_id=0
		while (clx_id < nr_clusters) and (len(clusters[clx_id])>=cores_per_llc):
			clx_id+=1
		assert clx_id<nr_clusters, "All clusters are full and there is no gap for a light sharing program"
		clusters[clx_id].append(app)
		app.llc_id=clx_id
		nr_light_sharing = nr_light_sharing - 1

	# Re-build solution with additional light apps
	per_app_masks = [None] * len(original_workload)
	per_app_ways = [None] * len(original_workload)
	cluster_id = [None] * len(original_workload)
	patched_workload = [None] * len(original_workload)

	for idxc,cluster in enumerate(clusters):
		patched_cluster = get_scaled_properties_cluster(cluster,nr_ways)
		for app in patched_cluster:
			app.llc_id=idxc
			orig_i = app.original_app.ix
			per_app_ways[orig_i] = cluster_way[idxc]
			per_app_masks[orig_i] = cluster_mask[idxc]
			patched_workload[orig_i] = app
			cluster_id[orig_i] = idxc

	return (patched_workload, per_app_ways, per_app_masks, cluster_id, total_branches, metadata)


def cpa(workload,nr_ways,uoptions={}):
	uoptions.setdefault("benchmark_categories","./data/classification_cpa_voltav2.csv")
	class_dict=get_benchmark_categories_from_file(uoptions["benchmark_categories"])
	classification=[class_dict[bench.name] for bench in workload]
	ad_hoc_cpa_clustering={(0,0):[0x7ff], (1,0): [0x3f,0x7f0], (0,1): [0x0FF,0x700], 
				(2,0): [0x01F,0x7f0,0x7f0] ,(1,1): [0x01F,0x7f0,0x700], (0,2):[0x01F,0x7f0,0x700]  , 
				(3,0): [0x00F,0x7F8,0x7F8,0x7F8], (2,1):[0x00F,0x7F8,0x7F8,0x600], (1,2):[0x00F,0x7F8,0x7F8,0x600],
				(0,3):[0x00F,0x7F8,0x600,0x600]  }
	SENSITIVE_CLASS=0
	MEDIUM_CLASS=1
	BULLY_CLASS=2
	SQUANDERER_CLASS=3
	NON_CRITICAL_CLASS=4
	STREAMING_CLASS=5
	critical_apps=[]
	sensitive_apps=[]
	medium_apps=[]
	other_apps=[]
	squanderer_apps=[]
	streaming_apps=[]
	cluster_id_app={}

	for i,app in enumerate(workload):
		app.bench_id=i ## for sorting
		app_class=classification[i]
		if app_class==SENSITIVE_CLASS:
			critical_apps.append(app)
			sensitive_apps.append(app)
		elif app_class==MEDIUM_CLASS:
			critical_apps.append(app)
			medium_apps.append(app)	
		elif app_class==BULLY_CLASS or app_class==NON_CRITICAL_CLASS:
			other_apps.append(app)
		elif app_class==SQUANDERER_CLASS:
			squanderer_apps.append(app)
		elif app_class==STREAMING_CLASS:
			streaming_apps.append(app)

	nr_critical=len(critical_apps)
	nr_sensitive=len(sensitive_apps)
	nr_medium=len(medium_apps)	
	nr_squanderers=len(squanderer_apps)
	nr_streaming=len(streaming_apps)

	clusters=[]
	nr_clusters=0

	if nr_critical==0 or nr_critical>3:
		index=(0,0)
		cluster_masks=list(ad_hoc_cpa_clustering[index])
		all_other_apps=other_apps+critical_apps
		clusters.append(all_other_apps)
		non_critical_cluster=all_other_apps
		## Update cluster id app
		for app in all_other_apps:
			cluster_id_app[app.bench_id]=nr_clusters
		non_critical_cluster_id=nr_clusters
		nr_clusters+=1

	else:
		index=(nr_sensitive,nr_medium)
		cluster_masks=list(ad_hoc_cpa_clustering[index])

		#First other apps.
		non_critical_cluster=other_apps
		clusters.append(other_apps)
		for app in other_apps:
			cluster_id_app[app.bench_id]=nr_clusters
		non_critical_cluster_id=nr_clusters
		nr_clusters+=1


		## Remaining clusters...	
		for app in sensitive_apps:
			clusters.append([app])
			cluster_id_app[app.bench_id]=nr_clusters
			nr_clusters+=1

		for app in medium_apps:
			clusters.append([app])
			cluster_id_app[app.bench_id]=nr_clusters
			nr_clusters+=1

	if nr_squanderers>0:
		cluster_masks.append(0x1)
		clusters.append(squanderer_apps)
		for app in squanderer_apps:
			cluster_id_app[app.bench_id]=nr_clusters
		nr_clusters+=1

	if nr_streaming>0:
		streaming_assigned=0
		possible_parts=[0x1,0x2]

		for i in range(min(nr_streaming,2-nr_squanderers)):
			app=streaming_apps[i]
			cluster_masks.append(possible_parts[i+nr_squanderers])
			clusters.append([app])
			streaming_assigned+=1
			cluster_id_app[app.bench_id]=nr_clusters
			nr_clusters+=1

		##The remaining ones go to non-critical partition
		for j in range(len(streaming_apps)-streaming_assigned):
			app=streaming_apps[j+streaming_assigned]
			non_critical_cluster.append(app)
			cluster_id_app[app.bench_id]=non_critical_cluster_id
		

	## Calculate ways in each cluster
	cluster_ways=[]
	for cluster_mask in cluster_masks:
		nr_ways_clust=0
		for i in range(nr_ways):
			if cluster_mask & (0x1<<i):
				nr_ways_clust+=1
		cluster_ways.append(nr_ways_clust)


	## Ahora hay que ordenar 
	per_app_masks = []
	per_app_ways = []
	apps_ix = []
	cluster_id = []
	for i_app,app in enumerate(workload):
		per_app_masks.append(hex(cluster_masks[cluster_id_app[i_app]]))
		per_app_ways.append(cluster_ways[cluster_id_app[i_app]])
		cluster_id.append(cluster_id_app[i_app])		

	return (per_app_ways,per_app_masks,cluster_id)


def apply_ucp_slowdown_cluster(clusters,nr_ways,max_bandwidth=float('Inf')):
	apps=[app for cluster in clusters for app in cluster.apps]
	norm_cluster=[cluster.apps for cluster in clusters]

	if len(clusters)==1:
		partitioning=[nr_ways]
	else:
		partitioning=ucp_clustering_solution(norm_cluster,nr_ways)	


	## Determine Per-application way-assignment
	per_app_partitioning=[]
	per_app_ways=[]
	slowdown_vector=[]
	slowdown_list=[] ## tuple list (slowdown, idx_cluster) 

	for idxc,nr_ways_cluster in enumerate(partitioning):
		## Traverse applications in a cluster to determine 
		cl=clusters[idxc]
		for idxa,app in enumerate(cl.apps):
			# Now stores the cluster ways so it can later retrieve the scaled value from properties
			per_app_partitioning.append((app,nr_ways_cluster,idxc))
			per_app_ways.append(nr_ways_cluster)
			slowdown=app.get_metric_table("slowdown")[nr_ways_cluster]

			slowdown_vector.append(slowdown)
			slowdown_list.append((slowdown,app.bench_id,idxc,idxa))


	#slowdown_vector=get_slowdown_vector(	apps,per_app_ways,max_bandwidth)

	if len(clusters) >  nr_ways:
		unfairness_value = (1000,-1000) ## (Unf,-STP)
	else:
		## Calculate throughput
		throughput_value=0
		for slowdown in slowdown_vector:
			throughput_value += 1/slowdown	
		
		unfairness = max(slowdown_vector) / min(slowdown_vector)

		unfairness_value=(unfairness,-throughput_value)


	## The third parameter is the only change (remains unused no BW)
	return (partitioning,per_app_partitioning,slowdown_list, unfairness_value)


def evaluate_clustering(clusters,nr_ways,partitioning, max_bandwidth=float('Inf')):
	apps=[app for cluster in clusters for app in cluster.apps]
	norm_cluster=[cluster.apps for cluster in clusters]

	## Determine Per-application way-assignment
	per_app_partitioning=[]
	per_app_ways=[]
	slowdown_vector=[]
	slowdown_list=[] ## tuple list (slowdown, idx_cluster) 

	for idxc,nr_ways_cluster in enumerate(partitioning):
		## Traverse applications in a cluster to determine 
		cl=clusters[idxc]
		for idxa,app in enumerate(cl.apps):
			# Now stores the cluster ways so it can later retrieve the scaled value from properties
			per_app_partitioning.append((app,nr_ways_cluster,idxc))
			per_app_ways.append(nr_ways_cluster)
			slowdown=app.get_metric_table("slowdown")[nr_ways_cluster]

			slowdown_vector.append(slowdown)
			slowdown_list.append((slowdown,app.bench_id,idxc,idxa))


	#slowdown_vector=get_slowdown_vector(	apps,per_app_ways,max_bandwidth)

	if len(clusters) >  nr_ways:
		unfairness_value = (1000,-1000) ## (Unf,-STP)
	else:
		## Calculate throughput
		throughput_value=0
		for slowdown in slowdown_vector:
			throughput_value += 1/slowdown	
		
		unfairness = max(slowdown_vector) / min(slowdown_vector)

		unfairness_value=(unfairness,-throughput_value)


	## The third parameter is the only change (remains unused no BW)
	return (partitioning,per_app_partitioning,slowdown_list, unfairness_value)



def pair_clustering_core(workload,nr_ways,max_bandwidth=float('Inf'),debugging=False,uoptions={}):
	curClusters = []
	total_apps = len(workload)	

	uoptions.setdefault("check_max_slowdown",False)
	uoptions.setdefault("use_ucp_slowdown",True)
	uoptions.setdefault("use_original_way_count",True)
	uoptions.setdefault("verbose",False)


#	check_max_slowdown=True
#	use_ucp_slowdown=True
#	use_original_way_count=True
	
	##Create distance matrix
	distance_matrix=np.full((total_apps,total_apps),float('Inf'))

	## Create lists for curves - triplets (combined_curve,buckets,partitioned_curve)
	curve_set=[]
	for i in range(total_apps):
		row=[]
		for j in range(total_apps):
			row.append(None)
		curve_set.append(row)


	## Initial clustering solution
	for i,app in enumerate(workload):
		buckets = [[nw] for nw in range(1,nr_ways + 1)]
		curve = app.get_metric_table("slowdown").values
		app.bench_id = i
		cl=Cluster([app],curve,buckets)
		curClusters.append(cl)

	## Calculate initial distances and curves
	## In this case distance matrix is complete for efficiency reasons
	## Don't forget to update it right away
	for i,clusterRef in enumerate(curClusters):
		idx_i=clusterRef.idx_cluster
		for j in range(i+1,len(curClusters)):
			cluster=curClusters[j]
			idx_j=cluster.idx_cluster
			(distance,combined_curve,buckets,partitioned_curve,distance_curve)=clusterRef.slowdown_distance(cluster)
		
			## Store in matrices (Symmetric stuff in matrix)
			distance_matrix[idx_i,idx_j]=distance
			distance_matrix[idx_j,idx_i]=distance

			## Just stored for i < j
			curve_set[idx_i][idx_j]=(combined_curve,buckets,partitioned_curve,distance_curve,clusterRef,cluster)


	(partitioning,per_app_partitioning,slowdown_list,unfairness_value)=apply_ucp_slowdown_cluster(curClusters,nr_ways,max_bandwidth)
	solutions = [(list(curClusters),(list(partitioning),per_app_partitioning,slowdown_list,unfairness_value))]

	if uoptions["verbose"]:

		print("***Initial assignment ***")
		print(curClusters, partitioning)
		print(slowdown_list)
		print(unfairness_value)
		print(distance_matrix)
		print("***********************************")

	it=0

	while len(curClusters) > 1:

		## Boolean vector that marks which clusters are left to explore
		nr_clusters=len(curClusters)
		merged=[False for i in range(nr_clusters)]
		sorted_slowdowns=sorted(slowdown_list,reverse=True)


		for (slowdown,idx_app,idxc,idxa) in sorted_slowdowns:
			it=it+1
			#norm_index=idxc*nr_clusters+idxa
			if merged[idxc]:
				continue

			# Step 1: find closest cluster in terms of distance and do not merge if distance >=0

			idx_min=distance_matrix[idxc].argmin()
			distance=distance_matrix[idxc,idx_min]


			if distance>=0:
				continue

			# Step 2: See if merging the clusters contributes to reducing the maximum slowdown 	
			# Normalize indexes
			(min_i,min_j)= (idxc,idx_min) if idx_min > idxc else (idx_min,idxc)

			(combined_curve,buckets,partitioned_curve,distance_curve,clusterA,clusterB)=curve_set[min_i][min_j]

			# How many ways have both clusters so far
			nr_combined_ways=partitioning[min_i]+partitioning[min_j]
 	
			if uoptions["check_max_slowdown"]:

				max_slowdown_combined=combined_curve[nr_combined_ways]

				max_slowdown_partitioned=partitioned_curve[nr_combined_ways]

				## It is not worth merging clusters
				if max_slowdown_combined > max_slowdown_partitioned:
					continue
					
			## Merge 2 closest clusters
			new_cluster=merge_clusters(clusterA,clusterB,combined_curve,buckets)

			if uoptions["verbose"]:
				print("***Iteration %d***" % it)
				print(distance_matrix)
				print(distance)
				print("Buckets -> ", [("%.3f" % v1,"%.3f" % v2)  for v1,v2 in  buckets]) 
				print("Combined -> ",[ ("%.3f" % v) for v in  combined_curve])
				print("Partitioned -> ",[ ("%.3f" % v) for v in  partitioned_curve])
				print("App0 : -> ",[ ("%.3f" % v) for v in  new_cluster.apps[0].original_app.get_metric_table("slowdown")])
				print("App0': -> ",[ ("%.3f" % v) for v in  new_cluster.apps[0].get_metric_table("slowdown")])
				print("App1 : -> ",[ ("%.3f" % v) for v in  new_cluster.apps[0].original_app.get_metric_table("slowdown")])
				print("App1': -> ",[ ("%.3f" % v) for v in  new_cluster.apps[0].get_metric_table("slowdown")])
				print("Selection:",clusterA,clusterB)


			## Clear items from distance matrix
			merged[min_i]=merged[min_j]=True

			# This way they will not be reconsidered as candidates
			distance_matrix[min_i,:]=float('Inf')
			distance_matrix[:,min_i]=float('Inf')
			distance_matrix[min_j,:]=float('Inf')
			distance_matrix[:,min_j]=float('Inf')

			## Create a copy of cur clusters to save solution (to see how unfairness goes)

			# Remove prev clusters and overwrite location
			curClusters[min_i]=None
			curClusters[min_j]=None
			curClusters[new_cluster.idx_cluster]=new_cluster

			## Filter
			potClustering=list(filter(lambda x: x is not None, curClusters))

			## Evaluate clustering (clusterA,cluster)
			if uoptions["use_ucp_slowdown"]:
				(ppartitioning,pper_app_partitioning,pslowdown_list,punfairness_value)=apply_ucp_slowdown_cluster(potClustering,nr_ways,max_bandwidth)
				## Update partitioning for assesment
				if not uoptions["use_original_way_count"]:
					for i in range(len(potClustering)):
						clust_idx=potClustering[i].idx_cluster
						partitioning[clust_idx]=ppartitioning[i]

			else:
				# Update partitioning vector (Overwrite only what matters)
				partitioning[new_cluster.idx_cluster]=nr_combined_ways

				merged_partitioning=[partitioning[cluster.idx_cluster] for cluster in potClustering]


				(ppartitioning,pper_app_partitioning,pslowdown_list,punfairness_value)=evaluate_clustering(potClustering,nr_ways,merged_partitioning,max_bandwidth)		

			if uoptions["verbose"]:
				print(pper_app_partitioning)
				print(potClustering, ppartitioning)
				print(pslowdown_list)
				print(punfairness_value)
				print("***********************************")

			solutions.append((potClustering,(ppartitioning,pper_app_partitioning,pslowdown_list,punfairness_value)))

		break	

	best_score=(float('Inf'),-float('Inf'))
	a=3
	best_k=0

	for (idx,sol) in enumerate(solutions):
		(clusters,data)=sol
		(partitioning,per_app_partitioning,slowdown_list,score)=data

		if score < best_score:
			best_score=score
			best_k=idx


	assert sum(solutions[best_k][1][0])==nr_ways, "The algorithm is assigning all available ways to clusters"

	return (solutions,best_k)



## Precondition !merged[idx_app]
def determine_slowdown_reductions(workload,idx_app,merged,curve_set,ucp_partitioning):
	slowdown_red=np.full(len(workload),float('Inf'))

	app_assigned_ways=ucp_partitioning[idx_app]
	for i,app in enumerate(workload):
		if merged[i]:
			slowdown_red[i]=float('Inf')
		elif i==idx_app:
			## What if it was possible to transfer one way
			slowdown_curve=app.properties.slowdown
			if app_assigned_ways==len(slowdown_curve): ## No more ways
				slowdown_red[i]=float('Inf')
			else:
				slowdown_red[i]=slowdown_curve[app_assigned_ways+1]-slowdown_curve[app_assigned_ways]
		else: ## Calculate local distance
			(min_i,min_j)= (i,idx_app) if idx_app > i else (idx_app,i)
			(combined_curve,buckets,partitioned_curve,distance_curve,clusterA,clusterB,mergedCluster)=curve_set[min_i][min_j]
			total_ways=ucp_partitioning[min_i]+ucp_partitioning[min_j]
			patched_apps=mergedCluster.apps
			slowdown_combined_single_part=patched_apps[0].properties["slowdown"][total_ways]+patched_apps[1].properties["slowdown"][total_ways]
			slowdown_partitioned=workload[min_i].properties["slowdown"][ucp_partitioning[min_i]]+workload[min_j].properties["slowdown"][ucp_partitioning[min_j]] 
			#slowdown_red[i]=slowdown_combined_single_part-partitioned_curve[total_ways-1]
			slowdown_red[i]=slowdown_combined_single_part-slowdown_partitioned
			#slowdown_red[i]=combined_curve[total_ways-1]-partitioned_curve[total_ways-1]

	return slowdown_red


def pair_clustering_core2(workload,nr_ways,max_bandwidth=float('Inf'),debugging=False,uoptions={}):
	curClusters = []
	total_apps = len(workload)	

	uoptions.setdefault("check_max_slowdown",False)
	uoptions.setdefault("use_ucp_slowdown",True)
	uoptions.setdefault("use_original_way_count",True)
	uoptions.setdefault("verbose",False)

#	check_max_slowdown=True
#	use_ucp_slowdown=True
#	use_original_way_count=True
	
	##Create distance matrix
	distance_matrix1w_nomask=np.full((total_apps,total_apps),float('Inf'))
	## Turn matrix into a masked array
	distance_matrix1w=np.ma.array(distance_matrix1w_nomask,mask=False)

	## Create lists for curves - triplets (combined_curve,buckets,partitioned_curve)
	curve_set=[]
	for i in range(total_apps):
		row=[]
		for j in range(total_apps):
			row.append(None)
		curve_set.append(row)


	## Initial clustering solution
	for i,app in enumerate(workload):
		buckets = [[nw] for nw in range(1,nr_ways + 1)]
		curve = app.get_metric_table("slowdown").values
		app.bench_id = i
		cl=Cluster([app],curve,buckets)
		curClusters.append(cl)

	## Calculate initial distances and curves
	## In this case distance matrix is complete for efficiency reasons
	## Don't forget to update it right away
	for i,clusterRef in enumerate(curClusters):
		idx_i=clusterRef.idx_cluster
		for j in range(i+1,len(curClusters)):
			cluster=curClusters[j]
			idx_j=cluster.idx_cluster
			(distance,combined_curve,buckets,partitioned_curve,distance_curve)=clusterRef.slowdown_distance(cluster)
			
			## Obtain merged cluster beforehand as in the C implementation
			new_cluster=merge_clusters(clusterRef,cluster,combined_curve,buckets)

			## What if we mix these applications in a 1-way partition?
			## It is just not possible to determine Slowdown 1-way ... we use slowdown 2-way
			## Unfortunately in C we do not have patched apps....
			slowdown_partitioned=workload[i].properties["slowdown"][1]+workload[j].properties["slowdown"][1] 
			patched_apps=new_cluster.apps
			slowdown_combined_single_part=patched_apps[0].properties["slowdown"][1]+patched_apps[1].properties["slowdown"][1]


			if uoptions["verbose"]:
				print("DISTANCE 1W is %f" % (slowdown_combined_single_part-slowdown_partitioned))
			
			distance_matrix1w[idx_i,idx_j]=distance_matrix1w[idx_j,idx_i]=slowdown_combined_single_part-slowdown_partitioned

			## Just stored for i < j
			curve_set[idx_i][idx_j]=(combined_curve,buckets,partitioned_curve,distance_curve,clusterRef,cluster,new_cluster)


	#
	ucp_unf_part=ucp_unfairness(workload,nr_ways)
	(partitioning,per_app_partitioning,slowdown_list,unfairness_value)=evaluate_clustering(curClusters,nr_ways,ucp_unf_part,max_bandwidth)
	solutions = [(list(curClusters),(list(partitioning),per_app_partitioning,slowdown_list,unfairness_value))]

	#(partitioning,per_app_partitioning,slowdown_list,unfairness_value)=apply_ucp_slowdown_cluster(curClusters,nr_ways,max_bandwidth)
	#solutions = [(list(curClusters),(list(partitioning),per_app_partitioning,slowdown_list,unfairness_value))]


	## KEEP TRACK OF INITIAL UCP PARTITIONING
	initial_partitioning=list(partitioning)


	## Clear 1w matrix with info from the freshly obtained partition
	for i,ways_assigned in enumerate(partitioning):
		if ways_assigned!=1:
			# This way they will not be reconsidered as candidates
			distance_matrix1w[i,:]=float('Inf')
			distance_matrix1w[:,i]=float('Inf')
	

	if uoptions["verbose"]:

		print("***Initial assignment PAIR CLUSTERING 2***")
		print(curClusters, partitioning)
		print(slowdown_list)
		print(unfairness_value)
		print(distance_matrix1w)
		print("***********************************")

	it=0


	## Boolean vector that marks which clusters are left to explore
	nr_clusters=len(curClusters)
	merged=[False for i in range(nr_clusters)]
	sorted_slowdowns=sorted(slowdown_list,reverse=True)


	for (slowdown,idx_app,idxc,idxa) in sorted_slowdowns:
		it=it+1
		clusters_created=0
		#norm_index=idxc*nr_clusters+idxa
		if merged[idxc]:
			continue

		# Step 1: Determine slowdown merge benefits for this application
		slowdown_reductions=determine_slowdown_reductions(workload,idx_app,merged,curve_set,partitioning)

		if uoptions["verbose"]:
			print("** Slowdown reduction for app %s with index %i***" % (workload[idxc].name,idxc))
			print(slowdown_reductions)
			print("***********************************")


		## Determinar mínimos 1w pero excluyendo idxc
		distance_matrix1w.mask[idxc,:]=True
		distance_matrix1w.mask[:,idxc]=True
		minval=np.amin(distance_matrix1w)

		if minval!=float('Inf'):
			result = np.where(distance_matrix1w == np.amin(distance_matrix1w)) 
			idx1w_a=result[0][0]
			idx1w_b=result[0][1]
			## Calculate cost
			slowdown_reductions[idx_app]+=minval
		else:
			slowdown_reductions[idx_app]=float('Inf') ## It is not possible to steal a way

		## What's best. Combine or add a way
		idx_min=slowdown_reductions.argmin()
		benefit=slowdown_reductions[idx_min]

		## Undo mask
		distance_matrix1w.mask[:,idxc]=False
		distance_matrix1w.mask[idxc,:]=False

		## Its  best to steal a way ... 
		if idx_min==idx_app and benefit < 0:
			## Let's steal a way
			## Merge 2 1-way clusters
			min_i=idx1w_a
			min_j=idx1w_b

			distance_matrix1w[min_i,:]=float('Inf')
			distance_matrix1w[:,min_i]=float('Inf')
			distance_matrix1w[min_j,:]=float('Inf')
			distance_matrix1w[:,min_j]=float('Inf')

			(combined_curve,buckets,partitioned_curve,distance_curve,clusterA,clusterB,new_cluster)=curve_set[min_i][min_j]

			merged[min_i]=merged[min_j]=True
			#merged[idx_app]=True

			# Remove prev clusters and overwrite location
			curClusters[min_i]=None
			curClusters[min_j]=None
			curClusters[new_cluster.idx_cluster]=new_cluster

			## Filter
			potClustering=list(filter(lambda x: x is not None, curClusters))

			# Update partitioning vector (Overwrite only what matters)
			partitioning[new_cluster.idx_cluster]=1
			partitioning[idx_app]=partitioning[idx_app]+1 ## New single-app cluster steals a new way
			clusters_created+=1
			slowdown_reductions[min_i]=float('Inf')
			slowdown_reductions[min_j]=float('Inf')

		## Can we still combine?
		slowdown_reductions[idx_app]=float('Inf')
		idx_min=slowdown_reductions.argmin()
		benefit=slowdown_reductions[idx_min]

		if benefit<0 and clusters_created==0:
			# Normalize indexes
			(min_i,min_j)= (idxc,idx_min) if idx_min > idxc else (idx_min,idxc)


			# How many ways have both clusters so far
			nr_combined_ways=partitioning[min_i]+partitioning[min_j]

			(combined_curve,buckets,partitioned_curve,distance_curve,clusterA,clusterB,new_cluster)=curve_set[min_i][min_j]

			## Merge 2 closest clusters

			## Clear items from distance matrix
			merged[min_i]=merged[min_j]=True

			# This way they will not be reconsidered as candidates
			distance_matrix1w[min_i,:]=float('Inf')
			distance_matrix1w[:,min_i]=float('Inf')
			distance_matrix1w[min_j,:]=float('Inf')
			distance_matrix1w[:,min_j]=float('Inf')			
			## Create a copy of cur clusters to save solution (to see how unfairness goes)

			# Remove prev clusters and overwrite location
			curClusters[min_i]=None
			curClusters[min_j]=None
			curClusters[new_cluster.idx_cluster]=new_cluster

			## Filter
			potClustering=list(filter(lambda x: x is not None, curClusters))

			# Update partitioning vector (Overwrite only what matters)
			partitioning[new_cluster.idx_cluster]=nr_combined_ways
			clusters_created+=1

		
		if clusters_created:
			if 	clusters_created==2:
				print("BIEN JUNTA2!")
			merged_partitioning=[partitioning[cluster.idx_cluster] for cluster in potClustering]

			(ppartitioning,pper_app_partitioning,pslowdown_list,punfairness_value)=evaluate_clustering(potClustering,nr_ways,merged_partitioning,max_bandwidth)		

			if uoptions["verbose"]:
				print(pper_app_partitioning)
				print(potClustering, merged_partitioning)
				print(slowdown_list)
				print(pslowdown_list)
				print(ppartitioning)
				print(punfairness_value)
				print("***********************************")

			solutions.append((potClustering,(ppartitioning,pper_app_partitioning,pslowdown_list,punfairness_value)))


	best_score=(float('Inf'),-float('Inf'))
	best_k=0

	for (idx,sol) in enumerate(solutions):
		(clusters,data)=sol
		(partitioning,per_app_partitioning,slowdown_list,score)=data

		if score < best_score:
			best_score=score
			best_k=idx


	assert sum(solutions[best_k][1][0])==nr_ways, "The algorithm is not assigning all available ways to clusters %d " % sum(solutions[best_k][1][0])

	return (solutions,best_k)

def pair_clustering(workload,nr_ways,max_bandwidth=float('Inf'),debugging=False,user_options={}):

	if "pc2" in user_options:
		(solutions,k)=pair_clustering_core2(workload,nr_ways,max_bandwidth,debugging,uoptions=user_options)
	else:
		(solutions,k)=pair_clustering_core(workload,nr_ways,max_bandwidth,debugging,uoptions=user_options)

	if debugging:
		for sol in solutions:
			clusters,(partitioning,per_app_partitioning,bw_app,score) = sol
			print(clusters, score)		
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


def slowdown_clustering_core(workload,nr_ways,max_bandwidth=float('Inf'),debugging=False,uoptions={}):
	curClusters = []
	total_apps = len(workload)	
	
	uoptions.setdefault("verbose",False)
	uoptions.setdefault("absolute_distance",False)
	
	## Create distance matrix
	distance_matrix=np.full((total_apps,total_apps),float('Inf'))

	## Create lists for curves - triplets (combined_curve,buckets,partitioned_curve)
	curve_set=[]
	for i in range(total_apps):
		row=[]
		for j in range(total_apps):
			row.append(None)
		curve_set.append(row)


	for i,app in enumerate(workload):
		buckets = [[nw] for nw in range(1,nr_ways + 1)]
		curve = app.get_metric_table("slowdown").values
		## Supercritical
		app.bench_id = i
		cl=Cluster([app],curve,buckets)
		curClusters.append(cl)

	## Calculate initial distances and curves
	for i,clusterRef in enumerate(curClusters):
		idx_i=clusterRef.idx_cluster
		for j in range(i+1,len(curClusters)):
			cluster=curClusters[j]
			idx_j=cluster.idx_cluster
			(distance,combined_curve,buckets,partitioned_curve, distance_curve)=clusterRef.slowdown_distance(cluster)

			##Store in matrices (Symmetric stuff in matrix)
			distance_matrix[idx_i,idx_j]=distance
			distance_matrix[idx_j,idx_i]=distance

			curve_set[idx_i][idx_j]=(combined_curve,buckets,partitioned_curve,clusterRef,cluster)

	(partitioning,per_app_partitioning,slowdown_list,unfairness_value)=apply_ucp_slowdown_cluster(curClusters,nr_ways,max_bandwidth)
	solutions = [(list(curClusters),(partitioning,per_app_partitioning,slowdown_list,unfairness_value))]

	if uoptions["verbose"]:
		print("***Initial assignment ***")
		print(curClusters, partitioning)
		print(slowdown_list)
		print(unfairness_value)
		print("***********************************")

	# (partitioning,per_app_partitioning,bw_app,sp_aggregate) <- determine_best_part()

	it=1
	max_iterations=len(curClusters)-1

	for i in range(max_iterations):
	
		## Find min distance
		## Absolute min distance
		if uoptions["absolute_distance"]:		
			min_i,min_j=np.unravel_index(distance_matrix.argmin(), distance_matrix.shape)
			distance=distance_matrix[min_i,min_j]
			(combined_curve,buckets,partitioned_curve,clusterA,clusterB)=curve_set[min_i][min_j]
		else:
			sorted_slowdowns=sorted(slowdown_list,reverse=True)

			## Search one that can be merged with others
			for (slowdown,idx_app,idxc,idxa) in sorted_slowdowns:

				# Step 1: find closest cluster in terms of distance and do not merge if distance >=0

				idx_min=distance_matrix[idxc].argmin()

				distance=distance_matrix[idxc,idx_min]

				if distance<0:
					break


			if distance >= 0:
				#Nothing can be merged (exit outer loop)
				break 

			# Step 2: Retrieve curves of clusters to be merged	
			# Normalize indexes
			(min_i,min_j)= (idxc,idx_min) if idx_min > idxc else (idx_min,idxc)

			(combined_curve,buckets,partitioned_curve,clusterA,clusterB)=curve_set[min_i][min_j]


		## Determine location of prev_clusters (traverse all)
		min_idx=[0,0]
		for idx,cluster in enumerate(curClusters):
			if cluster.idx_cluster == clusterA.idx_cluster:
				min_idx[0]=idx
			elif cluster.idx_cluster == clusterB.idx_cluster:
				min_idx[1]=idx

		if uoptions["verbose"]:		
			print("***Iteration %d***" % it)
			print(idx_min)
			print(min_idx)
			print(distance_matrix)
			print(distance)
			print("Selection:",clusterA,"(%d)" % min_idx[0],clusterB,"(%d)" % min_idx[1])
			print("***********************************")
		it=it+1

		## Merge 2 closest clusters
		new_cluster=merge_clusters(clusterA,clusterB,combined_curve,buckets)

		## Sorted is critical!
		min_idx=sorted(min_idx)

		## Remove prev clusters
		del curClusters[min_idx[0]]
		del curClusters[min_idx[1]-1]

		## One of the indexes is going away (the minimum one). Update the distance matrix to discard values in the future
		removed_index=clusterA.idx_cluster if clusterA.idx_cluster>clusterB.idx_cluster else clusterB.idx_cluster
		distance_matrix[removed_index,:]=float('Inf')
		distance_matrix[:,removed_index]=float('Inf')

		# Add merged cluster ...
		curClusters.insert(0,new_cluster)

		## Very important to update slowdown_list 
		(ppartitioning,pper_app_partitioning,slowdown_list,punfairness_value)=apply_ucp_slowdown_cluster(curClusters,nr_ways,max_bandwidth)

		if uoptions["verbose"]:
			print(pper_app_partitioning)
			print(curClusters, ppartitioning)
			print(slowdown_list)
			print(punfairness_value)
			print("***********************************")

		solutions.append((list(curClusters),(ppartitioning,pper_app_partitioning,slowdown_list,punfairness_value)))	

		## Calculate new distance if this is not the last iteration for sure
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
		
					(distance,combined_curve,buckets,partitioned_curve)=cluster_i.slowdown_distance(cluster_j)
		
					##Store in matrices
					distance_matrix[idx_i,idx_j]=distance
					distance_matrix[idx_j,idx_i]=distance
					curve_set[idx_i][idx_j]=(combined_curve,buckets,partitioned_curve,cluster_i,cluster_j)


	## Determine BEST K
	best_score=(float('Inf'),float('Inf'))
	
	best_k=0

	for (idx,sol) in enumerate(solutions):
		(clusters,data)=sol
		(partitioning,per_app_partitioning,bw_app,score)=data

		if score < best_score:
			best_score=score
			best_k=idx

	#assert sum(solutions[best_k][1][0])==nr_ways, "The algorithm is assigning all available ways to clusters"

	return (solutions,best_k)



def slowdown_clustering(workload,nr_ways,max_bandwidth=float('Inf'),debugging=False,user_options={}):

	(solutions,k)=slowdown_clustering_core(workload,nr_ways,max_bandwidth,debugging,uoptions=user_options)

	if debugging:
		for sol in solutions:
			clusters,(partitioning,per_app_partitioning,bw_app,score) = sol
			print(clusters, score)		
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


# Function that trivially maps applications to 'nr_clusters' clusters.
# For example, supplied with: workload = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']; nr_clusters = 3
# returns [['A', 'D', 'G'], ['B', 'E', 'H'], ['C', 'F']]
def initialize_clusters_dio(workload, nr_clusters):
	return  [workload[i::nr_clusters] for i in range(nr_clusters)]

# Returns a list of the same length as the number of clusters and containing the cumulative 
# sum of the llc miss rates of the apps in each cluster
def compute_cluster_llc_rate(cluster):
	nr_ways=11 # for now
	llc_miss_rate = 0
	for i,app in enumerate(cluster):
		llc_miss_rate += app.properties['llcmpkc'][nr_ways]

	return llc_miss_rate


def dio(workload, nr_ways, max_bandwidth, uoptions={}): 
	uoptions.setdefault("nr_core_groups", 4)
	uoptions.setdefault("cores_per_llc", 4)

	nr_core_groups = uoptions["nr_core_groups"]
	cores_per_llc = uoptions["cores_per_llc"]

	# Create an initial trivial clustering to get the patched workload (whose properties are scaled in order to consider contention)
	initial_clustering = initialize_clusters_dio(workload, nr_core_groups)
	patched_clustering, patched_workload = get_patched_clustering_and_workload(initial_clustering, nr_ways)
	'''
	print("Initial workload: ", workload, "\n")
	print("Initial clustering: ", initial_clustering, "\n")
	print("Patched Workload: ", patched_workload, "\n")
	'''
	clusters = [[] for i in range(nr_core_groups)]
	llc_miss_rates = [0.0 for cluster in clusters]

	# Sort apps in descending llcmpki
	patched_workload.sort(key=lambda x: x.properties['llcmpki'][nr_ways], reverse=True)
	'''
	print("Patched workload (sorted by llcmpki) ", patched_workload, "\n")
	for k, app in enumerate(patched_workload):
		print(k, ") ", app, " - llcmpki = ", app.properties['llcmpki'][nr_ways])
	'''

	for i, app in enumerate(patched_workload):
		#print(i, ") ", app, "\ninitial: clusters ", clusters, " llc_miss_rates ", llc_miss_rates)

		dst_id = llc_miss_rates.index(min(llc_miss_rates))		

		if len(clusters[dst_id]) < cores_per_llc:
			pass
		else:
			llc_miss_rates[dst_id] = float('inf')
			if min(llc_miss_rates) == float('inf'):
				break
			dst_id = llc_miss_rates.index(min(llc_miss_rates))

		clusters[dst_id].append(app.original_app)
		app.original_app.llc_id=dst_id

		llc_miss_rates[dst_id]= llc_miss_rates[dst_id] + app.properties['llcmpki'][nr_ways]
		if len(clusters[dst_id]) == cores_per_llc:
			llc_miss_rates[dst_id] = float('inf')

		#print("final: clusters ", clusters, " llc_miss_rates ", llc_miss_rates, "\n")

	#print("\nfinal clusters: ", clusters, "\n")

	clustering_sol = [nr_ways for cluster in clusters]

	return normalize_output_for_clustering_solution(workload,clusters,clustering_sol,nr_ways)

# Function which computes the number of ways of an app that results in a slowdown closest to a given one
def get_nr_ways_matching_slowdown(app, slowdown_value):
	# 'app.properties' is a pandas df where the index are the number of ways
	slowdown_array=app.properties["slowdown"]
	abs_diff = np.abs(slowdown_array-slowdown_value)
	nr_ways = abs_diff.loc[abs_diff == abs_diff.min()].index.values[0]

	return nr_ways


# Prototype of load-balancer which takes into account both memory BW and LLC usage to perform the mapping, and then applies LFOC+ to each CCX 
def llc_bw_balancer(workload, nr_ways, max_bandwidth, uoptions={}):
	uoptions.setdefault("nr_core_groups", 2)
	uoptions.setdefault("cores_per_llc", 4)
	uoptions.setdefault("cache_part", None)
	# Invoke test/sim.py with -O bechmark_categories=<path2categoryfile> to change the default
	uoptions.setdefault("benchmark_categories", "data/rome-classification.csv") 
	uoptions.setdefault("target_slowdown", 1.05) # target slowdown value for cache-sensitives apps to balance LLC contention 
	uoptions.setdefault("balancing_mode",2) #0-> sensitive, 1 -> streaming,  2-> auto  

	nr_core_groups = uoptions["nr_core_groups"]
	cores_per_llc = uoptions["cores_per_llc"]
	classification_file = uoptions["benchmark_categories"]
	cache_part=uoptions["cache_part"]
	target_slowdown = uoptions["target_slowdown"]
	balancing_mode = uoptions["balancing_mode"]

	# Create an initial trivial clustering to get the patched workload (whose properties are scaled in order to consider contention)
	initial_clustering = initialize_clusters_dio(workload, nr_core_groups)
	patched_clustering, patched_workload = get_patched_clustering_and_workload(initial_clustering, nr_ways)

	clusters = [[] for i in range(nr_core_groups)]
	min_nr_ways = [0 for cluster in clusters]
	bandwidths = [0.0 for cluster in clusters]

	# Extract benchmark categories
	apps_categories = get_benchmark_categories_from_file(classification_file, label_format=True)

	light_sharings = []
	streamings = []
	cache_sensitives = []

	for app in patched_workload:
		if apps_categories[app.name] == "light-sharing":
			light_sharings.append(app)
		elif apps_categories[app.name] == "streaming":
			streamings.append(app)
		else:
			cache_sensitives.append(app)

	# Sort streaming apps based on descending memory BW
	streamings.sort(key=(lambda x: x.properties["bandwidth_mbps"][nr_ways]), reverse=True)

	# Sort cache-sensitives apps based on descending nr_ways necessary to obtain a slowdown closest to the target
	cache_sensitives.sort(key=(lambda x: get_nr_ways_matching_slowdown(x.original_app, target_slowdown)), reverse=True)
	'''
	print("Cache-sensitives apps (number of ways):")
	for k, app in enumerate(cache_sensitives):
		print(k, ") ", app, " - nr_ways = ", get_nr_ways_matching_slowdown(app.original_app, target_slowdown), " - bw = ", app.properties["bandwidth_mbps"][nr_ways])
	
	print("Streaming apps (sorted by bw):")
	for k, app in enumerate(streamings):
		print(k, ") ", app, " - nr_ways = ", get_nr_ways_matching_slowdown(app.original_app, target_slowdown), " - bw = ", app.properties["bandwidth_mbps"][nr_ways])

	print("Light-sharing apps (not sorted):")
	for k, app in enumerate(light_sharings):
		print(k, ") ", app, " - nr_ways = ", get_nr_ways_matching_slowdown(app.original_app, target_slowdown), " - bw = ", app.properties["bandwidth_mbps"][nr_ways])
	'''


	# Apply a dio-like balancer and place apps by distributing BWs/LLCmiss rates for streaming/cache-sensitives apps among clusters 

	if balancing_mode==0:
		apps_to_balance = [cache_sensitives, streamings]
	elif balancing_mode==1:
		apps_to_balance = [streamings, cache_sensitives]
	else:
		if len(streamings) >= len(cache_sensitives):
			apps_to_balance = [streamings, cache_sensitives]
		else:
			apps_to_balance = [cache_sensitives, streamings]

	# Apply a dio-like balancer and place apps by distributing BWs/LLCmiss rates for streaming/cache-sensitives apps among clusters 
	for apps in apps_to_balance:
		if apps == streamings:
			aggregate_metric = bandwidths
		elif apps == cache_sensitives:
			aggregate_metric = min_nr_ways

		for app in apps:
			dst_id = aggregate_metric.index(min(aggregate_metric))      

			if len(clusters[dst_id]) < cores_per_llc:
				pass
			else:
				aggregate_metric[dst_id] = float('inf')
				if min(aggregate_metric) == float('inf'):
					break
				dst_id = aggregate_metric.index(min(aggregate_metric))

			clusters[dst_id].append(app.original_app)
			app.original_app.llc_id=dst_id

			min_nr_ways[dst_id] = min_nr_ways[dst_id] + get_nr_ways_matching_slowdown(app.original_app, target_slowdown)
			bandwidths[dst_id] = bandwidths[dst_id] + app.properties["bandwidth_mbps"][nr_ways]
			if len(clusters[dst_id]) == cores_per_llc:
				aggregate_metric[dst_id] = float('inf')

	'''
	print("clusters", clusters)
	print("Bandwidths per-cluster", bandwidths)
	print("min nr ways per-cluster", min_nr_ways)
	print("================================================================================")
	'''

	# Finally place light-sharing programs in the remaining available spots
	for app in light_sharings:
		for cluster_id,cluster in enumerate(clusters):
			if len(cluster) < cores_per_llc:
				cluster.append(app.original_app)
				app.original_app.llc_id = cluster_id
				break

	clustering_sol = [nr_ways for cluster in clusters]

	'''
	print("clusters", clusters)
	print("Bandwidths per-cluster", bandwidths)
	print("min nr ways per-cluster", min_nr_ways)
	print("================================================================================")
	'''

	if cache_part=="lfoc+":
		lfoc_uoptions=uoptions.copy()
		lfoc_uoptions["use_pair_clustering"]=True
		lfoc_uoptions["simple_output"]=True

		clusters_part=[]
		way_distribution=[]
		## Get the original apps and apply LFOC for each LLC_ID
		for dio_cluster in clusters:
			subworkload=[app.original_app for app in dio_cluster]
			## Invoke LFOC+ for that
			(subclusters,llc_way_distr)=lfoc(subworkload,nr_ways,max_bandwidth,lfoc_uoptions)	
			clusters_part.extend(subclusters)
			way_distribution.extend(llc_way_distr)

		clusters=clusters_part
		clustering_sol=way_distribution

	return normalize_output_for_clustering_solution(workload,clusters,clustering_sol,nr_ways)

# Prototype of load-balancer which takes into account both memory BW and LLC usage to perform the mapping, and then applies LFOC+ to each CCX 
def llc_bw_balancer_compositions(workload, nr_ways, max_bandwidth, uoptions={}):
	uoptions.setdefault("nr_core_groups", 4)
	uoptions.setdefault("cores_per_llc", 4)
	uoptions.setdefault("cache_part", None)
	# Invoke test/sim.py with -O bechmark_categories=<path2categoryfile> to change the default
	uoptions.setdefault("benchmark_categories", "data/rome-classification.csv") 
	uoptions.setdefault("target_slowdown", 1.05) # target slowdown value for cache-sensitives apps to balance LLC contention 

	nr_core_groups = uoptions["nr_core_groups"]
	cores_per_llc = uoptions["cores_per_llc"]
	classification_file = uoptions["benchmark_categories"]
	cache_part=uoptions["cache_part"]
	target_slowdown = uoptions["target_slowdown"]

	# Create an initial trivial clustering to get the patched workload (whose properties are scaled in order to consider contention)
	initial_clustering = initialize_clusters_dio(workload, nr_core_groups)
	patched_clustering, patched_workload = get_patched_clustering_and_workload(initial_clustering, nr_ways)

	clusters = [[] for i in range(nr_core_groups)]
	min_nr_ways = [0 for cluster in clusters]
	bandwidths = [0.0 for cluster in clusters]

	# Extract benchmark categories
	apps_categories = get_benchmark_categories_from_file(classification_file, label_format=True)

	light_sharings = []
	streamings = []
	cache_sensitives = []

	for app in patched_workload:
		if apps_categories[app.name] == "light-sharing":
			light_sharings.append(app)
		elif apps_categories[app.name] == "streaming":
			streamings.append(app)
		else:
			cache_sensitives.append(app)

	# Sort streaming apps based on descending memory BW
	streamings.sort(key=(lambda x: x.properties["bandwidth_mbps"][nr_ways]), reverse=True)

	# Sort cache-sensitives apps based on descending nr_ways necessary to obtain a slowdown closest to the target
	cache_sensitives.sort(key=(lambda x: get_nr_ways_matching_slowdown(x.original_app, target_slowdown)), reverse=True)

	'''
	print("Cache-sensitives apps (number of ways):")
	for k, app in enumerate(cache_sensitives):
		print(k, ") ", app, " - nr_ways = ", get_nr_ways_matching_slowdown(app.original_app, target_slowdown), " - bw = ", app.properties["bandwidth_mbps"][nr_ways])
	
	print("Streaming apps (sorted by bw):")
	for k, app in enumerate(streamings):
		print(k, ") ", app, " - nr_ways = ", get_nr_ways_matching_slowdown(app.original_app, target_slowdown), " - bw = ", app.properties["bandwidth_mbps"][nr_ways])

	print("Light-sharing apps (not sorted):")
	for k, app in enumerate(light_sharings):
		print(k, ") ", app, " - nr_ways = ", get_nr_ways_matching_slowdown(app.original_app, target_slowdown), " - bw = ", app.properties["bandwidth_mbps"][nr_ways])
	'''

	# Apply a dio-like balancer and place apps by distributing BWs/LLCmiss rates for streaming/cache-sensitives apps among clusters 
	if len(streamings) >= len(cache_sensitives):
		apps_to_balance = [streamings, cache_sensitives]
	else:
		apps_to_balance = [cache_sensitives, streamings]

	for apps in apps_to_balance:
		if apps == streamings:
			aggregate_metric = bandwidths
		elif apps == cache_sensitives:
			aggregate_metric = min_nr_ways

		for app in apps:
			dst_id = aggregate_metric.index(min(aggregate_metric))      

			if len(clusters[dst_id]) < cores_per_llc:
				pass
			else:
				aggregate_metric[dst_id] = float('inf')
				if min(aggregate_metric) == float('inf'):
					break
				dst_id = aggregate_metric.index(min(aggregate_metric))

			clusters[dst_id].append(app.original_app)
			app.original_app.llc_id=dst_id

			min_nr_ways[dst_id] = min_nr_ways[dst_id] + get_nr_ways_matching_slowdown(app.original_app, target_slowdown)
			bandwidths[dst_id] = bandwidths[dst_id] + app.properties["bandwidth_mbps"][nr_ways]
			if len(clusters[dst_id]) == cores_per_llc:
				aggregate_metric[dst_id] = float('inf')

	'''
	print("clusters", clusters)
	print("Bandwidths per-cluster", bandwidths)
	print("min nr ways per-cluster", min_nr_ways)
	print("================================================================================")
	'''

	# Finally place light-sharing programs in the remaining available spots
	for app in light_sharings:
		for cluster_id,cluster in enumerate(clusters):
			if len(cluster) < cores_per_llc:
				cluster.append(app.original_app)
				app.original_app.llc_id = cluster_id
				break

	clustering_sol = [nr_ways for cluster in clusters]

	'''
	print("clusters", clusters)
	print("Bandwidths per-cluster", bandwidths)
	print("min nr ways per-cluster", min_nr_ways)
	print("================================================================================")
	'''

	if cache_part=="lfoc+":
		lfoc_uoptions=uoptions.copy()
		lfoc_uoptions["use_pair_clustering"]=True
		lfoc_uoptions["simple_output"]=True

		clusters_part=[]
		way_distribution=[]
		## Get the original apps and apply LFOC for each LLC_ID
		for dio_cluster in clusters:
			subworkload=[app.original_app for app in dio_cluster]
			## Invoke LFOC+ for that
			(subclusters,llc_way_distr)=lfoc(subworkload,nr_ways,max_bandwidth,lfoc_uoptions)	
			clusters_part.extend(subclusters)
			way_distribution.extend(llc_way_distr)

		clusters=clusters_part
		clustering_sol=way_distribution

	return normalize_output_for_clustering_solution(workload,clusters,clustering_sol,nr_ways)



def dino(workload, nr_ways, max_bandwidth, uoptions={}): 
	uoptions.setdefault("nr_core_groups", 4)
	uoptions.setdefault("cores_per_llc", 4)
	uoptions.setdefault("high_llct",10)
	uoptions.setdefault("low_llct",5)
	uoptions.setdefault("cache_part",None)

	nr_core_groups = uoptions["nr_core_groups"]
	cores_per_llc = uoptions["cores_per_llc"]
	high_llct= uoptions["high_llct"]
	low_llct= uoptions["low_llct"]
	cache_part=uoptions["cache_part"]

	# Create an initial trivial clustering to get the patched workload (whose properties are scaled in order to consider contention)
	initial_clustering = initialize_clusters_dio(workload, nr_core_groups)
	patched_clustering, patched_workload = get_patched_clustering_and_workload(initial_clustering, nr_ways)

	clusters = [[] for i in range(nr_core_groups)]
	llc_miss_rates = [0.0 for cluster in clusters]

	# Sort apps in descending llcmpki
	patched_workload.sort(key=(lambda x: x.properties['llcmpki'][nr_ways]), reverse=True)
	'''
	print("Patched workload (sorted by llcmpki) ", patched_workload, "\n")
	for k, app in enumerate(patched_workload):
		print(k, ") ", app, " - llcmpki = ", app.properties['llcmpki'][nr_ways])
	'''

	turtles=[]
	devils=[]
	superdevils=[]

	for app in patched_workload:
		llcmpki=app.properties['llcmpki'][nr_ways]
		if llcmpki >= high_llct:
			superdevils.append(app)
		elif llcmpki >= low_llct:
			devils.append(app)
		else:
			turtles.append(app)

	left=True
	for apps in [superdevils,turtles,devils]:
		for i in range(len(apps)):
			if left:
				app=apps[i].original_app
			else:
				app=patched_workload[len(patched_workload)-i-1]	

			original_app=app.original_app
			dst_id = llc_miss_rates.index(min(llc_miss_rates))		

			if len(clusters[dst_id]) == cores_per_llc:
				llc_miss_rates[dst_id] = float('inf')
				if min(llc_miss_rates) == float('inf'):
					break
				dst_id = llc_miss_rates.index(min(llc_miss_rates))

			## It is important to add the original app to the cluster and not the patched one
			## To ensure that patching stuff is done correctly
			clusters[dst_id].append(app.original_app)
			original_app.llc_id=dst_id
			llc_miss_rates[dst_id]= llc_miss_rates[dst_id] + app.properties['llcmpki'][nr_ways]
			if len(clusters[dst_id]) == cores_per_llc:
				llc_miss_rates[dst_id] = float('inf')	

		left=not left

	clustering_sol = [nr_ways for cluster in clusters]


#	if cache_part=="lfoc+":
#		lfoc_uoptions=uoptions.copy()
#		lfoc_uoptions["use_pair_clustering"]=True
#		lfoc_uoptions["simple_output"]=True

#		clusters_part=[]
#		way_distribution=[]
		## Get the original apps and apply LFOC for each LLC_ID
#		for dio_cluster in clusters:
#			subworkload=[app.original_app for app in dio_cluster]
			## Invoke LFOC+ for that
#			(subclusters,llc_way_distr)=lfoc(subworkload,nr_ways,max_bandwidth,lfoc_uoptions)	
#			clusters_part.extend(subclusters)
#			way_distribution.extend(llc_way_distr)

#		clusters=clusters_part
#		clustering_sol=way_distribution


	return normalize_output_for_clustering_solution(workload,clusters,clustering_sol,nr_ways)

def dio_opt(workload, nr_ways, max_bandwidth, uoptions={}):
	uoptions.setdefault("nr_core_groups", 3)
	uoptions.setdefault("cores_per_llc", 4)

	nr_core_groups = uoptions["nr_core_groups"]
	cores_per_llc = uoptions["cores_per_llc"]

	patched_workload = get_scaled_properties_cluster(workload, nr_ways, metric='llcmpkc')

	nr_clusters = nr_core_groups
	apps_per_cluster = cores_per_llc

	clusters = [[] for n in range(nr_clusters)]

	clusters = initialize_clusters_dio(workload, apps_per_cluster)
	clusters_iterator = generate_possible_mappings(list(range(len(patched_workload))), 1, apps_per_cluster, nr_core_groups)

	distances = []
	distance = 0

	for k,clustering in enumerate(clusters_iterator):
		cluster_size = len(clustering)
		clustering = get_application_clustering_from_app_ids(clustering, workload)
		#print(k, ") ", clustering)
		for i in range(cluster_size-1):
			cluster = list(clustering[i])
			llc_miss_rate1 = compute_cluster_llc_rate(cluster)

			for j in range(i+1, cluster_size):
				cluster = list(clustering[j])
				cluster = get_application_clustering_from_app_ids(clustering[j], workload)
				llc_miss_rate2 = compute_cluster_llc_rate(cluster)
				distance += pow(llc_miss_rate2-llc_miss_rate1, 2)
				distances.append(distance)
		#print(k, ") ", distance)

	#print("Minimum distance ", min(distances))

	return (patched_workload,per_app_ways,per_app_masks,cluster_id)


def trivial_mapping(workload, nr_ways, max_bandwidth, uoptions={}):
	#uoptions are passed directly to ./test/sim.py, e.g. -O nr_core_groups 3 -O cores_per_llc 3
	uoptions.setdefault("nr_core_groups", 2) 
	uoptions.setdefault("cores_per_llc", 4)

	nr_core_groups = uoptions["nr_core_groups"]
	cores_per_llc = uoptions["cores_per_llc"]

	trivial_mapping = [workload[cores_per_llc*i:cores_per_llc*(i+1)] for i in range(nr_core_groups)]
	patched_clustering, patched_workload = get_patched_clustering_and_workload(trivial_mapping, nr_ways)

	clusters = patched_clustering

	clustering_sol = [nr_ways for cluster in trivial_mapping]

	for ccx,cluster in enumerate(clusters):
		for app in cluster:
			app.llc_id = ccx

	return normalize_output_for_clustering_solution(workload,clusters,clustering_sol,nr_ways)


def add_experimental_algorithms():
	add_extra_algorithms("dunn","bw-on-demand","class","class-stream","class-sens","lfoc","poptc","poptc-unf","poptc-stp","poptc-slow","pair-clustering",
		"scluster","lfoc+","lfoc-opt","ucp-unfairness","cpa","popt-map","popt-map-unf","popt-map-stp","popt-map-cov", "dio", "dino", "llc_bw_balancer", "llc_bw_balancer_compositions", "trivial_mapping")


def invoke_extra_algorithm(algorithm,workload,nr_ways,max_bandwidth,parallel=False,debugging=False,uoptions={},opt_spec=(unfairness,"False")):
	uoptions.setdefault("req_threshold",21)
	uoptions.setdefault("harmless_ratio",1)
	uoptions.setdefault("misses_gap_ratio_threshold",2)

	if algorithm=="dunn":
		patched_workload=workload
		(per_app_ways,per_app_masks,cluster_id)=get_schedule_dunn(workload, nr_ways)
	elif algorithm=="bw-on-demand":
		solution=on_demand_partitioning(workload,nr_ways,max_bandwidth)
		(patched_workload,per_app_ways,per_app_masks,cluster_id)=(workload,solution,get_partition_masks(solution),[i for i in range(1,len(workload)+1)])
	elif algorithm=="ucp-unfairness":
		solution=ucp_unfairness(workload,nr_ways,uoptions=uoptions)
		(patched_workload,per_app_ways,per_app_masks,cluster_id)=(workload,solution,get_partition_masks(solution),[i for i in range(1,len(workload)+1)])				
	elif algorithm=="class":
		(patched_workload,per_app_ways,per_app_masks,cluster_id)=classification_clustering(workload,nr_ways,max_bandwidth,0)
	elif algorithm=="class-stream":
		(patched_workload,per_app_ways,per_app_masks,cluster_id)=classification_clustering(workload,nr_ways,max_bandwidth,1)
	elif algorithm=="class-sens":
		(patched_workload,per_app_ways,per_app_masks,cluster_id)=classification_clustering(workload,nr_ways,max_bandwidth,2)		
	elif algorithm=="lfoc":
		(patched_workload,per_app_ways,per_app_masks,cluster_id)=lfoc(workload,nr_ways,max_bandwidth,uoptions)
	elif algorithm=="lfoc+":
		uoptions["use_pair_clustering"]=True
		(patched_workload,per_app_ways,per_app_masks,cluster_id)=lfoc(workload,nr_ways,max_bandwidth,uoptions)		
		uoptions["use_pair_clustering"]=False
	elif algorithm=="lfoc-opt":
		uoptions["opt"]=True
		(patched_workload,per_app_ways,per_app_masks,cluster_id)=lfoc(workload,nr_ways,max_bandwidth,uoptions)			
	elif (algorithm=="poptc-unf") or (algorithm=="poptc-stp") or (algorithm=="poptc-slow") or (algorithm=="poptc"):
		(patched_workload, per_app_ways, per_app_masks, cluster_id, total_branches, metadata)=clustering_pseudo_optimal(algorithm,workload,nr_ways,max_bandwidth,multiprocessing=parallel,user_options=uoptions,opt_spec=opt_spec)
	elif (algorithm=="popt-map-unf") or (algorithm=="popt-map-stp") or (algorithm=="popt-map-cov") or (algorithm=="popt-map"):
		(patched_workload, per_app_ways, per_app_masks, cluster_id, total_branches, metadata)=mapping_pseudo_optimal(algorithm,workload,nr_ways,max_bandwidth,multiprocessing=parallel,user_options=uoptions,opt_spec=opt_spec)		
	elif algorithm=="pair-clustering":
		(patched_workload,per_app_ways,per_app_masks,cluster_id)=pair_clustering(workload,nr_ways,max_bandwidth,debugging,user_options=uoptions) 
	elif algorithm=="scluster":
		(patched_workload,per_app_ways,per_app_masks,cluster_id)=slowdown_clustering(workload,nr_ways,max_bandwidth,debugging,user_options=uoptions) 
	elif algorithm=="cpa":
		patched_workload=workload
		(per_app_ways,per_app_masks,cluster_id)=cpa(workload, nr_ways,uoptions)
	elif algorithm=="dio":
		(patched_workload,per_app_ways,per_app_masks,cluster_id) = dio(workload, nr_ways, max_bandwidth,uoptions)
	elif algorithm=="dino":
		(patched_workload,per_app_ways,per_app_masks,cluster_id) = dino(workload, nr_ways, max_bandwidth,uoptions)
	elif algorithm=="llc_bw_balancer":
		(patched_workload,per_app_ways,per_app_masks,cluster_id) = llc_bw_balancer(workload, nr_ways, max_bandwidth,uoptions) 
	elif algorithm=="llc_bw_balancer_compositions":
		(patched_workload,per_app_ways,per_app_masks,cluster_id) = llc_bw_balancer_compositions(workload, nr_ways, max_bandwidth,uoptions)
	elif algorithm=="trivial_mapping":
		(patched_workload,per_app_ways,per_app_masks,cluster_id) = trivial_mapping(workload, nr_ways, max_bandwidth, uoptions)		
	else:
		print("algorithm not valid: %s" % algorithm)
		exit(1)		
	return (patched_workload,per_app_ways,per_app_masks,cluster_id)
