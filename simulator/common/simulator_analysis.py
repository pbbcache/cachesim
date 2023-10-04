# -*- coding: utf-8 -*-

#
# simulator_analysis.py
#
##############################################################################
#
# Copyright (c) 2020 Juan Carlos Saez<jcsaezal@ucm.es>
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
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import string
import os
from simulator_core import *
from simulator_results import *
import matplotlib.colors as colpackage
import matplotlib.cm as cmx


###############################################################
### AUXILIARY FUNCTIONS (NOT FOR PUBLIC AND SYSTEMATIC USE) ###
###############################################################

def apply_offset(row,offset,prefix_len=1):
	workload_name=row["W#"]
	number=int(workload_name[prefix_len:])+offset
	return workload_name[0:prefix_len]+str(number)

## Usage: df["W#"]=df.apply(prefix_workload,1,args=(new_prefix,old_prefix_len))
def prefix_workload(row,new_prefix,old_prefix_len=1):
	workload_name=row["W#"]
	return new_prefix+workload_name[old_prefix_len:] 

## Usage: df["Algorithm"]=df.apply(algorithm_suffix,1,args=(suffix,))
def algorithm_suffix(row,suffix):
	return row["Algorithm"]+suffix
	

#################################
### MODULE's PUBLIC FUNCTIONS ###
#################################

def simlog_to_df(filename,format="legacy"):
	"""Builds/Returns a Dataframe with the contents of the simulator logfile (filename) 

	Note that the logfile must comply with PBBCache's "table" output format (-f table)
	"""
	if format=="df":
		df=pd.read_csv(filename, delim_whitespace=True)
		df.rename(columns={"workload":"W#","scheme":"Algorithm","id":"BenchID","app":"Name","property":"Property","value":"Value"},inplace=True)
		return df
	else:
		return pd.read_csv(filename,sep='\\s+',engine='python')

def add_algorithm_suffix(df,suffix):
	"""Add a suffix to every partitioning algorithm listed in the Dataframe.

	df must have been read previously with simlog_to_df()
	"""
	df["Algorithm"]=df.apply(algorithm_suffix,1,args=(suffix,))
 

def add_prefix_workload(df,new_prefix,old_prefix_len=1):
	"""Add a prefix to every workload listed in the Dataframe.

	df must have been read previously with simlog_to_df()
	"""	
	df["W#"]=df.apply(prefix_workload,1,args=(new_prefix,old_prefix_len))
	
def apply_workload_offset(df,offset,prefix_len=1):
	"""Add a prefix to every workload listed in the Dataframe.

	df must have been read previously with simlog_to_df()
	"""	
	df["W#"]=df.apply(apply_offset,1,args=(offset,prefix_len))
	
	
def concat_simlog_dfs(dfs):
	"""Given a list of Simlog Dataframe returns a new Simlog Dataframe
	 with the resulting concatenation.

	The dataframes in the list must have been read previously with simlog_to_df()
	"""		
	## Concat
	new_df=pd.concat(dfs, axis = 0,ignore_index=True)
	## Sort (workload,algorithm,benchmark_id)
	new_df.sort_values(by=["W#","Algorithm","BenchID"],inplace=True)
	return new_df


def df_to_simlog(df,filename):
	"""Given a Simlog Dataframe, generate a textual/table representation 
	of it and store it in filename.
	"""		
	## Determine workload prefix automatically:
	workload_name=df["W#"].values[0]
	len_pref=len(workload_name)-1
	
	formatted_string="%-Xs %-19s %-9s %-17s %-13s %-15s %-17s %-17s %-17s\n"
	formatted_string=formatted_string.replace("X",str(len_pref+3))
	
	with open(filename, 'w+') as file:
		line=formatted_string % ("W#","Algorithm","BenchID","Name","Cluster","Mask/NR_WAYS","SlowdownNB/ANTT","STP","Slowdown")
		file.write(line)   
		for index, row in df.iterrows():
			line=formatted_string % (row["W#"],row["Algorithm"],row["BenchID"],row["Name"],row["Cluster"],row["Mask/NR_WAYS"],row["SlowdownNB/ANTT"],row["STP"],row["Slowdown"])
			file.write(line)


def build_chart_data(df,norm=None): 
	"""From a SimLog Dataframe, return 2 summary dataframes with STP and
	 Unfairness for the various algorithms. 

	The norm argument is used to specify the algorithm used to normalize the data
	 (e.g., whirlpool-c) if necessary. The returned dataframes --(stp,unf) -- are
	 already suitable for the creation of a bar chart (e.g. stp.plot.bar(rot=0,figsize=size,)) 
	"""	  
	overall=df[df.Name=="OVERALL"]
	## Get STP and Slowdown alone
	selection=overall[["W#","Algorithm","STP","Slowdown"]].set_index("W#")
	stp=selection.pivot_table(values='STP',columns='Algorithm',index=selection.index)
	unf=selection.pivot_table(values='Slowdown',columns='Algorithm',index=selection.index) #index=["W#"])

	if norm:
		for column in stp.columns:
		## Skip
			if column!=norm:
				stp[column]=stp[column]/stp[norm]
	
		stp[norm]=1.0
		for column in unf.columns:
			## Skip
			if column!=norm:
				unf[column]=unf[column]/unf[norm]

		unf[norm]=1.0        

	return (stp,unf)


def build_chart_data_metrics(df,metrics,norm=None): 
	"""From a SimLog Dataframe with the new (dF format), return a dict of dataframes 
	for the different metrics passed as a parameter

	The norm argument is used to specify the algorithm used to normalize the data
	 (e.g., whirlpool-c) if necessary. The dataframes in the returned dict -- are
	 already suitable for the creation of a bar chart
	  (e.g. dictionary["stp"].plot.bar(rot=0,figsize=size,)) 
	"""	  
	overall=df[df.Name=="OVERALL"]
	table_metrics=overall.Property.unique()
	
	cmetrics={}
 
	for metric in metrics:
		if not metric in table_metrics:
			print("No such metric: ",metric)
			return None
		
		dfm=overall[overall.Property==metric]
 
		## Get STP and Slowdown alone
		selection=dfm[["W#","Algorithm","Value"]].set_index("W#")
		## Critical as we used masks as property
		selection=selection.astype({'Value':'float'})
		metricdf=selection.pivot_table(values='Value',columns='Algorithm',index=selection.index)
		cmetrics[metric]=metricdf

		if norm:
			for column in metricdf.columns:
			## Skip
				if column!=norm:
					metricdf[column]=metricdf[column]/metricdf[norm]
			## Normalize the baseline too
			metricdf[norm]=1.0

	return cmetrics

tikz_template=r"""
\documentclass[tikz]{standalone}
\usepackage[utf8x]{inputenc}
\usepackage{lmodern,textcomp}
\usepackage{tikz}
\usepackage{spot}
\usetikzlibrary{calc,arrows}

\pgfdeclarelayer{background}
\pgfdeclarelayer{foreground}
\pgfsetlayers{background,main,foreground}

\begin{document}
\begin{tikzpicture}[scale=0.6,transform shape]
\def\boxH{0.8cm}
\def\boxW{1.3cm}
\def\wayH{0.6cm}
\def\wayW{1.01cm}
\def\coresep{0.5cm}
%\def\ccolor1{green!50!white, blue!20, red!20, yellow, gray!50!green, black, orange!20
$COLOR_DEFS$


\begin{scope}[draw=black,font=\small,
core/.style={draw=black,rectangle,rounded corners,minimum height=\boxH,minimum width=\boxW},
way/.style={draw=black,fill=gray!20,rectangle,minimum height=\wayH,minimum width=\wayW},
flecha/.style={>=stealth',black,semithick},
etiq/.style={text=black,align=center,font=\footnotesize}]


\coordinate (pointer) at (0,0);
\foreach \way in {$WAY_RANGE$} 
{
	\node (way\way) [way,above right] at (pointer) {};
	\node [etiq,below] at (way\way.south) {Way\way};
	\coordinate (pointer) at (way\way.south east);
}

\node [below,yshift=-0.5cm,text=black] at ($(way0.south)!0.5!(way$LAST_WAY$.south)$)  {Unfairness=$UNF$ -- STP=$STP$}; 

\coordinate (pointer) at ($(way0.north west)+(0,0.5cm)$);

%1/Name/\color
\foreach \core/\app/\acolor/\slowdown in {$BENCH_SPEC$} 
{
	\node (core\core) [core,above right,fill=\acolor] at (pointer) {};
	%\node [etiq,below] at (core\core.south) {Core \core};
	\node [etiq,above] at (core\core.north) {\app \\ (\slowdown)};
	\coordinate (pointer) at ($(core\core.south east)+(\coresep,0)$);
}

%0/\ccolor1%
\foreach \way/\wcolor in {$WAY_SPEC$} 
{
	\node (wayc\way) [way,fill=\wcolor] at (way\way) {};
}


\end{scope}


\end{tikzpicture}
\end{document}
"""

## Auxiliary function to turn DF clustering solution into 
## The necessary data structures to make the diagram
## generation easier 
def build_clustering_solution(df,max_ways=11):
	nr_apps=df
	diff_clusters=[int(data) for data in df.Cluster.unique()]
	nr_clusters=len(diff_clusters)
	cluster_zero=(0 in diff_clusters)
	#print nr_clusters
	clusters=[(i,0,0) for i in range(nr_clusters)] ##(id,way_mask,nr_ways)
	way_cluster_mapping=[-1 for i in range(max_ways)] ## ncluster each way i belongs to
	app_cluster_mapping=[] ## (appname,id_cluster)
	
	## Traverse dataframe 
	for index, row in df.iterrows():
		id_cluster=int(row.Cluster)
		## normalize
		if not cluster_zero:
			id_cluster-=1
		app_cluster_mapping.append((row.Name,id_cluster))
		## It is a new cluster?

		#print id_cluster
		(id,way_mask,nr_ways)=clusters[id_cluster]
		if nr_ways==0:
			cluster_info=row["Mask/NR_WAYS"]
			cluster_info=cluster_info[:-1] ## Remove last)
			tokens=cluster_info.split("(")
			assert(len(tokens)==2)
			mask=int(tokens[0],16)
			nr_ways=int(tokens[1])
			clusters[id_cluster]=(id_cluster,mask,nr_ways)
			## Traverse mask to populate way_cluster_mapping
			for i in range(max_ways):
				if (mask & (1<<i)):
					way_cluster_mapping[i]=id_cluster
		
	return (clusters,way_cluster_mapping,app_cluster_mapping)


def generate_colors_from_colormap(nr_colors,colormap):
	nr_colors=10
	values=[i for i in range(nr_colors)]
	jet = cm = plt.get_cmap(colormap) 
	cNorm  = colpackage.Normalize(vmin=0, vmax=values[-1])
	scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)
	print(scalarMap.get_clim())

	colors = []
	for idx in range(10):
		colorVal = scalarMap.to_rgba(values[idx])
		colorRGB=[int(255*v) for v in colorVal]
		del colorRGB[3]
		colors.append(colorRGB)

	return colors
			
def generate_clustering_chart(df,workload,algorithm,**kwargs):
	"""Given a Simlog Dataframe, generate a TIKZ diagram for the selected
	workload and algorithm.
	
	The diagram is generated in the "./charts" default directory (this setting 
	may be changed with the directory parameter). The number of ways in the LLC
	(max_ways)  must be also specified if !=1.

	The function returns a tuple (tex_file,pdf_file,png_file) with the path of
	the three files that have just been generated (3 versions of the diagram). In
	a Jupyter notebook png_file can be used to embed the image in the notebook as 
	follows: 
	
	from IPython.display import Image
	
	(tex,pdf,png)=generate_clustering_chart(all_sensitive,"S10-1","ucp-slowdown")
	Image(filename=png,width=500)
	"""		
	## Determine workload prefix automatically:	
	kwargs.setdefault("max_ways",11)
	kwargs.setdefault("directory","./charts")
	kwargs.setdefault("density",300)
	kwargs.setdefault("resize",100)
	kwargs.setdefault("colormap",'Paired')

	## Retrieve major params
	directory=kwargs["directory"] 
	max_ways=kwargs["max_ways"]

	assert(os.path.exists(directory))

	## Build filenames
	rawfilename="%s_%s.tex" % (workload,algorithm)
	rootfile="%s/%s_%s" % (directory,workload,algorithm)
	filename="%s.tex" % (rootfile)
	pdffile="%s.pdf"  % (rootfile)
	pngfile="%s.png"  % (rootfile)

	clustering_sol=df[(df["W#"]==workload) & (df["Algorithm"]==algorithm)]
	filtered_sol=clustering_sol[clustering_sol["Name"]!="OVERALL"]
	metrics=clustering_sol[clustering_sol["Name"]=="OVERALL"]
	unfairness=metrics["Slowdown"].values[0]
	slowdowns=filtered_sol["Slowdown"].values
	stp=metrics["STP"].values[0]

	#print clustering_sol
	## Get key data structures
	(clusters,way_cluster_mapping,app_cluster_mapping)=build_clustering_solution(filtered_sol,max_ways)

	## Define color_defs
	color_defs=""
	colors=generate_colors_from_colormap(len(clusters),kwargs["colormap"])
	#colors=[r"green!50!white",r"blue!20",r"red!20",r"yellow",r"gray!50!green",r"black",r"orange!20",r"blue!40",r"gray!25",r"white"]

	for i in range(len(clusters)):
		color_defs+=r'\definecolor{ccolor'
		color_str=','.join([str(c) for c in colors[i]])
		color_defs+=string.ascii_uppercase[i]+r"}{RGB}{" + color_str +"}\n"
	#print color_defs   
	
	## Define way ranges
	way_ranges=','.join([str(i) for i in range(len(way_cluster_mapping))])
	#print way_ranges
	
	## Last way
	last_way=len(way_cluster_mapping)-1
	
	## Bench spec
	bench_spec_list=[ r"%i/%s/ccolor%s/%.2f" % (i,app,string.ascii_uppercase[id_cluster],slowdowns[i]) for i,(app,id_cluster) in enumerate(app_cluster_mapping)]
	bench_spec=',\n'.join(bench_spec_list)
	#print bench_spec    
	
	## WAY SPECS
	way_spec_list=[ r"%i/ccolor%s" % (i,string.ascii_uppercase[id_cluster]) for i,id_cluster in enumerate(way_cluster_mapping)]
	way_spec=',\n'.join(way_spec_list)
	#print way_spec
	
	contents=tikz_template.replace(r'$COLOR_DEFS$',color_defs).replace(r'$WAY_RANGE$',way_ranges).replace(r'$LAST_WAY$',str(last_way)).replace(r'$BENCH_SPEC$',bench_spec).replace(r'$WAY_SPEC$',way_spec).replace(r'$UNF$',str(round(unfairness,3))).replace(r'$STP$',str(round(stp,3)))
	
	with open(filename, "w") as text_file:
		text_file.write(contents)

	## Compile figure 
	print("Building figure %s.{tex,pdf,png} ..." % rootfile,
	os.system("cd %s; pdflatex %s; cd - ; convert -density %i -resize %i%% %s %s" % (directory,rawfilename,kwargs["density"],kwargs["resize"],pdffile,pngfile)))
	print('Done!')
	return (filename,pdffile,pngfile)


## Sample function to generate STP/UNF Dataframes
def build_charts_norm(stp,unf,file_prefix=None,size=[15,4]):
	ax=stp.plot.bar(rot=0,figsize=size,) #ylim=[0.95,1.06],) #
	ax.set_ylabel("STP")
	ax.set_xlabel('')
	plt.grid(axis='y',linestyle=':')
	plt.tight_layout()
	plt.legend(bbox_to_anchor=(0.5, 1.1), loc='upper center',ncol=len(stp.columns))
	if file_prefix:
		plt.savefig("%s_STP.pdf" % file_prefix)
	
	ax=unf.plot.bar(rot=0,figsize=size,) #ylim=[0.6,1.2])
	ax.set_ylabel("Unfairness")
	plt.grid(axis='y',linestyle=':')
	ax.set_xlabel('')
	plt.tight_layout()
	plt.legend(bbox_to_anchor=(0.5, 1.1), loc='upper center',ncol=len(stp.columns))
	if file_prefix:
		plt.savefig("%s_UNF.pdf" % file_prefix)  


def select_clustering_solution(df,workload,algorithm=None):
	"""Given a Simlog Dataframe, generate subDataframe for the selected
	workload and algorithm."""
	if algorithm is None:
		return df[df["W#"]==workload]
	else:
		return df[(df["W#"]==workload) & (df["Algorithm"]==algorithm)]


#######################################################
### HIGH LEVEL API TO RUN SIMULATIONS FROM NOTEBOOK ###
#######################################################

def sim_read_application_info(summary_csv_path):
	"""Returns a dictionary of application objects from the file"""
	app_list=get_application_info_from_file(summary_csv_path)
	return { app.name:app for app in app_list }



def build_workload_from_str(workload_spec,info_apps,do_sort=False):
	"""
	Turn a list of application pattern names into a list of application objects

	The list can be optionally sorted alphabetically by application name
	"""
	workload=[]
	keys=[k for k in info_apps.keys()]
	for app_str in workload_spec:
		i=0
		## Parse tag
		tokens=app_str.split("/")
		while i<len(keys) and not (keys[i].startswith(tokens[0]) and (len(tokens)==1 or keys[i].endswith(tokens[1]))):
			i+=1
		if i<len(keys):
			app=info_apps[keys[i]]
			workload.append(App(app.name,app.properties))
		else:
			## Application not found
			return None

	if do_sort:
		workload.sort(key=lambda x: x.name.lower())
	return workload


## Auxiliary function to find an app pattern 
## and return an application object
def find_app(app_str,info_apps):
	tokens=app_str.split("/")
	keys=info_apps.keys()
	i=0
	while i<len(keys) and not (keys[i].startswith(tokens[0]) and (len(tokens)==1 or keys[i].endswith(tokens[1]))):
		i+=1
	if i<len(keys):
		app=info_apps[keys[i]]
		return App(app.name,app.properties)
	else:
		## Application not found
		return None 


def build_workload_from_clustering(cluster_spec,info_apps):
	"""
	Given a list of lists (with application pattern names)
	return the set of application objects in the cluster sorted
	in ascending order by name. The function also returns a list
	with the cluster ids each application belongs to

	The list can be optionally sorted alphabetically by application name
	"""	
	workload=[]
	id_alloc=0
	cluster_ids_unsorted=[]
	for pos_cluster,cluster in enumerate(cluster_spec):
		for app_str in cluster:
			app=find_app(app_str,info_apps)
			assert(app is not None)
			app.bench_id=id_alloc
			id_alloc+=1
			workload.append(app)
			cluster_ids_unsorted.append(pos_cluster)

	workload.sort(key=lambda x: x.name.lower()) 

	cluster_ids=[]
	for i,app in enumerate(workload):
		cluster_ids.append(cluster_ids_unsorted[app.bench_id])
		#New number
		app.bench_id=i

	return (workload,cluster_ids)
		
def sim_build_solution_df(sched_title,idx,sol_data,**kwargs):
	"""
	Turn a sol_data object (output from common clustering algorithms) into a dataframe
	"""
	kwargs.setdefault("workload_prefix","W")
	kwargs.setdefault("alg_suffix","")    
	kwargs.setdefault("max_bandwidth",float('Inf'))    

	workload_prefix=kwargs["workload_prefix"]
	alg_suffix=kwargs["alg_suffix"]
	max_bandwidth=kwargs["max_bandwidth"]

	metrics=compute_basic_metrics(sol_data,max_bandwidth)
	(sol_spec,statistics)=sol_data
	(workload,per_app_ways,per_app_masks,cluster_id)=sol_spec
	fix_intel_bug(per_app_ways,per_app_masks,cluster_id)


	## Determine slowdown with and widthout BW
	sbw=metrics["slowdown"]
	snbw=metrics["slowdown_no_bw"]

	workload_name="%s%i" % (workload_prefix,idx)
	alg_name="%s%s" % (sched_title,alg_suffix)

	columns=["W#","Algorithm","BenchID","Name","Cluster","Mask/NR_WAYS","SlowdownNB/ANTT","STP","Slowdown"]
	
	rows=[]
	for i,benchmark in enumerate(workload):	
		rows.append([workload_name,alg_name,i+1,benchmark.name,cluster_id[i],per_app_masks[i]+"("+str(per_app_ways[i])+")",round(snbw[i],5),round(1.0/sbw[i],5),round(sbw[i],5)])

	time_string= "%.6fs" % statistics["sim_time"]
	rows.append([workload_name,alg_name,len(workload)+1,"OVERALL","#"+str(statistics["total_branches"]),time_string,round(metrics["antt"],5),round(metrics["stp"],5),round(metrics["unfairness"],5)])

	return pd.DataFrame(rows,columns=columns)

def generate_clustering_chart_df(df,max_ways=11,**kwargs):
	"""
	Generate a TIKZ clustering chart from a solution dataframe
	"""
	algorithm=df["Algorithm"].values[0]
	workload=df["W#"].values[0]
	return generate_clustering_chart(df,workload,algorithm,max_ways=max_ways,**kwargs)


def apply_custom_partitioning(cluster_spec,ways_cluster,info_apps,nr_ways):
	"""
	Build a solution object out of the clustering proposed
	"""
	(workload,cluster_ids)=build_workload_from_clustering(cluster_spec,info_apps)

	per_clust_masks = get_partition_masks(ways_cluster)
	per_app_masks = [None]*len(cluster_ids)
	per_app_ways = [None]*len(cluster_ids)
	for i in range(len(cluster_ids)):
		per_app_ways[i] = ways_cluster[cluster_ids[i]]
		per_app_masks[i] = per_clust_masks[cluster_ids[i]]
	 
	patched_workload=[None]*len(per_app_ways)
	for cl_ix in range(len(cluster_ids)):
		# Get the corresponding indexes of this cluster to get the specific apps in the workload
		cl_apps_ixs = [i for i,clust_i in enumerate(cluster_ids) if clust_i==cl_ix]
		# Obtain the apps in the workload of this cluster
		clust_apps = [workload[app_ix] for app_ix in cl_apps_ixs]
		# Re build scaled apps to return the result and be able to obtain the performance metrics
		scaled_apps = get_scaled_properties_cluster(clust_apps, nr_ways)
		# Rebuild patched workload
		for i in range(len(cl_apps_ixs)):
			patched_workload[cl_apps_ixs[i]] = scaled_apps[i]

	stats={}
	stats["sim_time"]=0
	stats["total_branches"]=1
	solution_spec=(patched_workload,per_app_ways,per_app_masks,cluster_ids)

	return (solution_spec,stats) 


def apply_custom_partitioning_df(cluster_spec,ways_cluster,info_apps,nr_ways,**kwargs):
	"""
	Return a solution dataframe out of the clustering solution propsed
	"""
	sol=apply_custom_partitioning(cluster_spec,ways_cluster,info_apps,nr_ways)
	return sim_build_solution_df("custom",1,sol,**kwargs)


def apply_part_algorithm_df(algorithm,workload_spec,nr_ways,max_bandwidth=float('Inf'),parallel=False,debugging=False,uoptions={},workload_name="W"):
	"""
	Invoke the algorithm with PBBCache for the workload passed as a parameter and return a Pandas DataFrame with the results

	The workload can be specified in two ways:
	a) As a list of application objects
	b) As a tuple (list_of_str_app_patters , info_apps_dic)
		where info_apps_dic=sim_read_application_info("/Users/jcsaez/proyectos/pbbcache/simulator/data/volta2_turbo_off.300W.stalls.new.csv")
	"""
	if type(workload_spec) is list:
		workload=workload_spec
	else:
		(str_workload,info)=workload_spec
		workload=build_workload_from_str(str_workload,info,do_sort=False)

	sol=apply_part_algorithm(algorithm,workload,nr_ways,max_bandwidth,parallel,debugging,uoptions,workload_name)
	return sim_build_solution_df(algorithm,1,sol)