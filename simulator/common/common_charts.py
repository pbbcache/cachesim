# -*- coding: utf-8 -*-

#
# common_charts.py
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
import matplotlib.pyplot as plt
from matplotlib.markers import *

# Class to represent a Labeled point to use in creating the output charts of the simulator
class LabeledPoint:
	def __init__(self, label, x, y, z=0):
		self.label = label
		self.x = x
		self.y = y
		self.z = z

	def __repr__(self):
		return repr((self.label, self.x, self.y, self.z))



# Function to build an array of Labeled points to use in creating the charts
def buildArrays(labeledPoints):
	X=[]
	Y=[]
	Z=[]
	Labels=[]
	for point in labeledPoints:
		Labels.append(point.label)
		X.append(point.x)
		Y.append(point.y)
		Z.append(point.z)

	return (Labels,X,Y,Z)	

# Function that creates and shows a 2D scatter chart from labeled points of different algorithms
def scatter2D(labeledPointsSched,markerSpecs,xlabel,ylabel,windowTitle="Scatter plot 2D"):
	plt.rcParams["figure.figsize"] = [13.0, 9.0]
	fig = plt.gcf()
	fig.canvas.set_window_title(windowTitle)
	plt.xlabel(xlabel)
	plt.ylabel(ylabel)
	#plt.xlim((-0.05, 1.05))
	#plt.ylim((-0.05, 1.05))
	
	#markers = ['o', '^', 'D', '*', '+', 's', '1' ]
	

	for sched in labeledPointsSched.keys():
		schedPoints=labeledPointsSched[sched]
		(Labels,X,Y,Z)=buildArrays(schedPoints)

		## Hack para reconocer el color
		markerSpec=markerSpecs[sched]
		if 	len(markerSpec)>1:
			extraArgs=dict(facecolor=markerSpec[1])
		else:
			extraArgs=dict()
		plt.scatter(X, Y, s=70, facecolors='none', edgecolors='black', marker=markerSpec[0], label=sched, **extraArgs)
		for x, y, note in zip(X,Y,Labels):
			plt.annotate(
				note,
				xy = (x, y), xytext = (12, 5),
				textcoords = 'offset points', ha = 'right', va = 'bottom', color='blue')
	
	plt.legend(loc=0, scatterpoints=1)
	plt.grid(True)
	plt.show()

# Generic function to create and/or save a 2D scatter chart to compare different algorithms
def scatter2DGen(labeledPointsSched,markerSpecs,xlabel,ylabel,mode=0,figSize=[9.0,9.0],filename=None,windowTitle=None,legendLocation='best',axes_labelsize=None,show_window=False):
	plt.rcParams["figure.figsize"]=figSize

	plt.figure()
	fig = plt.gcf()
	if windowTitle!=None:
		fig.canvas.set_window_title(windowTitle)

	if axes_labelsize!=None:
		plt.xlabel(xlabel,fontsize=axes_labelsize)
		plt.ylabel(ylabel,fontsize=axes_labelsize)
	else:
		plt.xlabel(xlabel)
		plt.ylabel(ylabel)

	for sched in labeledPointsSched.keys():
		schedPoints=labeledPointsSched[sched]
		(Labels,X,Y,Z)=buildArrays(schedPoints)

		## Hack para reconocer el color
		markerSpec=markerSpecs[sched]
		if 	len(markerSpec)>1:
			extraArgs=dict(facecolor=markerSpec[1])
		else:
			extraArgs=dict()
		if mode==0:
			plt.scatter(X, Y, s=70, facecolors='none', edgecolors='black', marker=markerSpec[0], label=sched, **extraArgs)
			## X vs Y
			for x, y, note in zip(X,Y,Labels):
				plt.annotate(
					note,
					xy = (x, y), xytext = (12, 5),
					textcoords = 'offset points', ha = 'right', va = 'bottom', color='blue')
		elif mode==1:
			## X vs Z
			plt.scatter(X, Z, s=70, facecolors='none', edgecolors='black', marker=markerSpec[0], label=sched, **extraArgs)
			for x, y, note in zip(X,Z,Labels):
				plt.annotate(
					note,
					xy = (x, y), xytext = (12, 5),
					textcoords = 'offset points', ha = 'right', va = 'bottom', color='blue')
		else:
			print "Unsupported mode"
			return 1
	plt.legend(loc=legendLocation, scatterpoints=1)
	plt.grid(True)
	plt.draw()
	if filename!=None:
		plt.savefig(filename)

	if show_window:
		plt.show()
	return 0
