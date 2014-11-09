#generic python modules
import argparse
import operator
from operator import itemgetter
import sys, os, shutil, itertools
import os.path

################################################################################################################################################
# RETRIEVE USER INPUTS
################################################################################################################################################

#=========================================================================================
# create parser
#=========================================================================================
version_nb = "0.0.1"
parser = argparse.ArgumentParser(prog='bilayer_undulations', usage='', add_help=False, formatter_class=argparse.RawDescriptionHelpFormatter, description=\
'''
**************************************************
v''' + version_nb + '''
author: Jean Helie (jean.helie@bioch.ox.ac.uk)
git: https://github.com/jhelie/bilayer_undulations
**************************************************
	
[ DESCRITPION ]

This script records and plots the undulations of a bilayer along a given axis. 

Note: In each frame, the coordinates of the beads of interest are centered around the
leaflets center of geometry before binning.
	
[ REQUIREMENTS ]

The following python modules are needed :
 - MDAnalysis
 - matplotlib

[ NOTES ]

1. It's a good idea to trjconv the xtc first and only outputs the relevant beads (e.g. PO4)
   as the script will run MUCH faster.	

2. You might want to use pbc_fix first if the bilayer crosses the simulation box boundaries.
   information is only calculated for transmembrane proteins.
   
3. Identification of the bilayer leaflets is controlled via 3 ways:
   (a) beads
    By default, the particles taken into account to define leaflet are:e
    -> name PO4 or name PO3 or name B1A
   
    Note that only lipids which contain one of the beads mentioned in the selection string
    will be taken into account. If you wish to specify your own selection string (e.g. to
    choose different beads or add a bead not in the default list in order to take into
    account a particular lipid specie) you can do so by supplying a file via the --beads
    option. This file should contain a single line that can be passed as the argument
    to MDAnalysis selectAtoms() routine and should not contain any quotation marks, e.g.:
     -> name PO4 or name PO3 or name B1A or name AM1
        
   (b) leaflet finding method
    By default leaflets are identified using the MDAnalysis LeafletFinder routine and the
    the optimum cutoff to identify 2 lipids groups is determined using the optimize_cutoff
    routine.
    This optimisation process can take time in large systems and you can specify your own
    cutoff value to skip this step. For instance to use the default 15 Angstrom cutoff
    directly (without optimising):
     -> '--leaflets 15'
   
    In very large systems (more then ~50,000 phospholipids) LeafletFinder (or rather the
    networkX module that it relies on) can fail. To  avoid this you can choose not to use
    this routine by specifying:
     -> '--leaflets large'
    In this case lipids whose headgroups z value is above the average lipids z value will
    be considered to make up the upper leaflet and those whose headgroups z value is below
    the average will be considered to be in the lower leaflet.
    This means that the bilayer should be as flat as possible in the gro file supplied in
    order to get a meaningful outcome.

    Note that if --leaflets is set to 'no', bilayer leaflets will not be identified meaning
    the script will not be able to discriminate between surfacic and interfacial proteins.

   (c) flipflopping lipids
    In case lipids flipflop during the trajectory, a file listing them can be supplied
    with the --flipflops option. Each line of this file should follow the format:
     -> 'resname,resid,starting_leaflet,z_bead'
    where starting_leaflet is either 'upper' or 'lower' - e.g. 'POPC,145,lower,PO4'. The
    z_bead is used to track the position of the lipid.
    However in most scenarios flip-flopping can probably be ignored without introducing
    any bias in the outcome. Just here for comprehensiveness.


[ USAGE ]

Option	      Default  	Description                    
-----------------------------------------------------
-f			: structure file [.gro] (required)
-x			: trajectory file [.xtc] (required)
-o			: name of output folder
-b			: beginning time (ns)
-e			: ending time (ns)	
-t 		[10]	: process every t-frames

Binning details
-----------------------------------------------------
--axis		[z]	: axis along which to bin the coords
--bins		[1]	: distance between bins (Angstrom)
--min	[-100]	: minimum of range (Angstrom)
--max	[100]	: maxomum of range (Angstrom)

Lipids identification (see note 2)
-----------------------------------------------------
--bead		[PO4]	: lipids bead name
--flipflops		: input file with flipflopping lipids
--leaflets	optimise	: leaflet identification ('optimise', 'large' or float)

Graph details
-----------------------------------------------------
--xticks	[10]	: nb of ticks along x axis
--yticks	[10]	: nb of ticks along y axis
 
Other options
-----------------------------------------------------
--version		: show version number and exit
-h, --help		: show this menu and exit
  
''')

#data options
parser.add_argument('-f', nargs=1, dest='grofilename', default=['no'], help=argparse.SUPPRESS, required=True)
parser.add_argument('-x', nargs=1, dest='xtcfilename', default=['no'], help=argparse.SUPPRESS, required=True)
parser.add_argument('-o', nargs=1, dest='output_folder', default=['no'], help=argparse.SUPPRESS)
parser.add_argument('-b', nargs=1, dest='t_start', default=[-1], type=int, help=argparse.SUPPRESS)
parser.add_argument('-e', nargs=1, dest='t_end', default=[-1], type=int, help=argparse.SUPPRESS)
parser.add_argument('-t', nargs=1, dest='frames_dt', default=[10], type=int, help=argparse.SUPPRESS)

#binning details
parser.add_argument('--bins', nargs=1, dest='dbins', default=[1], type=float, help=argparse.SUPPRESS)
parser.add_argument('--axis', dest='axis', choices=['x','y','z'], default=['z'], help=argparse.SUPPRESS)
parser.add_argument('--min', nargs=-100, dest='bmin', default=[-100], type=float, help=argparse.SUPPRESS)
parser.add_argument('--max', nargs=100, dest='bmax', default=[100], type=float, help=argparse.SUPPRESS)

#lipids identification
parser.add_argument('--bead', nargs=1, dest='beadname', default=['PO4'], help=argparse.SUPPRESS)
parser.add_argument('--flipflops', nargs=1, dest='selection_file_ff', default=['no'], help=argparse.SUPPRESS)
parser.add_argument('--leaflets', nargs=1, dest='cutoff_leaflet', default=['optimise'], help=argparse.SUPPRESS)

#graph details
parser.add_argument('--xticks', nargs=1, dest='xticks', default=[10], type=int, help=argparse.SUPPRESS)
parser.add_argument('--yticks', nargs=1, dest='yticks', default=[10], type=int, help=argparse.SUPPRESS)

#other options
parser.add_argument('--version', action='version', version='%(prog)s v' + version_nb, help=argparse.SUPPRESS)
parser.add_argument('-h','--help', action='help', help=argparse.SUPPRESS)

#=========================================================================================
# store inputs
#=========================================================================================

#parse user inputs
#-----------------
args = parser.parse_args()
#data options
args.grofilename = args.grofilename[0]
args.xtcfilename = args.xtcfilename[0]
args.output_folder = args.output_folder[0]
args.t_start=args.t_start[0]
args.t_end=args.t_end[0]
args.frames_dt=args.frames_dt[0]
#binning details
args.dbins=args.dbins[0]
args.axis=args.axis[0]
args.bmin=args.bmin[0]
args.bmax=args.bmax[0]
#lipids identification
args.beadname = args.beadname[0]
args.cutoff_leaflet = args.cutoff_leaflet[0]
args.selection_file_ff = args.selection_file_ff[0]
#graph details
args.xticks=args.xticks[0]
args.yticks=args.yticks[0]
#leaflet identification
if args.cutoff_leaflet != "large" and args.cutoff_leaflet != "optimise":
	try:
		args.cutoff_leaflet = float(args.cutoff_leaflet)
	except:
		print "Error: the argument of the --leaflets option should set to a number, 'optimise' or 'large', see note 2"
		sys.exit(1)


#=========================================================================================
# import modules (doing it now otherwise might crash before we can display the help menu!)
#=========================================================================================
#generic science modules
try:
	import math
except:
	print "Error: you need to install the maths module."
	sys.exit(1)
try:
	import numpy as np
except:
	print "Error: you need to install the np module."
	sys.exit(1)
try:
	import scipy
except:
	print "Error: you need to install the scipy module."
	sys.exit(1)
try:
	import matplotlib as mpl
	mpl.use('Agg')
	import matplotlib.colors as mcolors
	mcolorconv = mcolors.ColorConverter()
	import matplotlib.cm as cm				#colours library
	import matplotlib.ticker
	from matplotlib.ticker import MaxNLocator
	from matplotlib.font_manager import FontProperties
	fontP=FontProperties()
except:
	print "Error: you need to install the matplotlib module."
	sys.exit(1)
try:
	import pylab as plt
except:
	print "Error: you need to install the pylab module."
	sys.exit(1)

#MDAnalysis module
try:
	import MDAnalysis
	from MDAnalysis import *
	import MDAnalysis.analysis
	import MDAnalysis.analysis.leaflet
	import MDAnalysis.analysis.distances
	#set MDAnalysis to use periodic boundary conditions
	MDAnalysis.core.flags['use_periodic_selections'] = True
	MDAnalysis.core.flags['use_KDTree_routines'] = False
except:
	print "Error: you need to install the MDAnalysis module first. See http://mdanalysis.googlecode.com"
	sys.exit(1)

#=========================================================================================
# sanity check
#=========================================================================================
if not os.path.isfile(args.grofilename):
	print "Error: file " + str(args.grofilename) + " not found."
	sys.exit(1)
if not os.path.isfile(args.xtcfilename):
	print "Error: file " + str(args.xtcfilename) + " not found."
	sys.exit(1)
if args.selection_file_ff != "no" and not os.path.isfile(args.selection_file_ff):
	print "Error: file " + str(args.selection_file_ff) + " not found."
	sys.exit(1)
if args.t_end != -1 and args.t_end < args.t_start:
	print "Error: the starting time (" + str(args.t_start) + "ns) for analysis is later than the ending time (" + str(args.t_end) + "ns)."
	sys.exit(1)
if args.dbins < 0:
	print "Error: the distance between bins must be greater than zero."
	sys.exit(1)

#=========================================================================================
# create folders and log file
#=========================================================================================
if args.output_folder == "no":
	args.output_folder="bilayer_undulation_" + args.xtcfilename[:-4]
if os.path.isdir(args.output_folder):
	print "Error: folder " + str(args.output_folder) + " already exists, choose a different output name via -o."
	sys.exit(1)
else:
	#create folder
	os.mkdir(args.output_folder)

	#create log
	filename_log=os.getcwd() + '/' + str(args.output_folder) + '/cluster_prot.log'
	output_log=open(filename_log, 'w')		
	output_log.write("[bilayer_undulations v" + str(version_nb) + "]\n")
	output_log.write("\nThis folder and its content were created using the following command:\n\n")
	tmp_log="python bilayer_undulations.py"
	for c in sys.argv[1:]:
		tmp_log+=" " + c
	output_log.write(tmp_log + "\n")
	output_log.close()

	#copy input files
	if args.selection_file_ff != "no":
		shutil.copy2(args.selection_file_ff,args.output_folder + "/")	

##########################################################################################
# FUNCTIONS DEFINITIONS
##########################################################################################

#=========================================================================================
# data loading
#=========================================================================================

def set_lipids_beads():

	global leaflet_sele_string

	#set default beads
	leaflet_sele_string = "name " + str(args.beadname)

	return
def load_MDA_universe():
	
	global U
	global all_atoms
	global nb_atoms
	global nb_frames_xtc
	global frames_to_process
	global frames_to_write
	global nb_frames_to_process
	global f_start
	global f_end	
	f_start = 0
	
	print "\nLoading trajectory..."
	U = Universe(args.grofilename, args.xtcfilename)
	U_timestep = U.trajectory.dt
	all_atoms = U.selectAtoms("all")
	nb_atoms = all_atoms.numberOfAtoms()
	nb_frames_xtc = U.trajectory.numframes		
	#sanity check
	if U.trajectory[nb_frames_xtc-1].time/float(1000) < args.t_start:
		print "Error: the trajectory duration (" + str(U.trajectory.time/float(1000)) + "ns) is shorted than the starting stime specified (" + str(args.t_start) + "ns)."
		sys.exit(1)
	if U.trajectory.numframes < args.frames_dt:
		print "Warning: the trajectory contains fewer frames (" + str(nb_frames_xtc) + ") than the frame step specified (" + str(args.frames_dt) + ")."

	#rewind traj (very important to make sure that later the 1st frame of the xtc will be used for leaflet identification)
	U.trajectory.rewind()
	
	#create list of index of frames to process
	if args.t_end != -1:
		f_end = int((args.t_end*1000 - U.trajectory[0].time) / float(U_timestep))
		if f_end < 0:
			print "Error: the starting time specified is before the beginning of the xtc."
			sys.exit(1)
	else:
		f_end = nb_frames_xtc - 1		
	if args.t_start != -1:
		f_start = int((args.t_start*1000 - U.trajectory[0].time) / float(U_timestep))
		if f_start > f_end:
			print "Error: the starting time specified is after the end of the xtc."
			sys.exit(1)
	if (f_end - f_start)%args.frames_dt == 0:
		tmp_offset = 0
	else:
		tmp_offset = 1
	frames_to_process = map(lambda f:f_start + args.frames_dt*f, range(0,(f_end - f_start)//args.frames_dt+tmp_offset))
	nb_frames_to_process = len(frames_to_process)
			
	#check for the presence of lipids
	test_beads = U.selectAtoms(leaflet_sele_string)
	if test_beads.numberOfAtoms() == 0:
		print "Error: invalid selection string '" + str(leaflet_sele_string) + "'"
		print "-> no lipid particles selected, check the --leaflets and/or --bead options."
		sys.exit(1)

	return
def identify_ff():
	print "\nReading selection file for flipflopping lipids..."
	
	#declare variables
	global lipids_ff_nb
	global lipids_ff_info
	global lipids_ff_resnames
	global lipids_ff_leaflet
	global lipids_ff_u2l_index
	global lipids_ff_l2u_index
	global lipids_sele_ff
	global lipids_sele_ff_bead
	global lipids_sele_ff_bonds
	global lipids_sele_ff_VMD_string
	global leaflet_sele_string
	lipids_ff_nb = 0
	lipids_ff_info = {}
	lipids_ff_resnames = []
	lipids_ff_leaflet = []
	lipids_ff_u2l_index = []
	lipids_ff_l2u_index = []
	lipids_sele_ff = {}
	lipids_sele_ff_bead = {}
	lipids_sele_ff_bonds = {}
	lipids_sele_ff_VMD_string={}
		
	with open(args.selection_file_ff) as f:
		lines = f.readlines()
	lipids_ff_nb = len(lines)
	print " -found " + str(lipids_ff_nb) + " flipflopping lipids"
	leaflet_sele_string = leaflet_sele_string + " and not ("
	for l_index in range(0,lipids_ff_nb):
		line = lines[l_index]
		if line[-1] == "\n":
			line = line[:-1]
		try:
			line_content = line.split(',')
			if len(line_content) != 4:
				print "Error: wrong format for line " + str(l_index+1) + " in " + str(args.selection_file_ff) + ", see note 4 in bilayer_perturbations --help."
				print " ->", line
				sys.exit(1)
			#read current lipid details
			lip_resname = line_content[0]
			lip_resnum = int(line_content[1])
			lip_leaflet = line_content[2]
			lip_bead = line_content[3]
			lipids_ff_info[l_index] = [lip_resname,lip_resnum,lip_leaflet,lip_bead]
						
			#update: starting leaflets
			if lip_leaflet not in lipids_ff_leaflet:
				lipids_ff_leaflet.append(lip_leaflet)

			#update: index in directional lists
			if lip_leaflet == "upper":
				lipids_ff_u2l_index.append(l_index)
			elif lip_leaflet == "lower":
				lipids_ff_l2u_index.append(l_index)
			else:
				print "->unknown starting leaflet '" + str(lip_leaflet) + "'."
				sys.exit(1)
			
			#update: resnames
			if lip_resname not in lipids_ff_resnames:
				lipids_ff_resnames.append(lip_resname)
	
			#update: leaflet selection string
			if l_index==0:
				leaflet_sele_string+="(resname " + str(lip_resname) + " and resnum " + str(lip_resnum) + ")"
			else:
				leaflet_sele_string+=" or (resname " + str(lip_resname) + " and resnum " + str(lip_resnum) + ")"

			#create selections
			lipids_sele_ff[l_index] = U.selectAtoms("resname " + str(lip_resname) + " and resnum " + str(lip_resnum))
			lipids_sele_ff_bead[l_index] = lipids_sele_ff[l_index].selectAtoms("name " + str(lip_bead))
			lipids_sele_ff_VMD_string[l_index]="resname " + str(lipids_ff_info[l_index][0]) + " and resid " + str(lipids_ff_info[l_index][1])
			if lipids_sele_ff[l_index].numberOfAtoms() == 0:
				print "Error:"
				print line
				print "-> no such lipid found."
				sys.exit(1)	
		except:
			print "Error: invalid flipflopping lipid selection string on line " + str(l_index+1) + ": '" + line + "'"
			sys.exit(1)
	leaflet_sele_string+=")"		

	return
def identify_leaflets():
	print "\nIdentifying leaflets..."
	
	#declare variables
	global leaflet_sele
	global leaflet_sele_atoms
	leaflet_sele = {}
	leaflet_sele_atoms = {}
	for l in ["lower","upper","both"]:
		leaflet_sele[l] = {}
		leaflet_sele_atoms[l] = {}
	
	#check the leaflet selection string is valid
	test_beads = U.selectAtoms(leaflet_sele_string)
	if test_beads.numberOfAtoms() == 0:
		print "Error: invalid selection string '" + str(leaflet_sele_string) + "'"
		print "-> no particles selected."
		sys.exit(1)

	#use LeafletFinder:
	if args.cutoff_leaflet != 'large':
		if args.cutoff_leaflet == 'optimise':
			print " -optimising cutoff..."
			cutoff_value = MDAnalysis.analysis.leaflet.optimize_cutoff(U, leaflet_sele_string)
			L = MDAnalysis.analysis.leaflet.LeafletFinder(U, leaflet_sele_string, cutoff_value[0])
		else:
			L = MDAnalysis.analysis.leaflet.LeafletFinder(U, leaflet_sele_string, args.cutoff_leaflet)
	
		if np.shape(L.groups())[0]<2:
			print "Error: imposssible to identify 2 leaflets."
			sys.exit(1)
		if L.group(0).centerOfGeometry()[2] > L.group(1).centerOfGeometry()[2]:
			leaflet_sele["upper"] = L.group(0)
			leaflet_sele["lower"] = L.group(1)
		else:
			leaflet_sele["upper"] = L.group(1)
			leaflet_sele["lower"] = L.group(0)
		leaflet_sele["both"] = leaflet_sele["lower"] + leaflet_sele["upper"]
		if np.shape(L.groups())[0] == 2:
			print " -found 2 leaflets: ", leaflet_sele["upper"].numberOfResidues(), '(upper) and ', leaflet_sele["lower"].numberOfResidues(), '(lower) lipids'
		else:
			other_lipids=0
			for g in range(2, np.shape(L.groups())[0]):
				other_lipids += L.group(g).numberOfResidues()
			print " -found " + str(np.shape(L.groups())[0]) + " groups: " + str(leaflet_sele["upper"].numberOfResidues()) + "(upper), " + str(leaflet_sele["lower"].numberOfResidues()) + "(lower) and " + str(other_lipids) + " (others) lipids respectively"
	#use cog:
	else:
		leaflet_sele["both"] = U.selectAtoms(leaflet_sele_string)
		tmp_lipids_avg_z = leaflet_sele["both"].centerOfGeometry()[2]
		leaflet_sele["upper"] = leaflet_sele["both"].selectAtoms("prop z > " + str(tmp_lipids_avg_z))
		leaflet_sele["lower"] = leaflet_sele["both"].selectAtoms("prop z < " + str(tmp_lipids_avg_z))
		print " -found 2 leaflets: ", leaflet_sele["upper"].numberOfResidues(), '(upper) and ', leaflet_sele["lower"].numberOfResidues(), '(lower) lipids'
		
	return

#=========================================================================================
# outputs
#=========================================================================================

def post_process_data():
	
	global bins_pos_avg
	
	#normalise
	z_hist["upper"] /= float(nb_frames_to_process)
	z_hist["lower"] /= float(nb_frames_to_process)
	
	return
def write_xvg_undulations():
	
	global bins_pos_avg
	
	filename_xvg = os.getcwd() + '/' + str(args.output_folder) + '/bilayer_undulations.xvg'
	output_xvg = open(filename_xvg, 'w')
	output_xvg.write("@ title \"Distribution of leaflet coodinates along the " + str(args.axis) + " axis\"\n")
	output_xvg.write("@ xaxis  label \"time (ns)\"\n")
	output_xvg.write("@ yaxis  label \"frequency\"\n")
	output_xvg.write("@ autoscale ONREAD xaxes\n")
	output_xvg.write("@ TYPE XY\n")
	output_xvg.write("@ view 0.15, 0.15, 0.95, 0.85\n")
	output_xvg.write("@ legend on\n")
	output_xvg.write("@ legend box on\n")
	output_xvg.write("@ legend loctype view\n")
	output_xvg.write("@ legend 0.98, 0.8\n")
	output_xvg.write("@ legend length 2\n")
	output_xvg.write("@ s0 legend \"upper\"\n")
	output_xvg.write("@ s1 legend \"lower\"\n")
	for b_index in range(0,bins_nb):
		results = str(bins_pos_avg[b_index]) + "	" + str(z_hist["upper"][b_index]) + "	" + str(z_hist["lower"][b_index])
		output_xvg.write(results + "\n")
	output_xvg.close()
	return
def graph_xvg_undulations():
	#create filenames
	#----------------
	filename_png = os.getcwd() + '/' + str(args.output_folder) + '/bilayer_undulations.png'
	filename_svg = os.getcwd() + '/' + str(args.output_folder) + '/bilayer_undulations.svg'

	#create figure
	#-------------
	fig=plt.figure(figsize=(8, 6.2))
	fig.suptitle("Distribution of leaflet coordinates along the " + str(args.axis) + " axis")

	#plot data: %
	#------------
	ax = fig.add_subplot(111)
	plt.plot(bins_pos_avg, z_hist["upper"], color = colours["upper"], linewidth = 1, alpha = 0.22, label = "upper")
	plt.plot(bins_pos_avg, z_hist["lower"], color = colours["lower"], linewidth = 1, alpha = 0.22, label = "lower")
	plt.fill_between(bins_pos_avg, np.zeros(bins_nb), z_hist["upper"], color = colours["upper"], edgecolor = colours["upper"], linewidth = 0, alpha = 0.22)
	plt.fill_between(bins_pos_avg, np.zeros(bins_nb), z_hist["lower"], color = colours["lower"], edgecolor = colours["lower"], linewidth = 0, alpha = 0.22)
	fontP.set_size("small")
	ax.legend(prop=fontP)
	plt.xlabel('z position ($\AA$)', fontsize="small")
	plt.ylabel('frequency (avg number of atoms)', fontsize="small")

	#save figure
	#-----------
	ax.set_xlim(args.bmin, args.bmax)
	ax.set_ylim(0, max(max(z_hist["lower"]),max(z_hist["upper"])))
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	ax.xaxis.set_ticks_position('bottom')
	ax.yaxis.set_ticks_position('left')
	ax.xaxis.set_major_locator(MaxNLocator(nbins=args.xticks))
	ax.yaxis.set_major_locator(MaxNLocator(integer=True, nbins=args.yticks))
	plt.setp(ax.xaxis.get_majorticklabels(), fontsize="small" )
	plt.setp(ax.yaxis.get_majorticklabels(), fontsize="small" )
	fig.savefig(filename_png)
	fig.savefig(filename_svg)
	plt.close()
	return

##########################################################################################
# ALGORITHM
##########################################################################################

#=========================================================================================
# pre-processing
#=========================================================================================

#declare variables
global z_hist
global colours
global bins_nb
global bins_boundaries
global bins_pos_avg
global axis_index
z_hist = {}
colours = {}
colours["upper"] = "b"
colours["lower"] = "r"

#process inputs
set_lipids_beads()
load_MDA_universe()
if args.selection_file_ff != "no":
	identify_ff()
identify_leaflets()

#update variables
if args.axis == "x":
	axis_index = 0
elif args.axis == "y":
	axis_index = 1
else:
	axis_index = 2
tmp_range = (args.bmax-args.bmin)
bins_nb = int(tmp_range/ float(args.dbins))
bins_boundaries = np.linspace(args.bmin,args.bmax,bins_nb+1)
bins_pos_avg = np.zeros(bins_nb)
for b_index in range(0,bins_nb):
	bins_pos_avg[b_index] = args.bmin + (b_index+0.5)*float(args.dbins)
print " -range of bins:", str(tmp_range), "Ansgtrom"
print " -nb of bins:", str(bins_nb)

#initiliase leaflet containers
z_hist["upper"] = np.zeros(bins_nb)
z_hist["lower"] = np.zeros(bins_nb)

#=========================================================================================
# generate data
#=========================================================================================
print "\nBrowsing trajectory..."

for f_index in range(0,nb_frames_to_process):
	
	ts = U.trajectory[frames_to_process[f_index]]

	#display progress
	progress = '\r -processing frame ' + str(f_index+1) + '/' + str(nb_frames_to_process) + ' (every ' + str(args.frames_dt) + ' frame(s) from frame ' + str(f_start) + ' to frame ' + str(f_end) + ' out of ' + str(nb_frames_xtc) + ')      '  
	sys.stdout.flush()
	sys.stdout.write(progress)

	#center coordinates
	tmp_upper = leaflet_sele["upper"].coordinates()[:,axis_index]
	tmp_lower = leaflet_sele["lower"].coordinates()[:,axis_index]
	tmp_mid = (np.average(tmp_upper) + np.average(tmp_lower)) / float(2)
	tmp_upper -= tmp_mid
	tmp_lower -= tmp_mid
	
	#bin them
	z_hist["upper"] += np.histogram(tmp_upper,bins_boundaries)[0]
	z_hist["lower"] += np.histogram(tmp_lower,bins_boundaries)[0]
	
print ""

#=========================================================================================
# post-processing
#=========================================================================================
post_process_data()

#=========================================================================================
# produce outputs
#=========================================================================================
print "\nWriting outputs..."
write_xvg_undulations()
graph_xvg_undulations()
		
#=========================================================================================
# exit
#=========================================================================================
print "\nFinished successfully! Check output in ./" + args.output_folder + "/"
print ""
sys.exit(0)
