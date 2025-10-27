import numpy as np
import shutil

import MDAnalysis as mda 
from MDAnalysis.analysis import pca, align
from MDAnalysis.coordinates.memory import MemoryReader
from MDAnalysis.analysis.base import AnalysisFromFunction
import tqdm
import uuid
import os
from sklearn.decomposition import PCA
import os 

import uuid

print('imported modules') 

import argparse 
#TODO. do more things in memory

import warnings
warnings.filterwarnings('ignore')

import pdb
import rlcompleter
pdb.Pdb.complete=rlcompleter.Completer(locals()).complete

parser = argparse.ArgumentParser()

parser.add_argument('--topology', '-t', dest='topology' , help='topology file',type=str) 
parser.add_argument('--trajectory', '-j', dest='traj', help='trajectory file',type=str) 
parser.add_argument('--out', '-o', dest='out_file',default="pca", help='Output file',type=str) 
parser.add_argument('--stride', '-dt', dest='stride' , default=1, help='Stride through trajectory skipping this many frames.',type=int) 
parser.add_argument('--no-std', dest='standardise_bool', action="store_false", help='Choose whether or not to standardise coordinates. Default behaviour is to standardise the cordinates.') 
parser.add_argument('--save-proj', dest='save_projection', action="store_true", help='Choose whether or not to standardise coordinates. Default behaviour is to standardise the cordinates.',default=False) 
parser.add_argument('--n-components', '-n', dest='n_components' , default = 2, help='calculate this many Principal Components of the trajectory.',type=int) 
parser.add_argument('--vis-multiply', dest='extend' , default = 1, help='multiply the eigenvectors by this magnitude for visualisation purposes.',type=float) 
parser.add_argument('--dir-root', dest='dir_root', help='Determines where the output will go. Default is to put it all into a temporary directory which is deleted after the run.') 

parser.add_argument('--ref', dest='reference',  help='Reference structure, used for alignment') 
parser.add_argument('--ref_top', dest='ref_top',  help='Reference topology, used for alignment') 
parser.add_argument('--symmetry-list', dest='symmetry_list',  help='List of selections which compose a single symmetry group. This script expects symmetric groups to have the same number of atoms when the selection string is applied. For example you might have 4 identical chains so you would pass the argument "--symmetry-list \'segid A\' \'segid B\' \'segid C\' \'segid D\'" to this script. The user is also warned that the groups will be treated cyclically. So in the previous example A will move to B, B to C, C to D and D to A. This will keep the relative arrangement of components so long as the user names the groups in the correct order. ',nargs="+") 

parser.add_argument('--ref_sel', dest='ref_sel', help='Selection string for fitting. Will be applied to both target and reference structures.',type=str, default="name CA") 
parser.add_argument('--analysis_selections', '-s', dest='analysis_selections', help='Selection string for fitting. Will be applied to both target and reference structures.',type=str, required=True, nargs="+") 
parser.add_argument('--n-frames', '-fn', dest='n_vis_frames', help='Number of frames to visualise',type=int, default=30) 
parser.add_argument('--in-mem', dest='in_mem', action="store_true", help='Do alignment processing in memory. I advise against doing this if the trajectory is large.') 

args = parser.parse_args()

print(args.standardise_bool)
selection_string = args.ref_sel 

universe = mda.Universe(args.topology,args.traj)
original_trajectory_n_frames = universe.trajectory.n_frames
#universe.trajectory = universe.trajectory[0::args.stride] 


selection_object = universe.select_atoms('all')

if args.dir_root == None : 
    dir_root = '/tmp/' + str(uuid.uuid4())  + '/'
else:
    dir_root = args.dir_root + '/'

print(dir_root)
os.makedirs(dir_root,exist_ok=True)

def copy_coords(ag):
        return ag.positions.copy()

selection_object.write(dir_root + 'temp_pdb_file.pdb')

analysis_universe = mda.Universe(dir_root + 'temp_pdb_file.pdb')
coordinates = np.array([selection_object.positions for ts in universe.trajectory[::args.stride]])
     
analysis_universe = analysis_universe.load_new (coordinates)

with mda.Writer((dir_root + 'analysis_universe_traj.dcd'), analysis_universe.atoms.n_atoms) as W:
    for ts in analysis_universe.trajectory:
        W.write(analysis_universe.atoms)

if args.reference != None: 
    if args.ref_top != None:
        ref_universe = mda.Universe(args.ref_top, args.reference)
    else:
        ref_universe = mda.Universe(args.topology, args.reference)
else:
    #ref_universe = analysis_universe.copy()
    ref_universe = mda.Universe(dir_root + 'temp_pdb_file.pdb')
    ref_universe.load_new(coordinates)



analysis_universe = mda.Universe(dir_root + 'temp_pdb_file.pdb')
#load like this to over write first frame
analysis_universe.load_new(dir_root + 'analysis_universe_traj.dcd')


#aligns by default
analysis_universe.atoms.write (dir_root + 'analysis.pdb')
print('Aligning Trajectory')

if args.in_mem==True:
    aligner = align.AlignTraj(analysis_universe, ref_universe, filename=dir_root + 'aligned.dcd',select=args.selection_string,verbose=True).run()
    analysis_universe = analysis_universe.load_new(dir_root + 'aligned.dcd')
    aligner = align.AlignTraj(analysis_universe, analysis_universe, filename=dir_root + 'aligned2.dcd',select = args.selection_string,verbose=True).run()
    os.remove(dir_root + 'aligned.dcd')
    analysis_universe = analysis_universe.load_new(dir_root + 'aligned2.dcd')
else:
    #need to fix this
    aligner = align.AlignTraj(analysis_universe, ref_universe, filename=dir_root + 'aligned.dcd',select = args.ref_sel,verbose=True).run()

    analysis_universe = analysis_universe.load_new(dir_root + 'aligned.dcd')
    aligner = align.AlignTraj(analysis_universe, analysis_universe, filename=dir_root + 'aligned2.dcd',select = args.ref_sel,verbose=True).run()
    analysis_universe = analysis_universe.load_new(dir_root + 'aligned2.dcd')

    #aligned_coords = AnalysisFromFunction(copy_coords, aligner_universe).run().results

print('Trajectory Aligned')

analysis_universe_coordinates = np.array([analysis_universe.atoms.positions for ts in analysis_universe.trajectory],dtype=np.float32)
analysis_universe_coordinates = analysis_universe_coordinates.reshape((analysis_universe.trajectory.n_frames,analysis_universe.atoms.n_atoms*3))

distances = np.zeros([analysis_universe.trajectory.n_frames,len(args.analysis_selections)])
reference_universe_selection = ref_universe.select_atoms(args.ref_sel)

for i in tqdm.tqdm(range(analysis_universe.trajectory.n_frames)):
    analysis_universe.trajectory[i]
    for j in range(len(args.analysis_selections)):
        subunits_distance_atoms = analysis_universe.select_atoms(args.analysis_selections[j])
        origin = reference_universe_selection.atoms.center_of_mass()

        #compare along the z axis.
        norming = (subunits_distance_atoms.center_of_mass()[2]-origin[2])
        distances[i,j] = norming

#data_matrix = np.hstack([], distances]
header_str = ''
for i in range(0,len(args.analysis_selections)):
    header_str = header_str + '\"' + str(args.analysis_selections[i]) +  '\"'
    header_str = header_str [:-2]

print(distances)
np.savetxt(args.out_file, distances, header=header_str)
