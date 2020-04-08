from ConsensusScoring.scoring.utils import naive_pose_gen, package_poses
import pandas as pd
import numpy as np
from tqdm import tqdm
import argparse

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', help='tab seperated file with SMILES NAME where with no header', required=True, type=str)
    parser.add_argument('-o', help='output of bundle of poses to be posted on leader board', required=True, type=str)
    parser.add_argument('--receptor_file', required=True, type=str, help='.oeb or .oeb.gz file for the protein and pocket target you are submitting for.')
    parser.add_argument('-n', required=False, default=10, type=float, help='max number of top scoring compounds to output.')
    return parser.parse_args()

#######
## How to take a ML hit on a SMILES to an actual 3D Pose
##
## the function naive_pose_gen generates 3D conformers and various sterochemistries
## and docks them to the single provided receptor. This is not the full consensus scoring score,
## but it will provide you the 3D pose to submit to the CS pipeline. The CS pipeline only accepts poses, not SMILES.
######

if __name__ == '__main__':
    args = get_args()
    df = pd.read_csv(args.i, header=None, sep='\t', names=['smiles', 'name'])

    my_poses = []
    for iter, smile in enumerate(tqdm(df.smiles)):
        highest_scoring_pose, oe_dock_score = naive_pose_gen(smile, args.receptor_file)
        print(highest_scoring_pose)
        my_poses.append((highest_scoring_pose, oe_dock_score))

        if iter == 10:
            break

    scores_ranking = np.argsort(np.array(zip(*my_poses)[1]))
    my_poses = [my_poses[i] for i in scores_ranking]
    my_poses = my_poses[:min(len(my_poses), args.n)]
    print(*list(map(lambda x : x[0].GetTitle(), my_poses)))





