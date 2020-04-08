from openeye import oechem
from openeye import oedocking

from ConsensusScoring import Scorer
import subprocess

subprocess.call("test1.py", shell=True)

protein = oechem.OEGraphMol()
oedocking.OEReadReceptorFile(protein, "/Users/austin/Downloads/ADRP.oeb.gz")
score = oedocking.OEScore(oedocking.OEScoreType_Chemgauss4)
score.Initialize(protein)
print(score.IsInitialized())
pose = oechem.OEMol()
imstr = oechem.oemolistream("/Users/austin/Downloads/adp_h.pdb")


for ligand in imstr.GetOEMols():
    Scorer.PrintScore(score, ligand)
    for atom in ligand.GetAtoms():
        Scorer.PrintAtomScore(score, ligand, atom)
