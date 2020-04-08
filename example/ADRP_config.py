from openeye import oechem, oedocking
from pymol import cmd
from ConsensusScoring import Target, Scorer
import subprocess

INPUT_PDB_CODE = '6w02'


## OE Score Prepare
cmd.fetch(INPUT_PDB_CODE)
cmd.do("save {}".format("example/6w02.pdb"))

cmd = "python /Users/austin/Downloads/proteinprep.py -ligandname APR -waterfilter nothing -waterprocessing ignore -verbose -proteinfilter chainA -in example/6w02.pdb -ligout example/6w02_lig.pdb -cplxout example/6w02_apo.pdb"
subprocess.call(cmd, shell=True)
cmd = "python /Users/austin/Downloads/proteinprep.py -ligandname APR -waterfilter nothing -waterprocessing ignore -verbose -proteinfilter chainA -in example/6w02.pdb -cplxout example/6w02.pdb"
subprocess.call(cmd, shell=True)


Target.prepare_receptor("example/6w02.pdb", "example/")


protein = oechem.OEGraphMol()
oedocking.OEReadReceptorFile(protein, "example/ADRP.oeb.gz")
score = oedocking.OEScore(oedocking.OEScoreType_Chemgauss4)
score.Initialize(protein)
print(score.IsInitialized())
pose = oechem.OEMol()
imstr = oechem.oemolistream("example/6w02_lig.sdf")

for ligand in imstr.GetOEMols():
    Scorer.PrintScore(score, ligand)
    for atom in ligand.GetAtoms():
        Scorer.PrintAtomScore(score, ligand, atom)