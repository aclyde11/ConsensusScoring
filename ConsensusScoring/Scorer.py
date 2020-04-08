from openeye import oechem
from openeye import oedocking

def PrintScore(score, pose):
    print(score.SystematicSolidBodyOptimize(pose))
    print("Total ligand score =  %f" % score.ScoreLigand(pose))
    print("Score components contributions to score:")
    for comp in score.GetComponentNames():
        print("%15s: %6.2f" % (comp, score.ScoreLigandComponent(pose, comp)))

def PrintAtomScore(score, pose, atom):
    print("\nAtom: %d  score: %6.2f " % (atom.GetIdx(), score.ScoreAtom(atom, pose)))
    print("Score components contribution to atom scores:")
    for comp in score.GetComponentNames():
        print("%15s: %.2f" % (comp, score.ScoreAtomComponent(atom, pose, comp)))

