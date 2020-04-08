from openeye import oequacpac, oechem, oeomega, oedocking

def GetSmilesFromFile(file_name):
    mol = oechem.OEGraphMol()
    ifs = oechem.oemolistream(file_name)
    while oechem.OEReadMolecule(ifs, mol):
        return oechem.OEMolToSmiles(mol)
    raise ValueError("OEReadMolecules must have not found any molecules in the file {}" % (file_name))
