import os
from openeye import oespruce, oechem

header = ""

def read_pdb_file(pdb_file):
    print(f'Reading receptor from {pdb_file}...')

    from openeye import oechem
    ifs = oechem.oemolistream()
    ifs.SetFlavor(oechem.OEFormat_PDB, oechem.OEIFlavor_PDB_Default | oechem.OEIFlavor_PDB_DATA | oechem.OEIFlavor_PDB_ALTLOC)  # noqa

    if not ifs.open(pdb_file):
        oechem.OEThrow.Fatal("Unable to open %s for reading." % pdb_file)

    mol = oechem.OEGraphMol()
    if not oechem.OEReadMolecule(ifs, mol):
        oechem.OEThrow.Fatal("Unable to read molecule from %s." % pdb_file)
    ifs.close()

    return (mol)

def prepare_receptor(complex_pdb_filename, output_basepath, dimer=False):
    """
    Parameters
    ----------
    complex_pdb_filename : str
        The complex PDB file to read in
    output_basepath : str
        Base path for output
    dimer : bool, optional, default=False
        If True, append a header to generate the dimer as a biological unit
    """
    basepath, filename = os.path.split(complex_pdb_filename)
    prefix, extension = os.path.splitext(filename)
    prefix = os.path.join(output_basepath, prefix)


    complex = read_pdb_file(complex_pdb_filename)

    print('Identifying design units...')
    design_units = list(oespruce.OEMakeDesignUnits(complex))
    if len(design_units) == 1:
        design_unit = design_units[0]
    elif len(design_units) > 1:
        print('More than one design unit found---using first one')
        design_unit = design_units[0]
    elif len(design_units) == 0:
        raise Exception('No design units found')

    # Prepare the receptor
    print('Preparing receptor...')
    from openeye import oedocking
    protein = oechem.OEGraphMol()
    design_unit.GetProtein(protein)
    ligand = oechem.OEGraphMol()
    design_unit.GetLigand(ligand)
    receptor = oechem.OEGraphMol()
    oedocking.OEMakeReceptor(receptor, protein, ligand)
    oedocking.OEWriteReceptorFile(receptor, basepath + "/ADRP.oeb.gz")

    with oechem.oemolostream(f'{prefix}-protein.pdb') as ofs:
        oechem.OEWriteMolecule(ofs, protein)
    with oechem.oemolostream(f'{prefix}-ligand.mol2') as ofs:
        oechem.OEWriteMolecule(ofs, ligand)
    with oechem.oemolostream(f'{prefix}-ligand.pdb') as ofs:
        oechem.OEWriteMolecule(ofs, ligand)
    with oechem.oemolostream(f'{prefix}-ligand.sdf') as ofs:
        oechem.OEWriteMolecule(ofs, ligand)

    # Filter out UNK from PDB files (which have covalent adducts)
    pdbfile_lines = [ line for line in open(f'{prefix}-protein.pdb', 'r') if 'UNK' not in line ]
    with open(f'{prefix}-protein.pdb', 'wt') as outfile:
        outfile.write(''.join(pdbfile_lines))
