from openeye import oequacpac, oechem, oeomega, oedocking

def package_poses(poses, receptor_name, fout):
    ofs = oechem.oemolostream(fout)
    for mol, score in poses:
        oechem.OEWriteMolecule(ofs, mol)

def naive_pose_gen(smile, receptor_file):
    mol = oechem.OEGraphMol()
    oechem.OEParseSmiles(mol, smile)
    tautomers = EnumerateTautomers(mol)
    high_score_pose = DockConformers(tautomers, receptor_file)
    s = high_score_pose.GetEnergy()
    return high_score_pose, s

def GetSmilesFromFile(file_name):
    mol = oechem.OEGraphMol()
    ifs = oechem.oemolistream(file_name)
    while oechem.OEReadMolecule(ifs, mol):
        return oechem.OEMolToSmiles(mol)
    raise ValueError("OEReadMolecules must have not found any molecules in the file {}" % (file_name))

# Enumerate tautomers
def EnumerateTautomers(molecule):
    tautomer_options = oequacpac.OETautomerOptions()
    tautomer_options.SetMaxTautomersGenerated(4096)
    tautomer_options.SetMaxTautomersToReturn(16)
    tautomer_options.SetCarbonHybridization(True)
    tautomer_options.SetMaxZoneSize(50)
    tautomer_options.SetApplyWarts(True)
    pKa_norm = True
    tautomers = [ oechem.OEMol(tautomer) for tautomer in oequacpac.OEGetReasonableTautomers(molecule, tautomer_options, pKa_norm) ]
    return tautomers


def DockConformers(tautomers, receptor_filename, return_min_over_tautomers=True, dockMethod=None ):
    receptor = oechem.OEGraphMol()
    if not oedocking.OEReadReceptorFile(receptor, receptor_filename):
        oechem.OEThrow.Fatal("Unable to read receptor")

    if dockMethod == oedocking.OEDockMethod_Hybrid2 and (not oedocking.OEReceptorHasBoundLigand(receptor)):
        raise ValueError("Receptor does not have bound ligand")
    elif dockMethod is None and oedocking.OEReceptorHasBoundLigand(receptor):
        dockMethod = oedocking.OEDockMethod_Hybrid2
    elif dockMethod is None:
        dockMethod = oedocking.OEDockMethod_Chemgauss4

    #print('Initializing receptor...')
    dockResolution = oedocking.OESearchResolution_High
    dock = oedocking.OEDock(dockMethod, dockResolution)
    success = dock.Initialize(receptor)
    if success != oedocking.OEDockingReturnCode_Success:
        # raise ValueError("dock initalize failed")
        pass

    omegaOpts = oeomega.OEOmegaOptions(oeomega.OEOmegaSampling_Dense)
    omegaOpts.SetMaxConfs(800)
    omegaOpts.SetMaxSearchTime(60.0)  # time out
    omega = oeomega.OEOmega(omegaOpts)
    omega.SetStrictStereo(False)  # enumerate sterochemistry if uncertain

    # Dock tautomers
    docked_molecules = list()

    min_id = None
    min_s = None
    for id, mol in enumerate(tautomers):
        print("tayt done", id, mol)

        dockedMol = oechem.OEMol()

        # Expand conformers
        omega.Build(mol)

        # Dock molecule
        retCode = dock.DockMultiConformerMolecule(dockedMol, mol, 1)
        if (retCode != oedocking.OEDockingReturnCode_Success):
            print("Docking Failed with error code " + oedocking.OEDockingReturnCodeGetName(retCode))
            continue

        # Store docking data
        sdtag = oedocking.OEDockMethodGetName(dockMethod)
        oedocking.OESetSDScore(dockedMol, dock, sdtag)
        dock.AnnotatePose(dockedMol)
        s = dock.ScoreLigand(dockedMol)
        print(s)
        if min_s is None or s < min_s:
            min_s = s
            min_id = id

        docked_molecules.append( oechem.OEMol(dockedMol ))


    if return_min_over_tautomers:
        print(docked_molecules, min_id, docked_molecules[min_id])
        return docked_molecules[min_id]
    else:
        return docked_molecules