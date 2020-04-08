from ConsensusScoring.utils import converter
from openeye import oechem



class Ligand:

    def __init__(self, name, pose_file, source, **data):
        pose = oechem.OEMol()

        imstr = oechem.oemolistream(pose_file)
        oechem.OEReadMolecule(imstr, pose)
        imstr.close()

        self.__name = name
        self.__smiles = converter.GetSmilesFromFile(pose)
        self.__source = source

        self.__pose = pose
        self.__CS_score = None

        self.__data = dict(data)

    def get_name(self):
        return self.__name

    def get_smiles(self):
        return self.__smiles

    def get_source(self):
        return self.__source

    def get_pose(self):
        return self.__pose

    def get_CS_score(self):
        return self.__CS_score

    def get_data(self, key = None):
        if key is None:
            return self.__data
        else:
            return self.__data[key]

    def set_data(self, key, value):
        self.__data[key] = value

