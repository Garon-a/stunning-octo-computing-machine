"""
@author: Arthur Garon
Department of Pharmaceutical Chemistry
University of Vienna
"""

import sys
import time

import CDPL.Biomol as Biomol
import CDPL.Chem as Chem
import CDPL.Pharm as Pharm
import argparse
import os
import pickle


def parseArguments():
    parser = argparse.ArgumentParser(description='>>> Generate the global trajectory object file from the cdf file.')
    parser.add_argument('-cdf',
                        dest='cdf',
                        required=True,
                        help='[Required] The path of the cdf file',
                        nargs=1)
    parser.add_argument('-lig',
                        dest='ligand_code',
                        required=True,
                        help='[Required] The 3-letters code of the ligand',
                        nargs=1)

    parser.add_argument('-n',
                        dest='gt_name',
                        help='[Optional] The name of the gt object file (Default: name of the cdf file)',
                        nargs=1,
                        default=None)
    parser.add_argument('-o',
                        dest='output',
                        help='[Optional] The output folder where the gt object file will be generated (Default: current directory)',
                        nargs=1,
                        default=None)


    parse_args = parser.parse_args()
    return parse_args


def loadCDFMolecule(fname):
    mol = Chem.BasicMolecule()
    cdf_reader = Chem.FileCDFMoleculeReader(fname)

    if not cdf_reader.read(mol):
        return None

    Chem.calcImplicitHydrogenCounts(mol, False)
    Chem.perceiveHybridizationStates(mol, False)
    Chem.setAtomSymbolsFromTypes(mol, False)
    Chem.perceiveSSSR(mol, False)
    Chem.setRingFlags(mol, False)
    Chem.setAromaticityFlags(mol, False)

    return mol


class Ph4InteractionInformation(object):
    def __init__(self, lig_feature, env_feature):
        ftype_names = {Pharm.FeatureType.H_BOND_ACCEPTOR: 'HBA', Pharm.FeatureType.H_BOND_DONOR: 'HBD',
                       Pharm.FeatureType.POS_IONIZABLE: 'PI', Pharm.FeatureType.NEG_IONIZABLE: 'NI',
                       Pharm.FeatureType.AROMATIC: 'AR', Pharm.FeatureType.HYDROPHOBIC: 'H',
                       Pharm.FeatureType.X_VOLUME: 'XV'}

        lig_feature_type = ftype_names[Pharm.getType(lig_feature)]
        lig_residue_code = Biomol.getResidueCode(Pharm.getSubstructure(lig_feature).atoms[0])
        lig_residue_number = Biomol.getResidueSequenceNumber(Pharm.getSubstructure(lig_feature).atoms[0])
        lig_residue_chain = Biomol.getChainID(Pharm.getSubstructure(lig_feature).atoms[0])

        env_feature_type = ftype_names[Pharm.getType(env_feature)]
        env_residue_code = Biomol.getResidueCode(Pharm.getSubstructure(env_feature).atoms[0])
        env_residue_number = Biomol.getResidueSequenceNumber(Pharm.getSubstructure(env_feature).atoms[0])
        env_residue_chain = Biomol.getChainID(Pharm.getSubstructure(env_feature).atoms[0])

        self.interaction_type = '{}-{}'.format(lig_feature_type, env_feature_type)
        self.lig_residue = '{}_{}_{}'.format(lig_residue_code, lig_residue_number, lig_residue_chain)
        self.env_residue = '{}_{}_{}'.format(env_residue_code, env_residue_number, env_residue_chain)

        atoms = []
        for atom in Pharm.getSubstructure(lig_feature).atoms:
            key_atom = '{}:{}'.format(Chem.getSymbol(atom), Biomol.getSerialNumber(atom))
            atoms.append(key_atom)

        self.lig_atom = sorted(atoms, key=lambda k: int(k.split(':')[1]))

        atoms = []
        for atom in Pharm.getSubstructure(env_feature).atoms:
            key_atom = '{}:{}'.format(Chem.getSymbol(atom), Biomol.getSerialNumber(atom))
            atoms.append(key_atom)

        self.env_atom = sorted(atoms, key=lambda k: int(k.split(':')[1]))

    def __str__(self):
        txt = ''
        for key in sorted(self.__dict__):
            txt += '> {:<20}: {:<20}\n'.format(key, self.__dict__[key])
        return txt

    def __eq__(self, other):
        return self.__dict__ == other.__dict__

    def __ne__(self, other):
        return not self == other

    def getInteractionType(self):
        return self.interaction_type

    def getLigand(self):
        return [self.lig_residue, str(self.lig_atom)] if self.lig_atom is not None else [self.lig_residue]

    def getEnvironment(self):
        return [self.env_residue, str(self.env_atom)] if self.env_atom is not None else [self.env_residue]


def getPh4Interactions(lig_pharm, interactions):
        features = {}
        for lig_feature in [x for x in lig_pharm if interactions.getValues(x) != []]:
            for env_feature in interactions.getValues(lig_feature):
                inf = Ph4InteractionInformation(lig_feature, env_feature)
                if inf.interaction_type not in features.keys():
                    features[inf.interaction_type] = []
                features[inf.interaction_type].append(inf)
        return features


def getPh4InteractionDictionary(cdf_path, ligand_code):

    ph4_interaction_dictionary = {}
    cdf_mol = loadCDFMolecule(cdf_path)
    num_confs = Chem.getNumConformations(cdf_mol)

    ligand = Chem.Fragment()
    for atom in cdf_mol.atoms:
        if Biomol.getResidueCode(atom) == ligand_code:
            Biomol.extractResidueSubstructure(atom, cdf_mol, ligand, False)
            break

    if ligand.numAtoms == 0:
        print('> Could not find ligand {}'.format(ligand_code))
        return 0

    Chem.perceiveSSSR(ligand, True)
    lig_env = Chem.Fragment()

    lig_pharm = Pharm.BasicPharmacophore()
    env_pharm = Pharm.BasicPharmacophore()
    pharm_gen = Pharm.DefaultPharmacophoreGenerator(True)

    analyzer = Pharm.DefaultInteractionAnalyzer()
    interactions = Pharm.FeatureMapping()

    for y in range(num_confs):
        lig_pharm.clear()
        env_pharm.clear()
        interactions.clear()
        lig_env.clear()

        coords_func = Chem.AtomConformer3DCoordinatesFunctor(y)
        pharm_gen.setAtom3DCoordinatesFunction(coords_func)
        Biomol.extractEnvironmentResidues(ligand, cdf_mol, lig_env, coords_func, 7)
        Chem.perceiveSSSR(lig_env, True)
        pharm_gen.generate(ligand, lig_pharm)

        pharm_gen.generate(lig_env, env_pharm)
        analyzer.analyze(lig_pharm, env_pharm, interactions)
        ph4_interaction_dictionary[y] = getPh4Interactions(lig_pharm, interactions)

    return ph4_interaction_dictionary


if __name__ == '__main__':
    args = parseArguments()

    initial_time = time.time()
    cdf = args.cdf[0]
    ligand_code = args.ligand_code[0]

    if args.gt_name is None:
        gt_name = os.path.basename(cdf).split('_chunk_')[0][:-4]
    else:
        gt_name = args.gt_name[0]

    if args.output is None:
        output = './'
    else:
        output = args.output[0]

    ph4_interaction_dictionary = getPh4InteractionDictionary(cdf, ligand_code)

    with open(output + gt_name + '.gt', 'wb') as handle:
        pickle.dump(ph4_interaction_dictionary, handle)

    calc_time = time.time() - initial_time
    print('> Gt object file generated in {}s'.format(int(calc_time)))




