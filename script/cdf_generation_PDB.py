"""
@author: Arthur Garon
Department of Pharmaceutical Chemistry
University of Vienna
"""

import sys
import time

import CDPL.Base as Base
import CDPL.Biomol as Biomol
import CDPL.Chem as Chem
import CDPL.Math as Math
import argparse
import os


def parseArguments():
    parser = argparse.ArgumentParser(description='>>> Generate a cdf file (CDPKIT format) and a text file listing all residues from a PDB file.')
    parser.add_argument('-pdb',
                        dest='pdb',
                        required=True,
                        help='[Required] The path of the PDB file',
                        nargs=1)

    parser.add_argument('-n',
                        dest='name',
                        help='[Optional] The name of the cdf file (Default: name of the pdb file)',
                        nargs=1,
                        default=None)
    parser.add_argument('-o',
                        dest='output',
                        help='[Optional] The output folder where the cdf file will be generated (Default: current directory)',
                        nargs=1,
                        default=None)

    parse_args = parser.parse_args()
    return parse_args


def getResidueID(res):
    if Biomol.getChainID(res) != ' ':
        return Biomol.getResidueCode(res) + '_' + str(Biomol.getResidueSequenceNumber(res)) + '_' + Biomol.getChainID(res)

    return Biomol.getResidueCode(res) + '_' + str(Biomol.getResidueSequenceNumber(res))


def cdfMol_pdb(pdb, output, name):
    initial_time = time.time()
    cdf_mol = Chem.BasicMolecule()
    pdb_mol = Chem.BasicMolecule()

    pdb_str = open(pdb, 'r').read().replace('WAT', 'HOH').replace('HIE', 'HIS')
    pdb_reader = Biomol.PDBMoleculeReader(Base.StringIOStream(pdb_str))

    Biomol.setPDBApplyDictAtomBondingToNonStdResiduesParameter(pdb_reader, True)
    if not pdb_reader.read(pdb_mol):
        return None

    Chem.calcImplicitHydrogenCounts(pdb_mol, False)
    Chem.perceiveHybridizationStates(pdb_mol, False)
    Chem.setAtomSymbolsFromTypes(pdb_mol, False)
    Chem.perceiveSSSR(pdb_mol, False)
    Chem.setRingFlags(pdb_mol, False)
    Chem.setAromaticityFlags(pdb_mol, False)

    cdf_mol.assign(pdb_mol)
    for atom in cdf_mol.atoms:
        Chem.set3DCoordinatesArray(atom, Math.Vector3DArray())

    i = 0
    while i < cdf_mol.numAtoms:
        Chem.get3DCoordinatesArray(cdf_mol.getAtom(i)).addElement(Chem.get3DCoordinates(pdb_mol.getAtom(i)))
        i += 1

    tmp_output = output + name + ".cdf"
    try:
        Chem.FileCDFMolecularGraphWriter(tmp_output).write(cdf_mol)
    except:
        print('> Cdf_mol writing failure.')
        raise

    residues = Biomol.ResidueList(cdf_mol)
    tmp_output = output + name + "_residue_info.txt"
    with open(tmp_output, 'w') as txt_writer:
        txt_writer.write('residue name_resid_chain\n')
        for res in residues:
            res_id = getResidueID(res)
            txt_writer.write('{}: \n'.format(res_id))

    calc_time = time.time() - initial_time
    print('> Cdf and amino acid residue number list files generated in {}s'.format(int(calc_time)))


if __name__ == '__main__':
    args = parseArguments()

    pdb = args.pdb[0]

    if args.name is None:
        name = os.path.basename(pdb)[:-4]
    else:
        name = args.name[0]

    if args.output is None:
        output = './'
    else:
        output = args.output[0]

    cdfMol_pdb(pdb, output, name)



