"""
@author: Arthur Garon, Thomas Seidel
Department of Pharmaceutical Sciences
University of Vienna
"""

import time
import math
import os
import sys
import copy

from collections import defaultdict

import MDAnalysis

import matplotlib
import matplotlib.pylab as plt
import numpy as np
import pandas as pd
import scipy.stats

from mpl_toolkits.axes_grid1 import make_axes_locatable

import CDPL.Base as Base
import CDPL.Biomol as Biomol
import CDPL.Chem as Chem
import CDPL.Math as Math
import CDPL.Pharm as Pharm
import CDPL.MolProp as MolProp
import CDPL.Vis as Vis

xrange = range

ALH_HYDROPHOBICIY_THRESH = 0.15


def cdfMol(psf, dcd, output, name, chunk_size=500):
    initial_time = time.time()
    u = MDAnalysis.Universe(psf, dcd)
    chunks_start = [0]
    chunks_end = []

    if chunk_size == 0 or chunk_size > len(u.trajectory):
        chunk_size = len(u.trajectory)

    for i in range(len(u.trajectory)):
        if i % chunk_size == 0 and i != 0 and i != (len(u.trajectory) - 1):
            chunks_end.append(i)
            chunks_start.append(i)
        elif i == (len(u.trajectory) - 1):
            chunks_end.append(i + 1)

    cdf_mol = Chem.BasicMolecule()
    waters = {}
    
    if psf.endswith('.pdb'):
        pdb_str = open(psf, 'r').read().replace('WAT', 'HOH').replace('HIE', 'HIS').replace('SPC X', 'HOH X')
        pdb_reader = Biomol.PDBMoleculeReader(Base.StringIOStream(pdb_str))

        Biomol.setPDBApplyDictAtomBondingToNonStdResiduesParameter(pdb_reader, True)

        try:
            if not pdb_reader.read(cdf_mol):
                sys.exit('!!! Could not read PDB file: ' + psf)
        except:
            sys.exit('!!! Could not read PDB file: ' + psf)

        for atom in cdf_mol.atoms:              
            Chem.setImplicitHydrogenCount(atom, 0)
            
            if Biomol.getResidueCode(atom) == 'HOH':
                res_id = Biomol.getResidueSequenceNumber(atom)
                 
                if res_id in waters:
                    waters[res_id].append(int(res_id))
                else:
                    waters[res_id] = [int(res_id)]
                    
            array = Math.Vector3DArray()
            array.resize(chunk_size, Math.Vector3D())

            Chem.set3DCoordinatesArray(atom, array)
    else:
        cdf_mol.reserveMemoryForAtoms(len(u.atoms))
        cdf_mol.reserveMemoryForBonds(len(u.bonds))

        # construct atoms
        for md_atom in u.atoms:
            atom = cdf_mol.addAtom()

            if md_atom.type.lower().startswith('cl') :
                Chem.setSymbol(atom, 'Cl')
            elif md_atom.type.lower().startswith('br') :
                Chem.setSymbol(atom, 'Br')
            else:
                Chem.setSymbol(atom, md_atom.name[0])

            Biomol.setChainID(atom, md_atom.segid)

            if md_atom.resname == 'WAT' or md_atom.resname == 'TIP3':
                Biomol.setResidueCode(atom, 'HOH')
            elif md_atom.resname == 'HSD' or md_atom.resname == 'HIE':
                Biomol.setResidueCode(atom, 'HIS')
            else:
                Biomol.setResidueCode(atom, md_atom.resname)

            Biomol.setResidueSequenceNumber(atom, int(md_atom.resid))
            Biomol.setResidueAtomName(atom, md_atom.name)
            Biomol.setSerialNumber(atom, int(md_atom.id))

            if Biomol.getResidueCode(atom) == 'HOH':
                if md_atom.resid in waters:
                    waters[md_atom.resid].append(int(md_atom.id))
                else:
                    waters[md_atom.resid] = [int(md_atom.id)]

            # fix positive charge on arginin nitrogen
            if md_atom.resname == 'ARG' and md_atom.name == 'NH2':
                Chem.setFormalCharge(atom, 1)

            Chem.set3DCoordinates(cdf_mol.getAtom(int(md_atom.id)), md_atom.position.tolist())
            coords_array = Math.Vector3DArray()
            coords_array.resize(chunk_size, Math.Vector3D())

            Chem.set3DCoordinatesArray(atom, coords_array)

        Chem.setAtomTypesFromSymbols(cdf_mol, True)

        # construct bonds
        for md_bond in u.bonds:
            if not Chem.getType(cdf_mol.atoms[int(md_bond.atoms[0].index)]) == Chem.AtomType.H == Chem.getType(
                    cdf_mol.atoms[int(md_bond.atoms[1].index)]):
                cdf_mol.addBond(int(md_bond.atoms[0].index), int(md_bond.atoms[1].index))
   
        # make sane biomolecule
        for a in cdf_mol.atoms:
            Chem.setImplicitHydrogenCount(a, 0)

    for water in waters.values():
        if len(water) < 2:
            continue

        for atom_idx in water:
            if Chem.getType(cdf_mol.atoms[atom_idx]) == Chem.AtomType.O:
                if cdf_mol.atoms[atom_idx].numBonds > 1:
                    break

                for atom_idx2 in water:
                    if Chem.getType(cdf_mol.atoms[atom_idx2]) == Chem.AtomType.H:
                        cdf_mol.addBond(atom_idx, atom_idx2)

                break

    Chem.perceiveSSSR(cdf_mol, True)
    Chem.setRingFlags(cdf_mol, True)
    Chem.perceiveBondOrders(cdf_mol, True)
    Chem.perceiveHybridizationStates(cdf_mol, True)
    Chem.setAromaticityFlags(cdf_mol, True)
    Chem.calcFormalCharges(cdf_mol, True)

    print('> Cdf atoms and bonds setup in {}s'.format(int(time.time() - initial_time)))

    # read timesteps
    tmp_time = time.time()
    u = MDAnalysis.Universe(psf, dcd)

    for i in xrange(len(chunks_start)):
        tmp_array = defaultdict(lambda: [])
        for ts in u.trajectory[chunks_start[i]: chunks_end[i]]:
            for md_atom in u.atoms:
                tmp_array[int(md_atom.index)].append([x for x in md_atom.position])

        for index in list(tmp_array.keys()):
            coords_array = Chem.get3DCoordinatesArray(cdf_mol.getAtom(index))
            for element_id, element in enumerate(tmp_array[index]):
                for dim in xrange(3):
                    coords_array[element_id][dim] = float(element[dim])

        if (i+1)*chunk_size > len(u.trajectory):
            over_size = (i+1)*chunk_size - len(u.trajectory)

            for md_atom in cdf_mol.atoms:
                coords_array = Chem.get3DCoordinatesArray(md_atom)
                for y in range(over_size):
                    coords_array.popLastElement()

        tmp_output = os.path.join(output, name + '_chunk_' + str(i) + '.cdf')

        try:
            if not Chem.FileCDFMolecularGraphWriter(tmp_output).write(cdf_mol):
                sys.exit('!!! Could not write CDF file: ' + tmp_output)
        except:
            sys.exit('!!! Could not write CDF file: ' + tmp_output)

        print('> Cdf chunk {} generated in {}s'.format(i, int(time.time() - tmp_time)))
        tmp_time = time.time()

    calc_time = time.time() - initial_time
    print('> Cdf file generated in {}s'.format(int(calc_time)))
    return len(chunks_start)

def mergeCDFMolecule(fname, chunk_number):
    initial_time = time.time()
    main_cdf = loadCDFMolecule(fname)
    
    for chunk_index in range(1, int(chunk_number)):
        tmp_cdf = loadCDFMolecule(os.path.join(os.path.dirname(fname), os.path.basename(fname).split('_chunk_')[0] + '_chunk_' + str(chunk_index) + '.cdf'))
        coords_func = Chem.AtomConformer3DCoordinatesFunctor(0)
        for atom_index in range(tmp_cdf.numAtoms):
            main_coords_array = Chem.get3DCoordinatesArray(main_cdf.getAtom(atom_index))
            tmp_coords_array = Chem.get3DCoordinatesArray(tmp_cdf.getAtom(atom_index))
            for position_index in range(tmp_coords_array.size):
                main_coords_array.addElement(tmp_coords_array[position_index])

    main_file = fname.split('_chunk_')[0] + '.cdf'
    cdf_writer = Chem.FileCDFMolecularGraphWriter(main_file)

    try:
        if not cdf_writer.write(main_cdf):
            sys.exit('!!! Could not write merged CDF file: ' + main_file)
    except:
        sys.exit('!!! Could not write merged CDF file: ' + main_file)
        
    for chunk_index in range(chunk_number):
        os.remove(os.path.join(os.path.dirname(fname), os.path.basename(fname).split('_chunk_')[0] + '_chunk_' + str(chunk_index) + '.cdf'))

    calc_time = time.time() - initial_time
    print('> {} Cdf files merged in {}s'.format(chunk_number, int(calc_time)))

def loadCDFMolecule(fname):
    mol = Chem.BasicMolecule()
    cdf_reader = Chem.FileCDFMoleculeReader(fname)

    try:
        if not cdf_reader.read(mol):
            sys.exit('!!! Could not load CDF file: ' + fname)
    except:
        sys.exit('!!! Could not load CDF file: ' + fname)
            
    Chem.calcImplicitHydrogenCounts(mol, False)
    Chem.perceiveHybridizationStates(mol, False)
    Chem.setAtomSymbolsFromTypes(mol, False)
    Chem.perceiveSSSR(mol, False)
    Chem.setRingFlags(mol, False)
    Chem.setAromaticityFlags(mol, False)

    return mol

def cdfMol_pdb(pdb, output, name):
    initial_time = time.time()
    pdb_mol = Chem.BasicMolecule()

    pdb_str = open(pdb, 'r').read().replace('WAT', 'HOH').replace('HIE', 'HIS').replace('SPC X', 'HOH X')
    pdb_reader = Biomol.PDBMoleculeReader(Base.StringIOStream(pdb_str))

    Biomol.setPDBApplyDictAtomBondingToNonStdResiduesParameter(pdb_reader, True)

    try:
        if not pdb_reader.read(pdb_mol):
            sys.exit('!!! Could not load PDB file: ' + pdb)
    except:
        sys.exit('!!! Could not load PDB file: ' + pdb)
            
    Pharm.prepareForPharmacophoreGeneration(pdb_mol, False)
    Biomol.setHydrogenResidueSequenceInfo(pdb_mol, False)
    
    for atom in pdb_mol.atoms:
        array = Math.Vector3DArray()
        array.addElement(Chem.get3DCoordinates(atom))
        Chem.set3DCoordinatesArray(atom, array)
        
    tmp_output = os.path.join(output, name + ".cdf")
    
    try:
        if not Chem.FileCDFMolecularGraphWriter(tmp_output).write(pdb_mol):
            sys.exit('!!! Could not write CDF file: ' + tmp_output)
    except:
        sys.exit('!!! Could not write CDF file: ' + tmp_output)

    calc_time = time.time() - initial_time
    
    print('> Cdf file generated in {}s'.format(int(calc_time)))

def drawLigand(cdf_path, lig_code, output):
    cdf_mol = loadCDFMolecule(cdf_path)
    ligand = Chem.Fragment()

    for atom in cdf_mol.atoms:
        if Biomol.getResidueCode(atom) == lig_code:
            Biomol.extractResidueSubstructure(atom, cdf_mol, ligand, True)
            break

    for atom in ligand.atoms:
        if Biomol.hasSerialNumber(atom):
            serial = Biomol.getSerialNumber(atom)
        else:
            serial =  Biomol.getSerialNumber(atom.atoms[0])
                
        Chem.setAtomMappingID(atom, serial)

    Chem.calcTopologicalDistanceMatrix(ligand, True)
    Chem.perceiveComponents(ligand, True)
    Chem.perceiveSSSR(ligand, True)
    Chem.calc2DCoordinates(ligand, True)
    
    svg_writer = Vis.FileSVGMolecularGraphWriter(output + lig_code + '.svg')

    Vis.setAtomColorTableParameter(svg_writer, Vis.AtomColorTable.ELEMENT_COLORS_2D)
    
    try:
        if not svg_writer.write(ligand):
            sys.exit('!!! Could not generate ligand structure depiction')
    except:
        sys.exit('!!! Could not generate ligand structure depiction')
    
    svg_writer.close()

class Ph4InteractionInformation(object):
    
    def __init__(self, lig_feature, env_feature):
        ftype_names = {Pharm.FeatureType.H_BOND_ACCEPTOR: 'HBA', Pharm.FeatureType.H_BOND_DONOR: 'HBD',
                       Pharm.FeatureType.POSITIVE_IONIZABLE: 'PI', Pharm.FeatureType.NEGATIVE_IONIZABLE: 'NI',
                       Pharm.FeatureType.AROMATIC: 'AR', Pharm.FeatureType.HYDROPHOBIC: 'H',
                       Pharm.FeatureType.HALOGEN_BOND_ACCEPTOR: 'XBA', Pharm.FeatureType.HALOGEN_BOND_DONOR: 'XBD',
                       Pharm.FeatureType.EXCLUSION_VOLUME: 'XV'}

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
            if Biomol.hasSerialNumber(atom):
                serial = Biomol.getSerialNumber(atom)
            else:
                serial =  Biomol.getSerialNumber(atom.atoms[0])
                
            key_atom = '{}:{}'.format(Chem.getSymbol(atom), serial)
            atoms.append(key_atom)

        self.lig_atom = sorted(atoms, key=lambda k: int(k.split(':')[1]))

        atoms = []

        for atom in Pharm.getSubstructure(env_feature).atoms:
            if Biomol.hasSerialNumber(atom):
                serial = Biomol.getSerialNumber(atom)
            else:
                serial = Biomol.getSerialNumber(atom.atoms[0])
                
            key_atom = '{}:{}'.format(Chem.getSymbol(atom), serial)
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

def getPh4InteractionDictionary(cdf_path, ligand_code, alh):
    ph4_interaction_dictionary = {}
    cdf_mol = loadCDFMolecule(cdf_path)
    num_confs = Chem.getNumConformations(cdf_mol)
    ligand = Chem.Fragment()

    for atom in cdf_mol.atoms:
        if Biomol.getResidueCode(atom) == ligand_code:
            Biomol.extractResidueSubstructure(atom, cdf_mol, ligand, False)
            break

    if ligand.numAtoms == 0:
        sys.exit('!!! Could not find ligand {}'.format(ligand_code))

    Chem.perceiveSSSR(ligand, True)
    MolProp.calcAtomHydrophobicities(ligand, True)

    lig_env = Chem.Fragment()
    lig_pharm = Pharm.BasicPharmacophore()
    env_pharm = Pharm.BasicPharmacophore()
    pharm_gen = Pharm.DefaultPharmacophoreGenerator(Pharm.DefaultPharmacophoreGenerator.STATIC_H_DONORS)

    pharm_gen.enableFeature(Pharm.FeatureType.HALOGEN_BOND_ACCEPTOR, True);

    if alh:
        h_gen = Pharm.HydrophobicAtomFeatureGenerator()
        h_gen.setHydrophobicityThreshold(ALH_HYDROPHOBICIY_THRESH)

        pharm_gen.setFeatureGenerator(Pharm.FeatureType.HYDROPHOBIC, h_gen)
    
    analyzer = Pharm.DefaultInteractionAnalyzer()
    interactions = Pharm.FeatureMapping()

    for y in range(num_confs):
        lig_pharm.clear()
        env_pharm.clear()
        interactions.clear()

        coords_func = Chem.AtomConformer3DCoordinatesFunctor(y)

        if y == 0:
            Biomol.extractEnvironmentResidues(ligand, cdf_mol, lig_env, coords_func, 7)
            Chem.perceiveSSSR(lig_env, True)
            MolProp.calcAtomHydrophobicities(lig_env, True, alh)
            
        pharm_gen.setAtom3DCoordinatesFunction(coords_func)
        pharm_gen.generate(ligand, lig_pharm)
        pharm_gen.generate(lig_env, env_pharm)

        analyzer.analyze(lig_pharm, env_pharm, interactions)
        ph4_interaction_dictionary[y] = getPh4Interactions(lig_pharm, interactions)

        if (y + 1) % 500 == 0:
            print('... Processed %s frames ...' % str(y + 1))
        
    return ph4_interaction_dictionary

def getGlobalPh4InteractionList(ph4_interaction_dictionary):
    global_ph4_interaction_list = []
    global_count = []

    for dict_int in ph4_interaction_dictionary.values():
        for list_int in dict_int.values():
            for tmp_int in list_int:
                if tmp_int not in global_ph4_interaction_list:
                    global_ph4_interaction_list.append(tmp_int)
                    global_count.append(1)
                else:
                    global_count[global_ph4_interaction_list.index(tmp_int)] += 1

    tmp_list = [(global_ph4_interaction_list[x], global_count[x]) for x in range(len(global_ph4_interaction_list))]

    return sorted(tmp_list, key=lambda k: (k[0].interaction_type, int(k[0].env_residue.split('_')[1])))


def getDataframeIM(global_ph4_interaction_list):
    col = sorted(
        list({tuple([x[0].getInteractionType()] + x[0].getLigand()) for x in global_ph4_interaction_list}),
        key=lambda k: (k[0], k[1]))
    ind = sorted(list({tuple(x[0].getEnvironment()) for x in global_ph4_interaction_list}),
                 key=lambda k: int(k[0].split('_')[1]))
    matrix = np.zeros((len(ind), len(col)))
    
    for x in global_ph4_interaction_list:
        matrix[ind.index(tuple(x[0].getEnvironment()))][col.index(tuple([x[0].getInteractionType()] + x[0].getLigand()))] = x[1]

    return pd.DataFrame(matrix, index=ind, columns=col)

def plotInteractionMap(df, number_frames, output=None):
    SMALL_FONT_SIZE = 18
    LARGE_FONT_SIZE = 25
    
    plt.clf()
    plt.cla()
    plt.close()

    ftype_colors = {'AR-AR':'blue', 'AR-PI':'lightblue', 'PI-AR':'purple', 'PI-NI':'azure', 'NI-PI':'pink',
                    'H-H':'y', 'HBA-HBD':'red', 'HBD-HBA': 'green', 'XBD-XBA': 'magenta', 'XBA-XBD': 'magenta'}

    fig, ax = plt.subplots()
    heatmap = ax.pcolor(df, cmap=plt.cm.Blues, alpha=0.9, vmax=number_frames)
    ind_label = []
    labels = [x for x in df.columns]
    stock = labels[0][0]

    for i, v in enumerate(labels):
        if v[0] == stock:
            continue
        else:
            ind_label.append(i)
            stock = v[0]

    ind_label.append(0)
    yl, yh = ax.get_ylim()
    left = yl - (yh - yl) * -0.000001
    right = yh + (yh - yl) * -0.000001

    for x in ind_label:
        ax.vlines(x, left, right, color=ftype_colors[labels[x][0]], linewidth=5)
        plt.text(x+0.1, -0.5, labels[x][0], rotation=75, color=ftype_colors[labels[x][0]], size=LARGE_FONT_SIZE)

    coef = (len(str(max(df.max()))) + 5) * 0.2
    len_col = len(df.columns)
    len_row = len(df.index)

    if len_col < 8:
        len_col = 8
    if len_row < 8:
        len_row = 8

    fig = plt.gcf()
    fig.set_size_inches(coef * len_col, coef * len_row)
    ax.set_frame_on(False)
    ax.set_yticks(np.arange(df.shape[0]) + 0.5, minor=False)
    ax.set_xticks(np.arange(df.shape[1]) + 0.5, minor=False)
    ax.invert_yaxis()
    ax.xaxis.tick_top()
    ax.set_xticklabels(df.columns, minor=False, ha='left', fontsize=SMALL_FONT_SIZE)
    ax.set_yticklabels(df.index, minor=False, fontsize=SMALL_FONT_SIZE)
    plt.xticks(rotation=75)
    ax.grid(False)

    for y in range(df.shape[0]):
        for x in range(df.shape[1]):
            if df[df.columns[x]][df.index[y]] == 0:
                    continue
            else:
                val = int(round(df[df.columns[x]][df.index[y]] * 100.0 / number_frames))
                co = 'black'
                if val == 0:
                    val = '*'
                elif  val > 75:
                    co = 'white'

                plt.text(x + 0.5, y + 0.5, val, horizontalalignment='center', verticalalignment='center', color=co, size=SMALL_FONT_SIZE)

    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size=coef, pad=0.05)
    c = fig.colorbar(heatmap, cax=cax, ticks=np.linspace(0, number_frames, 11, dtype=int, endpoint=True))
    c.set_label('Number of frames')
    labels = np.linspace(0, 100, 11, dtype=int, endpoint=True)
    labels = labels.astype('str')
    labels[0] = '0<*<1'
    c.ax.set_yticklabels(labels, fontsize=SMALL_FONT_SIZE)
    c.set_label('% of frames', fontsize=SMALL_FONT_SIZE)

    if output is not None:
        plt.savefig(output, bbox_inches="tight")

def getPh4FingerprintDictionary(ph4_interaction_dictionary, global_ph4_interaction_list):
    tmp_list = [x[0] for x in global_ph4_interaction_list]
    global_ph4_fingerprint_dictionary = defaultdict(lambda:[0 for y in tmp_list])

    for frame_number, dict_int in ph4_interaction_dictionary.items():
        for list_int in dict_int.values():
            for tmp_int in list_int:
                global_ph4_fingerprint_dictionary[frame_number][tmp_list.index(tmp_int)] = 1
    return global_ph4_fingerprint_dictionary

def getPh4TimeSeries(ph4_fingerprint_dictionary, global_ph4_interaction_list):
    global_ph4_fingerprint_dictionary = {}
    nb_int = len(global_ph4_interaction_list)
    nb_frame = len(ph4_fingerprint_dictionary)

    for i in range(nb_int):
        global_ph4_fingerprint_dictionary[i] = [0]* nb_frame

    for i in range(nb_int):
        for j in range(nb_frame):
            global_ph4_fingerprint_dictionary[i][j] = ph4_fingerprint_dictionary[j][i]

    df = pd.DataFrame.from_dict(global_ph4_fingerprint_dictionary)
    df.columns=[(x[0].interaction_type, x[0].env_residue, str(x[0].env_atom), x[0].lig_residue, str(x[0].lig_atom)) for x in global_ph4_interaction_list]

    return [(global_ph4_interaction_list[x][0], global_ph4_interaction_list[x][1], global_ph4_fingerprint_dictionary[x]) for x in range(nb_int)], df

def getDataframeIM2(ph4_time_series):
    col = sorted(list({tuple([x[0].getInteractionType()] + x[0].getEnvironment() + x[0].getLigand()) for x in ph4_time_series[0]}),key=lambda k: (k[0], k[1]))

    matrix = np.zeros((len(col), len(col)))
    for x in range(len(col)):
        for y in range(len(col)):
            matrix[x][y] = scipy.stats.pearsonr(ph4_time_series[0][x][2], ph4_time_series[0][y][2])[0]
    return pd.DataFrame(matrix, index=col, columns=col)

def plotCorrelationMap(df, output=None):
    plt.clf()
    plt.cla()
    plt.close()

    fig, ax = plt.subplots()
    df = df.fillna(0)
    cmap_copy = copy.copy(plt.cm.Blues)
    heatmap = ax.pcolor(df, cmap=cmap_copy, alpha=0.9, vmax=1, vmin=-1)
    cmap_copy.set_over('black')

    coef = (len(str('{:.2f}'.format(max(df.max())))) + 2) * 0.2
    len_col = len(df.columns)*coef
    len_row = len(df.index)*coef
    if len_col < 5:
        len_col = 5
    if len_row < 5:
        len_row = 5

    fig = plt.gcf()
    fig.set_size_inches(len_col, len_row)
    ax.set_frame_on(False)

    ax.set_yticks(np.arange(df.shape[0]) + 0.5, minor=False)
    ax.set_xticks(np.arange(df.shape[1]) + 0.5, minor=False)
    ax.invert_yaxis()
    ax.xaxis.tick_top()

    ax.set_xticklabels(df.columns, minor=False, ha='left')
    ax.set_yticklabels(df.index, minor=False)
    plt.xticks(rotation=75)
    ax.grid(False)

    divider = make_axes_locatable(ax)

    correlation_class = {"Strong negative correlation [-1,-0.6]":set(), "Negative correlation ]-0.6,-0.4]":set(), "Positive correlation [0.4,0.6[":set(), "Strong positive correlation [0.6,1]":set()}
    for y in range(df.shape[0]):
        for x in range(df.shape[1]):
            co = 'black'
            val = round(df[df.columns[x]][df.index[y]], 2)
            if val >= 0.75:
                co = 'white'
            tmp = None
            plt.text(x + 0.5, y + 0.5, val, horizontalalignment='center', verticalalignment='center', color=co)
            if val<=-0.6:
                tmp = "Strong negative correlation [-1,-0.6]"
            elif -0.6<val<=-0.4:
                tmp = "Negative correlation ]-0.6,-0.4]"
            elif 0.4<=val<0.6:
                tmp = "Positive correlation [0.4,0.6["
            elif 0.6<=val:
                tmp = "Strong positive correlation [0.6,1]"

            if y>x and tmp is not None and df.columns[x]!= df.columns[y]:
                correlation_class[tmp].add((df.columns[x], df.columns[y], val))

    cax = divider.append_axes('right', size=coef, pad=0.05)
    c = fig.colorbar(heatmap, cax=cax, ticks=np.linspace(-1, 2, 11, dtype=int, endpoint=True))
    c.set_label('Pearson correlation coefficient')

    with open(output[:-4]+'.txt', 'w') as text_file:
        for key, value in sorted(correlation_class.items()):
            text_file.write('{}\n{:^50}\n{}\n\n'.format('*'*50, key, '*'*50))
            for x in sorted(value, key=lambda k: k[2]):
                if len(x[0]) == 5:
                    text_file.write(
                    "Column:\n{:35}: {}\n{:35}: {:12} {}\n{:35}: {:12} {}\n\nRow:\n{:35}: {}\n{:35}: {:12} {}\n{:35}: {:12} {}\n\n{:35}: {}\n{}\n\n".format(
                        ' ' * 4 + 'Type of interaction', x[0][0], ' ' * 4 + 'Ligand information', x[0][3], x[0][4],
                        ' ' * 4 + 'Environment information', x[0][1], x[0][2], ' ' * 4 + 'Type of interaction', x[1][0],
                        ' ' * 4 + 'Ligand information', x[1][3], x[1][4], ' ' * 4 + 'Environment information', x[1][1],
                        x[1][2], 'Pearson correlation coefficient', x[2], '*' * 5))
                elif len(x[0]) == 4:
                    text_file.write(
                        "Column:\n{:35}: {}\n{:35}: {:12} {}\n{:35}: {:12}\n\nRow:\n{:35}: {}\n{:35}: {:12} {}\n{:35}: {:12}\n\n{:35}: {}\n{}\n\n".format(
                            ' ' * 4 + 'Type of interaction', x[0][0], ' ' * 4 + 'Ligand information', x[0][2], x[0][3],
                            ' ' * 4 + 'Environment information', x[0][1], ' ' * 4 + 'Type of interaction', x[1][0],
                            ' ' * 4 + 'Ligand information', x[1][2], x[1][3], ' ' * 4 + 'Environment information',
                            x[1][1], 'Pearson correlation coefficient', x[2], '*' * 5))

            text_file.write('\n')

    if output is not None:
        plt.savefig(output, bbox_inches="tight")

def renameAa(ph4_interaction_dictionary, water_name = "HOH"):
    tmp_dictionary = dict(ph4_interaction_dictionary)

    for frame_number, dict_int in ph4_interaction_dictionary.items():
        for interaction_type, list_int in dict_int.items():
            tmp_dictionary[frame_number][interaction_type] = []
            tmp_check_list = []
            wat_check_list = []
            for tmp_int in list_int:
                if tmp_int.env_residue.split('_')[0]==water_name:
                        if (tmp_int.interaction_type, tmp_int.lig_atom) not in wat_check_list:
                            tmp_dictionary[frame_number][interaction_type].append(tmp_int)
                            tmp_dictionary[frame_number][interaction_type][-1].env_atom = None
                            tmp_dictionary[frame_number][interaction_type][-1].env_residue = "WAT_1_" + tmp_dictionary[frame_number][interaction_type][-1].env_residue.split('_')[2]
                            wat_check_list.append((tmp_int.interaction_type, tmp_int.lig_atom))

                else:
                    if (tmp_int.interaction_type, tmp_int.env_residue, tmp_int.lig_atom) not in tmp_check_list:
                        tmp_dictionary[frame_number][interaction_type].append(tmp_int)
                        tmp_dictionary[frame_number][interaction_type][-1].env_atom = None
                        tmp_check_list.append((tmp_int.interaction_type, tmp_int.env_residue, tmp_int.lig_atom))

    return tmp_dictionary
