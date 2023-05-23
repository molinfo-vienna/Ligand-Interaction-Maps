"""
@author: Arthur Garon, Thomas Seidel
Department of Pharmaceutical Sciences
University of Vienna
"""

import sys
import time
import argparse
import pickle

from common import *


def parseArguments():
    parser = argparse.ArgumentParser(description='>>> Generates the full ligand-receptor interactions map output for an input MD trajectory.')
    parser.add_argument('-trj',
                        dest='trajectory',
                        required=True,
                        help='[Required] The path of the topology file',
                        nargs="*")
    parser.add_argument('-top',
                        dest='topology',
                        required=True,
                        help='[Required] The path of the trajectory file, can also handle multiple trajectories of the same system',
                        nargs=1)
    parser.add_argument('-lig',
                        dest='ligand_code',
                        required=True,
                        help='[Required] The 3-letters code of the ligand',
                        nargs=1)
    parser.add_argument('-cs',
                        dest='chunk_size',
                        help='[Optional] The number of frames to consider per chunk (default: all frames)',
                        nargs=1,
                        default=None)
    parser.add_argument('-o',
                        dest='output',
                        help='[Optional] The output folder where the cdf file will be generated (Default: current directory)',
                        nargs=1,
                        default=None)
    parser.add_argument('-alh',
                        dest='alh',
                        required=False,
                        action='store_true',
                        default=False,
                        help='Generate atom-level hydrophobic features')
    
    return parser.parse_args()

if __name__ == '__main__':
    args = parseArguments()

    initial_time = time.time()
    trajectory = args.trajectory
    topology = args.topology[0]
    name = os.path.basename(trajectory[0])[:-4]
    ligand_code = args.ligand_code[0]
   
    if args.chunk_size is None:
        chunk_size = 500
    else:
        chunk_size = int(args.chunk_size[0])

    if args.output is None:
        output = './'
    else:
        output = args.output[0] + '/'

    chunk_number = cdfMol(topology, trajectory, output, name, chunk_size=chunk_size)
    mergeCDFMolecule(output + name + '_chunk_0.cdf', chunk_number)
    
    print('> Global trajectory object generation ...')
    ph4_interaction_dictionary = getPh4InteractionDictionary(output + name + '.cdf', ligand_code, args.alh)
 
    with open(output + name + '.gt', 'wb') as handle:
        pickle.dump(ph4_interaction_dictionary, handle, pickle.HIGHEST_PROTOCOL)

    with open(output + name + '.gt', 'rb') as handle:
        ph4_interaction_dictionary = pickle.load(handle)

    frame_list = [min(ph4_interaction_dictionary.keys()), max(ph4_interaction_dictionary.keys())]
    
    print('> Generation of images ...')
    
    drawLigand(output + name + '.cdf', ligand_code, output)

    ph4_interaction_dictionary = {key:ph4_interaction_dictionary[key] for key in xrange(frame_list[0],frame_list[1])}

    """
    global_ph4_interaction_list = getGlobalPh4InteractionList(ph4_interaction_dictionary)
    df = getDataframeIM(global_ph4_interaction_list)
    plotInteractionMap(df, number_frames=frame_list[1] - frame_list[0],
                       output=output + name + '_full_interaction_map.pdf')

    ph4_fingerprint_dict = getPh4FingerprintDictionary(ph4_interaction_dictionary, global_ph4_interaction_list)
    ph4_time_series = getPh4TimeSeries(ph4_fingerprint_dict, global_ph4_interaction_list)
    df = getDataframeIM2(ph4_time_series)
    plotCorrelationMap(df, output=output + name + '_full_correlation_map.pdf')
    """
    ph4_interaction_dictionary = renameAa(ph4_interaction_dictionary)
    global_ph4_interaction_list = getGlobalPh4InteractionList(ph4_interaction_dictionary)
    df = getDataframeIM(global_ph4_interaction_list)
    
    plotInteractionMap(df, number_frames=frame_list[1] - frame_list[0],
                       output=output + name + '_interaction_map.pdf')

    ph4_fingerprint_dict = getPh4FingerprintDictionary(ph4_interaction_dictionary, global_ph4_interaction_list)
    ph4_time_series = getPh4TimeSeries(ph4_fingerprint_dict, global_ph4_interaction_list)
    df = getDataframeIM2(ph4_time_series)
    
    plotCorrelationMap(df, output=output + name + '_correlation_map.pdf')

    calc_time = time.time() - initial_time
    
    print('> Interaction and Correlation maps generated in {}s'.format(int(calc_time)))
