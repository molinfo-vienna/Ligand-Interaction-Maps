"""
@author: Arthur Garon, Thomas Seidel
Department of Pharmaceutical Sciences
University of Vienna
"""

import time
import os
import pickle
import argparse

from common import *


def parseArguments():
    parser = argparse.ArgumentParser(description='>>> Generate a picture of the ligand with atom numbering, 2 interaction maps, a correlation map and a text file that summarize the positive and negative correlations of the correlation map.')
    parser.add_argument('-cdf',
                        dest='cdf',
                        required=True,
                        help='[Required] The path of the cdf file',
                        nargs=1)
    parser.add_argument('-gt',
                        dest='gt_object',
                        required=True,
                        help='[Required] The path of the gt object file',
                        nargs=1)
    parser.add_argument('-lig',
                        dest='ligand_code',
                        required=True,
                        help='[Required] The 3-letters code of the ligand',
                        nargs=1)
    parser.add_argument('-fl',
                        dest='frame_list',
                        help='[Optional] The list of frames to consider, e.g. [frame_start,frame_end] (Default: all frames)',
                        nargs=2,
                        default=None)
    parser.add_argument('-o',
                        dest='output',
                        help='[Optional] The output folder where the pictures will be generated (Default: current directory)',
                        nargs=1,
                        default=None)

    return parser.parse_args()

if __name__ == '__main__':
    args = parseArguments()

    initial_time = time.time()
    cdf = args.cdf[0]
    gt_object = args.gt_object[0]
    ligand_code = args.ligand_code[0]

    with open(gt_object, 'rb') as handle:
        ph4_interaction_dictionary = pickle.load(handle)
        
    if args.frame_list is None:
        frame_list = [min(ph4_interaction_dictionary.keys()), max(ph4_interaction_dictionary.keys())]
    else:
        frame_list = [int(args.frame_list[0]), int(args.frame_list[1])]
        
    if args.output is None:
        output = os.getcwd()
    else:
        output = args.output[0]

    drawLigand(cdf, ligand_code, output)

    ph4_interaction_dictionary = {key:ph4_interaction_dictionary[key] for key in xrange(frame_list[0],frame_list[1])}
    """
    global_ph4_interaction_list = getGlobalPh4InteractionList(ph4_interaction_dictionary)
    df = getDataframeIM(global_ph4_interaction_list)

    plotInteractionMap(df, number_frames=frame_list[1]-frame_list[0], output=output + os.path.basename(cdf)[:-4] + '_full_interaction_map.pdf')

    ph4_fingerprint_dict = getPh4FingerprintDictionary(ph4_interaction_dictionary, global_ph4_interaction_list)
    ph4_time_series = getPh4TimeSeries(ph4_fingerprint_dict, global_ph4_interaction_list)

    try:
        df = getDataframeIM2(ph4_time_series)
        plotCorrelationMap(df, output=output + os.path.basename(cdf)[:-4] + '_full_correlation_map.pdf')
    except:
        print('!!! Plotting correlation map failed')
        pass
    """    
    ph4_interaction_dictionary = renameAa(ph4_interaction_dictionary)
    global_ph4_interaction_list = getGlobalPh4InteractionList(ph4_interaction_dictionary)
    df = getDataframeIM(global_ph4_interaction_list)

    im_fname = os.path.join(output, os.path.basename(cdf)[:-4] + '_interaction_map.pdf')
    
    plotInteractionMap(df, number_frames=frame_list[1]-frame_list[0], output=im_fname)

    ph4_fingerprint_dict = getPh4FingerprintDictionary(ph4_interaction_dictionary, global_ph4_interaction_list)
    ph4_time_series = getPh4TimeSeries(ph4_fingerprint_dict, global_ph4_interaction_list)

    try:
        df = getDataframeIM2(ph4_time_series)
        cm_fname = os.path.join(output, os.path.basename(cdf)[:-4] + '_correlation_map.pdf')
        
        plotCorrelationMap(df, output=cm_fname)
        
    except:
        print('!!! Plotting correlation map failed')
        pass
 
    calc_time = time.time() - initial_time
    
    print('> Interaction and Correlation maps generated in {}s'.format(int(calc_time)))
