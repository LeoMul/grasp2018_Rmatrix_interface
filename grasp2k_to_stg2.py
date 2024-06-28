import read_grasp2k
import stg2_lib
import numpy as np 
import argparse

parser = argparse.ArgumentParser()

parser.add_argument('-f', '--file',  help='Specify path of grasp2k csfs.out file')
parser.add_argument('-l','--num_levels',  help='Number of close-coupling levels',type= int)
parser.add_argument('-p', '--num_partial_waves',  help='Number of N+1 partial waves',type=int)
args = parser.parse_args()

def main(num_partial_waves,num_levels_cc,csf_file):

    f = open(csf_file,'r')

    csf_file_read = f.readlines()
    peel_subshells = read_grasp2k.find_peel_subshells(csf_file_read)

    core_subshells = read_grasp2k.find_core_subshells(csf_file_read)

    csf_index = read_grasp2k.find_csfs(csf_file_read)

    rel_csfs = read_grasp2k.collect_all_rel_valance_csf(csf_file_read,peel_subshells,csf_index)


    nr_peel,position_in_r_peel,needs_adding = read_grasp2k.find_nr_peel_subshells(peel_subshells)

    nr_csfs = np.zeros([np.shape(rel_csfs)[0],len(nr_peel)])


    for ii in range(0,len(nr_peel)):
        nr_csfs[:,ii] = rel_csfs[:,position_in_r_peel[ii]]
        if needs_adding[ii]:
            nr_csfs[:,ii] += rel_csfs[:,position_in_r_peel[ii]-1] 

    nr_csfs = np.unique(nr_csfs,axis=0)
    nr_csfs = np.flip(nr_csfs,axis=0)
    #num_orbitals = len(nr_peel)

    #print_out_orbitals(num_orbitals,nr_csfs,nr_peel,core_subshells)

    stg2_lib.write_dstg2(num_levels_cc,num_partial_waves,core_subshells,nr_peel,nr_csfs)

    return 0

if not args.file:
    print("No csf file specified - assuming rcsfs.out")
    file = 'rcsfs.out'
else:
    file = args.file 

if not args.num_partial_waves:
    print("no partial waves requested")
    parser.print_help()
    exit()

if not args.num_levels:
    print("no cc levels requested")
    parser.print_help()
    exit()

main(num_partial_waves=args.num_partial_waves,
     num_levels_cc=args.num_levels,
     csf_file=file)
