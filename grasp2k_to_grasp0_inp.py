import read_grasp2k
import numpy as np

def main():
    csf_file ='rcsf.out'
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
    num_orbitals = len(nr_peel)
    print('GRASP: INPUT')
    print(np.shape(nr_csfs)[0],np.shape(nr_csfs)[1]+len(core_subshells),0)
    read_grasp2k.print_out_orbitals(num_orbitals,nr_csfs,nr_peel,core_subshells)
    
    print('ANG  11 10') 
    print('-1')                
    print('MCP')
    print('MCBP')
    print('MCT  1  0')         
    print('MCDF')
    print('0 0 0 ')            
    print('52')
    print('FIX 0')
    print('EAL 30') 
    print('BENA 20')
    print('')
    print('OSCL 20 9 10 11')
    print('')
    print('STOP')


    return 0

main()
