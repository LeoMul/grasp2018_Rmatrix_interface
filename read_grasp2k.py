import numpy as np
import stg2_lib 


angular_symbols  = np.array(['s','p','d','f','g','h'])
angular_kappa    = np.array([ -1, -2, -3, -4, -5, -6])
angular_core_occ = np.array([  2,  6, 10, 14, 18, 22])

def find_core_subshells(csf_file_read):
    for ii in range(0,len(csf_file_read)):
        current_line_split = csf_file_read[ii].replace('\n','').split()
        if len(current_line_split) > 0:
            if current_line_split[0] == 'Core':
                position =  ii+1 
                break 

    core_subshells = csf_file_read[position].replace('\n','').split()

    for core in core_subshells:
        if '-' in core:
            core_subshells.remove(core)

    return np.array(core_subshells)

def find_peel_subshells(csf_file_read):
    
    for ii in range(0,len(csf_file_read)):
        current_line_split = csf_file_read[ii].replace('\n','').split()
        if len(current_line_split) > 0:
            if current_line_split[0] == 'Peel':
                position =  ii+1 
                break 

    peel_subshells = csf_file_read[position].replace('\n','').split()

    return np.array(peel_subshells)

def find_csfs(csf_file_read):
    for ii in range(0,len(csf_file_read)):
        current_line_split = csf_file_read[ii].replace('\n','').split()
        if len(current_line_split) > 0:
            if current_line_split[0] == 'CSF(s):':
                return ii+1 
    return 0

def construct_rel_valence_csf(jj_csf,peel_subshells,relativistic_csfs):
    
    jj_csf_copy = jj_csf.replace('\n','')
    num_occupations = int(len(jj_csf_copy) / 9)

    start = 0
    end = 9
    occupation_numbers = np.zeros(len(peel_subshells))
    for ii in range(0,num_occupations):
        orbital_plus_occupation = jj_csf_copy[start:end]
        #print(orbital_plus_occupation)
        orbital = orbital_plus_occupation[0:5].replace(' ','')
        occupation = orbital_plus_occupation[6:8].replace(' ','')

        index = np.argwhere(peel_subshells == orbital)[0][0]
        occupation_numbers[index] = int(occupation)
        start = end 
        end += 9
    return occupation_numbers
    

def collect_all_rel_valance_csf(csf_file_read,peel_subshells,csf_index):


    relativistic_csfs = []
    offset = 0

    for ii in range(csf_index,len(csf_file_read),3):
        
        #print(offset)
        jj = ii + offset
        if jj > len(csf_file_read)-1:
            break
        #print(jj,csf_file_read[jj])
        if csf_file_read[jj].replace('\n','') ==' *':
            #print(jj)
            offset = offset + 1
            jj += 1 
            #print('hello,',jj,csf_file_read[jj])

        relativistic_csfs.append(construct_rel_valence_csf(csf_file_read[jj],peel_subshells,relativistic_csfs))


    return np.array(relativistic_csfs)

def find_nr_peel_subshells(peel_subshells):

    nr_peel = []
    position_in_r_peel = []

    needs_adding = []

    for ii in range(0,len(peel_subshells)):
        if '-' not in peel_subshells[ii]:
            nr_peel.append(peel_subshells[ii])
            position_in_r_peel.append(ii)

            if 's' in peel_subshells[ii]:
                needs_adding.append(False)
            else:
                needs_adding.append(True)


    #print(nr_peel)
    #print(position_in_r_peel)
    #print(needs_adding)

    return nr_peel,position_in_r_peel,needs_adding


def average_occupation_numbers(rel_csfs):
    return np.average(rel_csfs,axis=0)

def print_out_orbitals(num_orbitals,nr_csfs,nr_peel,core_subshells):

    for ii in range(0,len(core_subshells)):
        string = str(core_subshells[ii]).upper()+' '
        angular_symbol = str(core_subshells[ii])[-1]
        index_of_angular_symbol = np.where(angular_symbols == angular_symbol)
        occupation = angular_core_occ[index_of_angular_symbol][0]
        string += str(occupation)
        print(string)


    for ii in range(0,num_orbitals):
        orbital_occ = nr_csfs[:,ii]

        string = str(nr_peel[ii]).upper()+' '
        for jj in range(0,len(orbital_occ)):
            string =string + str(int(orbital_occ[jj]))+' '

        print(string)

def main():
    #for testing, ignore
    f = open('rcsf.out')

    csf_file_read = f.readlines()
    peel_subshells = find_peel_subshells(csf_file_read)

    core_subshells = find_core_subshells(csf_file_read)

    csf_index = find_csfs(csf_file_read)

    rel_csfs = collect_all_rel_valance_csf(csf_file_read,peel_subshells,csf_index)


    nr_peel,position_in_r_peel,needs_adding = find_nr_peel_subshells(peel_subshells)

    nr_csfs = np.zeros([np.shape(rel_csfs)[0],len(nr_peel)])


    for ii in range(0,len(nr_peel)):
        nr_csfs[:,ii] = rel_csfs[:,position_in_r_peel[ii]]
        if needs_adding[ii]:
            nr_csfs[:,ii] += rel_csfs[:,position_in_r_peel[ii]-1] 

    nr_csfs = np.unique(nr_csfs,axis=0)
    nr_csfs = np.reverse(nr_csfs,axis=1)

    num_orbitals = len(nr_peel)

    print_out_orbitals(num_orbitals,nr_csfs,nr_peel,core_subshells)

    stg2_lib.write_dstg2(150,24,core_subshells,nr_peel,nr_csfs)

#main()