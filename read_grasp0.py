import numpy as np 

def read_grasp0_inp(file_path):
    
    graspinp = open(file_path,'r')
    orbitals = []

    graspinp_read = graspinp.readlines()
    
    csf_orbs = graspinp_read[1].split()
    num_csfs = int(csf_orbs[0])
    num_orbs = int(csf_orbs[1]) 

    occupations = np.zeros([num_orbs,num_csfs])

    for ii in range(2,2+num_orbs):
        current_line = graspinp_read[ii].split()

        #print(current_line)
        orbitals.append(current_line[0].lower())
    
        if len(current_line) == 2:  
            occupations[ii-2,:] = int(current_line[1]) 
        else:
            for i in range(1,len(current_line)):
                occupations[ii-2,i-1] = int(current_line[i])


        #print(occupations[ii-2] )

    graspinp.close()

    return orbitals,occupations

#read_grasp0_inp('GRASP.INP')