import numpy as np
angular_symbols = np.array(['s','p','d','f','g','h'])
angular_kappa   = np.array([ -1, -2, -3, -4, -5, -6])


def write_many_orbitals(core_orbitals,peel_orbitals,occupations,f):
    #first do core
    #f = open('DSTG2.INP','w')
    for ii in range(0,len(core_orbitals)):
        princ_n = int(core_orbitals[ii][0:len(core_orbitals[ii])-1])
        orbital_symbol = core_orbitals[ii][-1]
        index_of_angular_symbol = np.where(angular_symbols == orbital_symbol)[0]
        kappa = angular_kappa[index_of_angular_symbol][0]
        orbital_string = make_orbital_string(princ_n,kappa,[])
        f.write(orbital_string)
    for ii in range(0,len(peel_orbitals)):
        princ_n = int(peel_orbitals[ii][0:len(peel_orbitals[ii])-1])
        orbital_symbol = peel_orbitals[ii][-1]
        index_of_angular_symbol = np.where(angular_symbols == orbital_symbol)[0]
        kappa = angular_kappa[index_of_angular_symbol][0]
        orbital_string = make_orbital_string(princ_n,kappa,occupations[:,ii])
        f.write(orbital_string)

    return 0

def make_orbital_string(princ_n,kappa,occupations):

    orbital = '&ORB PRINC={} KAPPA={} '

    string = orbital.format(princ_n,kappa)

    if len(occupations) > 1:
        string += 'CSF= '
        counter = 0
        for ii in range(0,len(occupations)):
            counter += 1
            string += str(int(occupations[ii])) +' '
            if (counter > 5) and (ii < len(occupations)-1):
                 string +=', \n CSF={}* '.format(ii+1)
                 counter = 0

    string += '/ \n'
   # print(string)

    return string



def write_partial_wave():
    return 0

def write_dstg2(num_levels,num_partial_waves,core_orbitals,peel_orbitals,occupations):
    total_orbtials = len(core_orbitals) + len(peel_orbitals)

    num_csfs = np.shape(occupations)[0]

    f = open('DSTG2.INP','w')

    f.write('DSTG2: DARC\n')
    f.write('&PREINP NPW={} NP_PER_SUBWORLD=1 IDIMCHECK=0 IANGULAR=-1 /\n'.format(num_partial_waves))
    f.write('&DSTG2 NWM={} NMAN={} INAST={} NAST={}/\n'.format(total_orbtials,num_csfs,num_partial_waves,num_levels))

    write_many_orbitals(core_orbitals,peel_orbitals,occupations,f)

    f.write(' &ANGOPT /\n')
    f.write(' &JVALUE /\n')

    num_electrons = int(np.sum(occupations[0,:]))

    if (num_electrons%2) == 0:
        jstart = 0.5
    else:
        jstart = 0.0
    j = jstart
    for ii in range(0,int(num_partial_waves/2)):
        f.write( '&SYM JTOT={}  NPTY= 1 /\n'.format(j))
        f.write( '&SYM JTOT={}  NPTY=-1 /\n'.format(j))
        j += 1

#write_orbital(1,-1,[1,1,1,1,1,1,1,1,1])