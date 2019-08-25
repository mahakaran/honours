###########################################
#   Written by Mahakaran Sandhu, 2019.    #
#          Jackson Lab, ANU               #
###########################################


import numpy as np
import os
import argparse




parser = argparse.ArgumentParser(description='Perform hydrogen mass repartioning on a GROMACS include-topology file (ITP).')

parser.add_argument('infile', help='Input ITP file on which to perform HMR')
parser.add_argument('outfile', help='Name of output file')
parser.add_argument('hmr', help='Mass value to increase hydrogen masses by')

args = parser.parse_args()






#####################                           FUNCTION DEFINITIONS                           ##########################

def read_itp(itpfile): 
    
    """ Reads a GROMACS include topology (ITP) file to return an list of lists that contains the 
    atom data per residue."""
    with open(itpfile) as top: 
        top_data = top.readlines()
    begin_ndx    = [top_data.index(i) for i in top_data if '[ atoms ]' in i]
    end_ndx      = [top_data.index(i) for i in top_data if '[ bonds ]' in i]
    start_read   = begin_ndx[0]+2
    end_read     = end_ndx[0]
    atom_data    = top_data[start_read:end_read]
    res_ndx      = [atom_data.index(i) for i in atom_data if 'residue' in i]   
   
    residue_list = []
    for i in range(len(res_ndx)-1):
        if i != (len(res_ndx)):
            residue_data = atom_data[res_ndx[i]:res_ndx[i+1]]
            residue_list.append(residue_data)
            

    residue_list.append(atom_data[res_ndx[-1]:end_ndx[0]])
        
        
    sorted_residue_list = []
    for residue in residue_list:
        resdata_sorted = []
        for atom in residue:
            atomdata_sorted = atom.split()
            resdata_sorted.append(atomdata_sorted)
        sorted_residue_list.append(resdata_sorted)

    return (residue_list, sorted_residue_list, top_data, start_read, end_read)
        



def doHMR(sorted_res_list, hmr_mass):
    """Breaks down the input into individual atoms and their properties; then performs the hydrogen mass 
    repartitioning on relevant atoms. Returns a list of lists. HMR_mass is the mass with which you want to 
    increase the hydrogen mass by."""
    
    
    #definitions (how many hydrogens attached to this atomtype)
    n_term = ['NMET', 'NARG', 'NLYS', 'NASP', 'NGLU', 'NHIS', 'NGLY', 'NALA', 'NVAL', 'NLEU',
              'NILE', 'NPHE', 'NTYR', 'NSER', 'NTHR', 'NGLN', 'NASN', 'NCYS', 'NPRO']
    defs = {
     'MET':{'N':1, 'CA':1, 'CB':2,'CG':2,'CE':3  },
     'ARG':{'N':1, 'CA':1, 'CB':2,'CG':2,'CD':2,'NE':1, 'NH1':2, 'NH2':2}, 
     'LYS':{'N':1, 'CA':1, 'CB':2,'CG':2,'CD':2,'CE':2,'NZ':3},
     'ASP':{'N':1, 'CA':1, 'CB':2},
     'GLU':{'N':1, 'CA':1, 'CB':2,'CG':2},
     'HIS':{'N':1, 'CA':1, 'CB':2,'CE1':1,'NE2':1,'CD2':1},
     'GLY':{'N':1, 'CA':2},
     'ALA':{'N':1, 'CA':1, 'CB':3},
     'VAL':{'N':1, 'CA':1, 'CB':1,'CG1':3,'CG2':3},
     'LEU':{'N':1, 'CA':1, 'CB':2,'CG':1,'CD1':3,'CD2':3},
     'ILE':{'N':1, 'CA':1, 'CB':1,'CG2':3,'CG1':2,'CD':3},
     'PHE':{'N':1, 'CA':1, 'CB':2,'CD1':1,'CE1':1,'CZ':1,'CE2':1,'CD2':1},
     'TYR':{'N':1, 'CA':1, 'CB':2,'CD1':1,'CE1':1,'OH':1,'CE2':1,'CD2':1},
     'SER':{'N':1, 'CA':1, 'CB':2,'OG':1},
     'THR':{'N':1, 'CA':1, 'CB':1,'CG2':3,'OG1':1},
     'GLN':{'N':1, 'CA':1, 'CB':2,'CG':2,'NE2':2},
     'ASN':{'N':1, 'CA':1, 'CB':2,'ND2':2},
     'CYS':{'N':1, 'CA':1, 'CB':2,'SG':1},
     'PRO':{'CD':2,'CG':2, 'CB':2,'CA':1},
     'ACK':{'N':1, 'CA':1, 'CB':2,'CG':2,'CD':2,'CE':2,'NZ':1,'C2':3}
    }    
    
    init_masses = []
    fin_masses = []
    hmr_all = []
    
    
    
    for residue in sorted_res_list:  
        resname = residue[0][5]

        hmr_res = []
        hmr_res.append(residue[0])          
        for atom_dat in residue[1:]:
            if len(atom_dat) !=0:
                #print (atom_dat)
                atom_num    = atom_dat[0]
                atom_name   = atom_dat[1]
                res_num     = atom_dat[2]
                res_name    = atom_dat[3]            
                atomtype    = atom_dat[4]
                atom_num2   = atom_dat[5]
                atom_charge = atom_dat[6]
                atom_mass   = float(atom_dat[7])
                init_masses.append(atom_mass)

              
                if atomtype[0]=='H':
                    new_mass = atom_mass+hmr_mass
                    
                elif defs.get(res_name) != None:
                    hmr_factor = defs.get(res_name).get(atomtype)
                    if (resname=='NPRO') and atomtype[0]=='N':
                        new_mass = atom_mass-(2*hmr_mass)
                    elif (resname in n_term) and atomtype[0]=='N':
                        new_mass = atom_mass-(3*hmr_mass)
                    elif hmr_factor != None:
                        new_mass = atom_mass-(hmr_factor*hmr_mass)                                         
                    else:
                        new_mass = atom_mass
                        
                elif defs.get(res_name) == None:
                    new_mass = atom_mass           
                    
                fin_masses.append(new_mass)
                hmr_atom = [atom_num,atom_name,res_num,res_name,atomtype, atom_num2,atom_charge,new_mass]

            hmr_res.append(hmr_atom)
        hmr_all.append(hmr_res)

    if round(np.sum(init_masses)) == round(np.sum(fin_masses)):
        print ('REPARTITIONING SUCCESSFUL:\n ')
        print ('INITIAL MASS: '+str(np.sum(init_masses)))
        print ('FINAL MASS: ' + str(np.sum(fin_masses)))
        return hmr_all
    else:
        print ('MASS DIFFERENCE:')
        print ('INITIAL MASS:'+str(np.sum(init_masses)))
        print ('FINAL MASS:' + str(np.sum(fin_masses)))
        return (hmr_all, init_masses, fin_masses)
        
        #raise Exception('REPARTITIONING FAILED: initial mass of molecule does not equal final mass. Initial mass is: {} Da.'.format(np.sum(init_masses))+ ' Final mass is: {} Da'.format(np.sum(fin_masses)))



def writeITP(readITP_output, HMRed_dat, outname):
    """Does the reverse of readITP()"""
    with open(outname, 'w+') as outfile:
      outfile.write(';THIS FILE HAS BEEN GENERATED BY GROMACS_HMR. USE WITH CAUTION.\n')
      for i in readITP_output[2][:topology[-2]]:
          outfile.write(i)
      for i in HMRed_dat:
          for j in i:
              for k in j:
                  outfile.write( (str(k)+'   '))
              outfile.write('\n')
      for i in readITP_output[2][topology[-1]:]:
          outfile.write(i)
      outfile.close()
    
################################################################################################################################################################################################################


##################################                                      OUTPUT WORKFLOW                                                                                 ########################################


topology = read_itp(args.infile)
HMR_topology = doHMR(topology[1], int(args.hmr))
writeITP(topology,HMR_topology, args.outfile)

print('gromacsHMR has successfully transformed the input topology. Great Success!\n Remember to remove the final atom line in the output topology before use!')
