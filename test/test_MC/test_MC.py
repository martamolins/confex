import unittest
import os


from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
import numpy as np
import json
from MC.MC import *


class TestMD(unittest.TestCase):
    def test_md(self):
        with open(os.getcwd()+'/test/test_MC/input/test_MC.json') as file:
            Information = json.load(file)
        
        if os.path.isfile(os.getcwd()+'/test/test_MC/output/Trajectory_MC.pdb'):
            os.remove(os.getcwd()+'/test/test_MC/output/Trajectory_MC.pdb') 
        if os.path.isfile(os.getcwd()+'/test/test_MC/output/MC_states.txt'):
            os.remove(os.getcwd()+'/test/test_MC/output/MC_states.txt')  
        if os.path.isfile(os.getcwd()+'/test/test_MC/output/'+Information["Structure_name"]+'_intermedi.pdb'):
            os.remove(os.getcwd()+'/test/test_MC/output/'+Information["Structure_name"]+'_intermedi.pdb')
        if os.path.isfile(os.getcwd()+'/test/test_MC/output/'+Information["Structure_name"]+'_intermedi2.pdb'):
            os.remove(os.getcwd()+'/test/test_MC/output/'+Information["Structure_name"]+'_intermedi2.pdb')

        pdb_input = os.getcwd()+'/test/test_MC/input/Trajectory_MC.pdb'
        pdb_in = parsePDB(pdb_input)
        a1_in = pdb_in.getNames()
        atoms_in = []
        for a_in in a1_in:
            atoms_in.append(a_in)   


        Metropolis(Information, 'P', 43)
        
        states_out = Information["Output_address"]+'/MC_states.txt'
        pdb_out = Information["Output_address"]+'/Trajectory_MC.pdb'
        
        # Number of states
        with open(states_out) as f:
            states = f.readlines()
            
        # Atoms pdb out
        pdb = parsePDB(pdb_out) 
        a1 = pdb.getNames()
        atoms = []
        for a in a1:
            atoms.append(a)
            
         
        assert len(states) == Information["N_conf"]
        np.testing.assert_array_equal(atoms_in, atoms)
        
if __name__ == '__main__':
    unittest.main()
    
