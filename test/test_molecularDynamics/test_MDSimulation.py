import unittest
import os


from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
import numpy as np
import json
from MolecularDynamics.MD_Simulation import *
import csv
from prody import *





class TestMD(unittest.TestCase):
    def test_md(self):
        with open(os.getcwd()+'/test/test_molecularDynamics/input/test_md.json') as file:
            Information = json.load(file)
            
        address_states_out = Information["Output_address"]+'/States_MD.csv'
        address_trajectory_out =  Information["Output_address"]+'/Trajectory_MD.pdb'
        
        if os.path.isfile(address_states_out):
            os.remove(address_states_out)
        if os.path.isfile(address_trajectory_out):
            os.remove(address_trajectory_out)

         # MD
        
        MolecularDynamicsSimulation(Information, 'MD', 1)
                
        
        # Input Sequence
        pdb_in = parsePDB(os.getcwd()+'/test/test_molecularDynamics/input/Trajectory_MD.pdb')
        a1_in = pdb_in.getNames()
        atoms_in = []
        for a_in in a1_in:
            atoms_in.append(a_in)        
        
        
        address_states_out = Information["Output_address"]+'/States_MD.csv'
        address_trajectory_out =  Information["Output_address"]+'/Trajectory_MD.pdb'
        
        # Number of states
        states = []
        with open(address_states_out, newline='') as csvfile:
            file = csv.reader(csvfile, delimiter=' ', quotechar='|')
            for row in file:
                states.append(row)
                
        # Atoms obtained
        pdb = parsePDB(address_trajectory_out) 
        a1 = pdb.getNames()
        atoms = []
        for a in a1:
            atoms.append(a)
            


        assert len(states) == Information["Conformations"]+1
        np.testing.assert_array_equal(atoms_in, atoms)
        
        
if __name__ == '__main__':
    unittest.main()
    
