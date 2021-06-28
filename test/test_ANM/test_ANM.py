import unittest
import os

os. chdir('/home/marta/Documentos/TFM/ProgramaTFM')
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
import numpy as np
import json
from ANM.ANM_Simulation import  *



class TestANM(unittest.TestCase):
    def test_anm(self):
        with open('/home/marta/Documentos/TFM/ProgramaTFM/test/test_ANM/input/test_anm.json') as file:
            Information = json.load(file)
        if os.path.isfile('/home/marta/Documentos/TFM/ProgramaTFM/test/test_ANM/output/Trajectory_ANM.pdb'):
            os.remove('/home/marta/Documentos/TFM/ProgramaTFM/test/test_ANM/output/Trajectory_ANM.pdb') 
        S = Information["Simulation"]
        Sim = "P"
        ANMSimulation(S)
        
        
        f1 = "/home/marta/Documentos/TFM/ProgramaTFM/test/test_ANM/input/Trajectory_ANM.pdb"
        f2 = "/home/marta/Documentos/TFM/ProgramaTFM/test/test_ANM/output/Trajectory_ANM.pdb"
          
        positions1 = np.array(PDBFile(f1).positions.value_in_unit(nanometer))
        positions2 = np.array(PDBFile(f2).positions.value_in_unit(nanometer))
        
        np.testing.assert_array_equal(positions1, positions2)
        
        
if __name__ == '__main__':
    unittest.main()
    

    
 