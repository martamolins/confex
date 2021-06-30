import unittest
import os

from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
import numpy as np
import json
from StructuralAnalysis.SimulationAnalysis import *

with open(os.getcwd()+'/test/test_structuralAnalysis/input/test_analysis.json') as file:
            Information = json.load(file)


class TestAnalysis(unittest.TestCase):
    def test_RMSD(self):
        rmsd1 = ([0.        , 0.85736664, 1.03940437, 1.27388586, 1.39749654,
        1.44948356, 1.45832844, 1.46058589, 1.48945237, 1.42424423,
        1.42553873, 1.41181371, 1.43091458, 1.52640545, 1.49551714,
        1.58217925, 1.57001111, 1.63571363, 1.5444326 , 1.50952326,
        1.54386284, 1.68161048, 1.60238063, 1.79657304, 1.7088092 ,
        1.61083145, 1.53842833, 1.51701644, 1.59208137, 1.56808092,
        1.59568148, 1.56267781, 1.61251131, 1.57057054, 1.62876329,
        1.65162585, 1.58404059, 1.65828061, 1.5189949 , 1.63232734,
        1.53525086, 1.61125903, 1.6314651 , 1.60818587, 1.60225947,
        1.65572109, 1.6007152 , 1.75948687, 1.7378134 , 1.69665129])
        

        A = SAnalysis(Information)
        rmsds = A.RMSD()
        
        rmsd2 = np.array(rmsds[0:50]).T
        rmsd2 = np.matrix.round(rmsd2,decimals = 8)
             
        np.testing.assert_array_equal(rmsd1, rmsd2)


    def test_RMSF(self):
        
        rmsf1 = ([0.54771155, 0.43360187, 0.36982576, 0.38782814, 0.38277796,
       0.4345859 , 0.58675408, 0.89493947, 1.22244483, 1.12547281,
       0.61516159, 0.53365482, 0.5076679 , 0.45849707, 0.42521871,
       0.47265563, 0.46538987, 0.61244563, 0.68667521, 0.6987402 ,
       0.58809972, 0.53412254, 0.44200366, 0.57316169, 0.47310176,
       0.38963909, 0.43791775, 0.46443025, 0.4384698 , 0.46667172,
       0.50906809, 0.60329003, 0.72707903, 0.62539378, 0.62879232,
       0.59633517, 0.59247447, 0.60858219, 0.73057446, 0.59206898,
       0.45873884, 0.45554263, 0.43154487, 0.49920721, 0.54229234,
       1.00049586, 1.01845807, 0.74509023, 0.53481202, 0.56786323])
        

        A = SAnalysis(Information)
        rmsfs = A.RMSF()
        
        rmsf2 = ((rmsfs[0:50]).T)
        rmsf2 = np.matrix.round(rmsf2,decimals = 8)
                
        np.testing.assert_array_equal(rmsf1, rmsf2)
                
        
    def test_PCA(self):
        
        eigen_in = -0.0021692607122784776
        
        A = SAnalysis(Information)
        pca = A.PCA()
        eigen = pca.getEigvecs()
        
        np.testing.assert_array_equal(eigen_in, eigen[1000][1])
                        
if __name__ == '__main__':
    unittest.main()
    

    
 
