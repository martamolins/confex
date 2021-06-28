from prody import *
from pylab import *
from matplotlib.pylab import *
import os
        
ion()



class ANMSim(object):
    def __init__(self, S):
        self.address = S["Input_address"]
        self.address_output_pdb_ANM = S["Output_address"]+'/Trajectory_ANM.pdb'
        self.address_output= S["Output_address"]
        self.N_conf = S["N_conf"]
        self.structure_name = S["Structure_name"]
        self.rmsd = S["RMSD"]
        
    def ANMPreparation(self):
        if not os.path.exists(self.address_output):
            os.makedirs(self.address_output)
        pdb = parsePDB(self.address) 
        return pdb
    
    def ANMCalculation(self, pdb):
        
        print('PDB loaded')        
        anm = ANM('%s' % self.structure_name)    
        anm.buildHessian(pdb.select('calpha'))
        
        anm.calcModes()
        print('Modes calculated')
        
        anm, atoms = extendModel(anm, pdb.calpha, pdb.protein, norm=True)
        return anm   
        
    def ANMConformations(self, pdb, anm, conformations = None, address_out = None):
        if conformations is None:
            conformations = self.N_conf
        if address_out is None:
            address_out = self.address_output_pdb_ANM
            
        ensemble = sampleModes(anm[:3], pdb.protein, n_confs = conformations, rmsd=self.rmsd)

        pdb2 = (pdb.protein).copy()    
        pdb2.addCoordset(ensemble)
        
        writePDB(address_out, pdb2)
        print('Done')
        return anm
        
def ANMSimulation(S):
    ANM_S = ANMSim(S)
    pdb = ANM_S.ANMPreparation()
    anm = ANM_S.ANMCalculation(pdb)
    ANM_S.ANMConformations(pdb, anm)
