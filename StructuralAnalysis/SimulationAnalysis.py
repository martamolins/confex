from prody import *
from pylab import *
ion()
import numpy as np
import matplotlib.pyplot as plt
import os

class SAnalysis(object):
    def __init__(self, Info):
        self.address_inicial = Info["Input_address"]
        # self.address_final = Info["Output_address"]
        self.struct_name = Info["Structure_name"]
        self.address_out = Info["Output_address"]
        self.address_dcd = Info["Output_address"]+Info["Structure_name"]+'.dcd'
        self.Name = Info["Name"]

    def preparation(self):
        if not os.path.exists(self.address_out):
            os.makedirs(self.address_out)

        pdb_out = parsePDB(self.address_inicial)
        # Alineament prot
        pdb_out.setACSIndex(0)
        alignCoordsets(pdb_out.protein)
        plt.ioff()
        return pdb_out

    def RMSD(self): 
        pdb_out = self.preparation()
        rmsds = calcRMSD(pdb_out.protein)
        # Plot
        # plt.figure()
        plt.plot(range(pdb_out.numCoordsets()), rmsds)
        plt.ylabel('RMSD ($\AA$)')
        plt.xlabel('Conformations')
        plt.ylim((0,5.5))
        plt.title('RMSD '+self.struct_name)
        plt.savefig(self.address_out+'RMSD_'+self.struct_name+'.png', dpi=300, bbox_inches='tight')
        plt.clf()
        

        np.savetxt(self.address_out+'RMSD_'+self.struct_name+'.txt', rmsds)
        return rmsds
      
    def RMSF(self):
        pdb_out = self.preparation()
        rmsfs = calcRMSF(pdb_out.calpha)
        # plt.figure()
        plt.plot(range(pdb_out.numAtoms('calpha')), rmsfs)
        plt.ylabel('RMSF ($\AA$)')
        plt.xlabel('Residue number')
        plt.ylim((0,7.5))
        plt.title('RMSF '+self.struct_name)
        plt.savefig(self.address_out+'RMSF_'+self.struct_name+'.png', dpi=300, bbox_inches='tight')
        plt.clf()
        
        np.savetxt(self.address_out+'RMSF_'+self.struct_name+'.txt',rmsfs)

        return rmsfs

    def PCA(self):
        pdb_out = self.preparation() 
        writeDCD(self.address_dcd,pdb_out.protein)
        dcd_out = parseDCD(self.address_dcd) 
        dcd_out.superpose()
        pca = PCA(self.struct_name)
        pca.buildCovariance(dcd_out)
        pca.calcModes()
        # plt.figure()
        showProjection(dcd_out, pca[:2], color='red', marker='.')
        plt.title('PCA '+self.struct_name)   
        plt.ylim((-2,2))
        plt.xlim((-2,2))
        plt.savefig(self.address_out+'PCA_'+self.struct_name+'.png', dpi=300, bbox_inches='tight')
        plt.clf()
        
        saveModel(pca, filename = self.address_out+'PCA_'+self.struct_name+'.anm.npz')

        return pca
        
def simulationAnalysis(Info):
    A = SAnalysis(Info)
    for a in A.Name:
        if a == 'RMSD':
            rmsds = A.RMSD()
        if a == 'RMSF':
            rmsfs = A.RMSF()
        if a == 'PCA':
            pca = A.PCA()
    
    
    
