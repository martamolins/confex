from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout

from prody import *
from pylab import *
from pdbfixer import *

import os
import numpy as np

ion()

from MolecularDynamics.MD_Simulation import *
from ANM.ANM_Simulation import *
import random
import time


class MetropolisCriteria(object):
    def __init__(self, S, Sim, seedset):
        self.stepsize = S["stepsize"]
        self.Temperature = S["Temperature"]
        self.Pressure = S["Pressure"]
        self.FF = S["Force_Field"]["Force_Field_name"]
        self.FF_WM = S["Force_Field"]["Water_Model"] 
        self.Name = S["Name"]
        self.SimType = Sim     
        self.rmsd = S["RMSD"]
        self.tol = S["Tolerance"]
        self.maxIt = S["maxIterations"]
        self.structure_name = S["Structure_name"]
        
        self.number_steps = S["N_conf"]        
        self.address_inicial = S["Input_address"]
        self.address_trajectory = S["Output_address"]+'/Trajectory_MC.pdb'
        self.address_pdb = S["Output_address"]+'/'+S["Structure_name"]+'_intermedi.pdb'
        self.address_scwrl = S["Output_address"]+'/'+S["Structure_name"]+'_intermedi2.pdb'
        self.address_states = S["Output_address"]+'/MC_states.txt'
        
        self.seed = seedset
        self.scwrl_bin = S["SCWRL_BIN"]
        self.structure_name = S["Structure_name"]
        
           
    def FirstConformation(self):
        pdb, system, modeller = MolecularDynamics.SystemPreparation(self)
        integrator, simulation = MolecularDynamics.SimulationPreparation(self, modeller, system, self.tol, self.maxIt)
        Energy_in = MolecularDynamics.EnergyCalculation(self, simulation)
        
        PDBFile.writeFile(modeller.topology, modeller.getPositions(), open(self.address_pdb,'w'))

        pdb_in = parsePDB(self.address_pdb) 
        os.remove(self.address_pdb) 
        
        return pdb_in, Energy_in


    def ANM_S(self, pdb, anm):                
        anm = ANMSim.ANMConformations(self, pdb, anm, 1, self.address_pdb)
        
        pdb2 = parsePDB(self.address_pdb) 
        os.remove(self.address_pdb) 
        
      
        pdb2.setACSIndex(1)
        pdb3 = (pdb.protein).copy()
        pdb3.delCoordset(0)
        pdb3.addCoordset(pdb2.getCoords())
        
        writePDB(self.address_pdb, pdb3)

        
        pdb_i = PDBFile(self.address_pdb)
        os.remove(self.address_pdb) 
        PDBFile.writeFile(pdb_i.topology, pdb_i.getPositions(), open(self.address_pdb,'w'))   

        
    def SideChains(self):
        cmd = '%s -i %s -o %s -v -h -t> /dev/null 2>&1' % (self.scwrl_bin, self.address_pdb, self.address_scwrl)
        os.system(cmd)
        
        self.fix(self.address_scwrl)
        
    
    def Intermediate(self, pdb_in, anm):
        self.ANM_S(pdb_in, anm)     
        self.SideChains()

        # Energy
        pdb_intermediate, system, pdb_int = MolecularDynamics.SystemPreparation(self, self.address_scwrl)

        try:
            integrator, simulation = MolecularDynamics.SimulationPreparation(self, pdb_int, system, self.tol, self.maxIt)
        except:
            return 1, 1
                

        Energy_fin = MolecularDynamics.EnergyCalculation(self, simulation)
        os.remove(self.address_scwrl) 
        PDBFile.writeFile(pdb_int.topology, pdb_int.getPositions(), open(self.address_pdb,'w'))
        
        return pdb_int, Energy_fin
    
    def fix(self, address):
        
        fixer = PDBFixer(filename=address)
        fixer.findMissingResidues()
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()

        
        os.remove(address)
        PDBFile.writeFile(fixer.topology, fixer.positions, open(address, 'w'))
        
    def Accepted(self, pdb_in, Energy_in, pdb_int, Energy_fin, steps_done, Pas):
        pdb_in = parsePDB(self.address_pdb) 
        Energy_in = Energy_fin

        file1 = open(self.address_trajectory, "a")
        file2 = open(self.address_pdb, "r")
        PDB_text = file2.read()
        file2.close()
        file1.write("Model %d \n" % (steps_done) + PDB_text )
        file1.close()
        
        file_states = open(self.address_states, "a")
        file_states.write(str(steps_done)+" "+  str(Pas) +" "+ str(Energy_in/kilojoule_per_mole)+'\n')
        file_states.close()
        
        steps_done = steps_done + 1
        print('Conformation accepted')
        return pdb_in, Energy_in, steps_done
            
    def MHAlgorithm(self, pdb_in, Energy_in, pdb_int, Energy_fin, steps_done, Pas):
        DifE = (Energy_fin - Energy_in)/kilojoule_per_mole
        print('E_in: '+ str(Energy_in) + '\nEfin: '+ str(Energy_fin)+ '\nDif: '+ str(DifE))

        if DifE <= 0:
            if DifE < -15000:
                pass
            else:
                pdb_in, Energy_in, steps_done = self.Accepted(pdb_in, Energy_in, pdb_int, Energy_fin, steps_done, Pas)
        if DifE > 0:
            P = np.exp(-DifE/self.Temperature[0])
            if self.seed != None:
                random.seed(self.seed)             
            rand_num = random.random()
            print('P: '+ str(P) + '   -   Random number: ' + str(rand_num))
            if P > rand_num:
                pdb_in, Energy_in, steps_done = self.Accepted(pdb_in, Energy_in, pdb_int, Energy_fin, steps_done, Pas)
            else:
                print('Conformation rejected')
        return pdb_in, Energy_in, steps_done

def Metropolis(S, Sim, seedset = None):
    MC = MetropolisCriteria(S, Sim, seedset)
    pdb_in, Energy_in = MC.FirstConformation()
    anm = ANMSim.ANMCalculation(MC, pdb_in)
    Pas = 0
    steps_done = 0
    error_atoms = 0

    while steps_done < MC.number_steps:
        MC.error = 0
        pdb_int, Energy_fin = MC.Intermediate(pdb_in, anm)

        a1 = pdb_in.numAtoms()
        pdbi = parsePDB(MC.address_pdb)
        a2 = pdbi.numAtoms() 
        if a1 != a2:
            print('Different atoms')
            error_atoms += 1
            continue
        print('Same atoms')

        # Metropolis-Hastings algorithm
        pdb_in, Energy_in, steps_done = MC.MHAlgorithm(pdb_in, Energy_in, pdb_int, Energy_fin, steps_done, Pas)
                
        os.remove(MC.address_pdb) 
        Pas = Pas + 1
        print('Number of conformations: '+ str(steps_done))
        print('Number of iterations: ' + str(Pas))

        if Energy_in/kilojoule_per_mole > 200000:
            break
        
