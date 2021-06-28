from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
import numpy as np
import random
import os

    

class MolecularDynamics(object):
    def __init__(self, S, Sim, seed = None):
        self.address_inicial = S["Input_address"]
        self.address_output_pdb_MD = S["Output_address"]+'/Trajectory_MD.pdb'
        self.address_output_states = S["Output_address"]+'/States_MD.csv'
        self.Time = S["Time"]
        self.stepsize = S["stepsize"]
        self.Temperature = S["Temperature"]
        self.Pressure = S["Pressure"]
        self.tol = S["Tolerance"]
        # self.Integrator = S["Integrator"]
        self.FF = S["Force_Field"]["Force_Field_name"]
        self.FF_WM = S["Force_Field"]["Water_Model"] 
        self.Name = S["Name"]
        self.SimType = Sim
        self.nsteps =  self.Time/self.stepsize
        self.seed = seed
        self.structure_name = S["Structure_name"]
        self.conformations = S["Conformations"]
        
    def SystemPreparation(self, address = None):
        if address is None:
            address = self.address_inicial
        pdb = PDBFile(address)     
        print('PDB loaded')    
        forcefield = ForceField(self.FF, self.FF_WM)  
        print('Force field created')   
        modeller = Modeller(pdb.topology, pdb.positions)
        if self.seed != None:
            random.seed(self.seed) 
        modeller.addHydrogens(forcefield) 
        print('Hydrogens addded')
        if self.SimType == "MD":
            modeller.addSolvent(forcefield, model = 'tip3p',padding=1.0*nanometers, neutralize = True)        
            system = forcefield.createSystem(modeller.topology, nonbondedMethod=PME, nonbondedCutoff=1*nanometer, constraints=HBonds, switchDistance=0.9*nanometer) 
            system.addForce(MonteCarloBarostat(self.Pressure*bar, self.Temperature[0]*kelvin))
            print('System created')
        if self.SimType == "P":
            system = forcefield.createSystem(modeller.topology)
            print('System created')
        return pdb, system, modeller
    
    def SimulationPreparation(self, modeller, system, tolerance = None, maxIt = None):   
        integrator = LangevinMiddleIntegrator(self.Temperature[0]*kelvin, 1/picosecond, self.stepsize*picoseconds)
        print('Integrator created')
        if self.seed != None:
            integrator.setRandomNumberSeed(self.seed)
            print('seeed')
        simulation = Simulation(modeller.topology, system, integrator)
        print('Simulation created')
        simulation.context.setPositions(modeller.positions)
        if self.seed != None:
            simulation.context.setVelocitiesToTemperature(self.Temperature[0]*kelvin, self.seed)
            print('seed')
        else:
            simulation.context.setVelocitiesToTemperature(self.Temperature[0]*kelvin)
        print('Minimizing energy...')
        if tolerance == '' and maxIt == '':
            simulation.minimizeEnergy()
        if tolerance != '' and maxIt == '':
            print('Tolerance = %s' %(tolerance))
            simulation.minimizeEnergy(tolerance = (tolerance*kilojoule/mole))
        if maxIt != '' and tolerance == '':
            print('Iterations = %s' %(maxIt))
            simulation.minimizeEnergy(maxIterations = int(maxIt))
        if maxIt != '' and tolerance != '':
            print('Tolerance = %s \nIterations = %s' %(tolerance, maxIt))
            simulation.minimizeEnergy(tolerance = (tolerance*kilojoule/mole), maxIterations = maxIt)
        print('Energy minimized') 
        return integrator, simulation
    
        
    def RunSimulation(self, integrator, simulation):
        if self.nsteps/self.conformations/10 < 250:
            rep = 250
        else:
            rep = int(self.nsteps/self.conformations/10)
            
        simulation.reporters.append(PDBReporter(self.address_output_pdb_MD, self.nsteps/self.conformations))
        simulation.reporters.append(StateDataReporter(self.address_output_states, rep, step=True, time=True,
                                                    potentialEnergy=True, temperature=True))
        simulation.reporters.append(StateDataReporter(stdout, 500, step=True, time=True, temperature=True))
        print('Starting simulation...')
        for i in self.Temperature:
            integrator.setTemperature(i*kelvin)
            simulation.step(self.nsteps/np.size(self.Temperature))
        print('Simulation done')
        
    def EnergyCalculation(self, simulation):

        state = simulation.context.getState(getEnergy=True)
        print(state.getPotentialEnergy())
        
        Energy = state.getPotentialEnergy()
        return Energy
            
def MolecularDynamicsSimulation(S, Sim, seed = None):
    MD = MolecularDynamics(S, Sim, seed)
    pdb, system, modeller = MD.SystemPreparation()
    integrator, simulation = MD.SimulationPreparation(modeller, system)
    MD.RunSimulation(integrator, simulation)
    
    
    
