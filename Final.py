import os

from MolecularDynamics.MD_Simulation import *
from ANM.ANM_Simulation import *
from MC.MC_simulation import *
from StructuralAnalysis.SimulationAnalysis import *
import json

path = sys.argv[1]
with open(path) as file:
    Information = json.load(file)
    

if Information["Type"] == 'Simulation':
    if Information["Name"] == "MD":
        print("Initializing Molecular dynamics...")
        MolecularDynamicsSimulation(Information, 'MD')
    if Information["Name"] == "ANM":
        print("Initializing ANM...")
        ANMSimulation(Information)        
    if Information["Name"] == "P":
        print("Initializing MC...")
        Metropolis(Information, 'P')

if Information["Type"] == 'Analysis':
    print("Initializing Analysis...")
    simulationAnalysis(Information)

