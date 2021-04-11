import os
from ase.eos import EquationOfState
from ase.optimize import BFGS
from gpaw import GPAW
from ase.io.trajectory import Trajectory
from ase.io import read,write
from fractions import Fraction
import numpy as np
from ase.dft.bee import BEEFEnsemble

def optimize_bulk(atoms,step=0.05,fmax=0.01,location=' ',extname=' '):
    cell=atoms.get_cell()
    name=atoms.get_chemical_formula(mode='hill')
    vol=atoms.get_volume()
    volumes=[]
    energies=[]
    ## TO-DO: how to not output error when making directory
    #os.makedirs(location)
    #os.makedirs(location+'/'+'eos_fit')
    for x in np.linspace(1-2*step,1+2*step,5):
        atoms.set_cell(cell*x,scale_atoms=True)
        atoms.calc.set(txt=location+'/'+'eos_fit'+'/'+name+'_'+str(x)+'-'+extname+'.txt')
        energies.append(atoms.get_potential_energy())
        volumes.append(atoms.get_volume())
    eos=EquationOfState(volumes,energies,eos='birchmurnaghan')
    v0,e0,B=eos.fit()
    x0=(v0/vol)**Fraction('1/3')
    atoms.set_cell(x0*cell,scale_atoms=True)
    atoms.calc.set(txt=location+'/'+name+'_'+str(x0)+'-'+extname+'.txt')
    atoms.calc.attach(atoms.calc.write,5,location+'/'+name+'_'+str(x0)+'-'+extname+'.gpw')
    dyn=BFGS(atoms=atoms,trajectory=location+'/'+name+'_'+str(x0)+'-'+extname+'.traj',
            logfile=location+'/'+name+'_'+str(x0)+'-'+extname+'.log')
    dyn.run(fmax=fmax)
    atoms.calc.write(location+'/'+name+'_'+str(x0)+'-'+extname+'.gpw')
    ## TO-DO: add ensemble energies to file

def surf_relax(surf, name, fmax=0.01, maxstep=0.04, replay_traj=None):
    #calc = surf.calc
    #name = surf.get_chemical_formula(mode='hill')
    gpwname=name+'/'+'slab'
    surf.calc.set(txt=gpwname+'.txt')
    #atoms.calc.attach(atoms.calc.write, 5, gpwname+'.gpw')
    dyn=BFGS(atoms=surf,trajectory=gpwname+'.traj',logfile = gpwname+'.log',restart=gpwname+'qn.pckl',maxstep=maxstep)
    if(replay_traj):
        print('Replaying trajctory file: {}\n'.format(replay_traj))
        dyn.replay_trajectory(replay_traj)
    dyn.run(fmax=fmax)
    surf.calc.write(gpwname+'.gpw')
    #Writing ensemble energies to file 
    if surf.calc.get_xc_functional()=='BEEF-vdW':
        ens = BEEFEnsemble(surf.calc)
        ens_material = ens.get_ensemble_energies()
        with open(str(gpwname)+'_Ensemble_Energies.txt','w+') as file:
            file.write(gpwname+'\n')
            for i in range(len(ens_material)):
                file.write(str(ens_material[i])+'\n')