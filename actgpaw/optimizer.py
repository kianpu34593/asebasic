import os
from ase.eos import EquationOfState
from ase.optimize import BFGS
from gpaw import GPAW
from ase.io.trajectory import Trajectory
from ase.io import read,write
from fractions import Fraction
import numpy as np
from ase.dft.bee import BEEFEnsemble

def optimize_bulk(atoms,step=0.05,fmax=0.01,location='',extname=''):
    cell=atoms.get_cell()
    name=atoms.get_chemical_formula(mode='hill')
    vol=atoms.get_volume()
    volumes=[]
    energies=[]
    for x in np.linspace(1-2*step,1+2*step,5):
        atoms.set_cell(cell*x,scale_atoms=True)
        atoms.calc.set(txt=location+'/'+'eos_fit'+'/'+name+'_'+str(np.round(x,decimals=2))+'-'+str(extname)+'.txt')
        energies.append(atoms.get_potential_energy())
        volumes.append(atoms.get_volume())
    eos=EquationOfState(volumes,energies,eos='birchmurnaghan')
    v0=eos.fit()[0]
    x0=(v0/vol)**Fraction('1/3')
    atoms.set_cell(x0*cell,scale_atoms=True)
    file_name=location+'/'+name+'-'+str(extname)
    atoms.calc.set(txt=file_name+'.txt')
    dyn=BFGS(atoms=atoms,trajectory=file_name+'.traj',
            logfile=file_name+'.log') ## TO-DO: add maxstep control
    dyn.run(fmax=fmax)
    atoms.calc.write(file_name+'.gpw')
    ## TO-DO: add ensemble energies to file

# def relax(atoms, name, fmax=0.01, maxstep=0.04):
#     gpwname=name+'/'+'slab'
#     atoms.calc.set(txt=gpwname+'.txt')
#     atoms.calc.attache(atoms.calc.write, 10, 'interm.gpaw')
#     dyn=BFGS(atoms=atoms,trajectory=gpwname+'.traj',
#                 logfile = gpwname+'.log',maxstep=maxstep)
#     dyn.run(fmax=fmax)
#     atoms.calc.write(gpwname+'.gpw')
#     # TO-DO: add ensemble energies to file


def relax(atoms, name, fmax=0.01, maxstep=0.04):
    gpwname=name+'/'+'slab'
    atoms.calc.set(txt=gpwname+'.txt')
    atoms.calc.attach(atoms.calc.write, 10, "interm.gpw")

    def _check_file_exists(filename):
        """Check if file exists and is not empty"""
        if os.path.isfile(filename):
            return os.path.getsize(filename) > 0
        else:
            return False

    # check if it is a restart
    if _check_file_exists(gpwname+".traj"):
        latest = read(gpwname+".traj", index=":")
        # check if already restarted previously and extend history if needed
        if _check_file_exists(gpwname+'_history'+".traj"):
            hist = read(gpwname+'_history'+".traj", index=":")
            hist.extend(latest)
            hist.write(gpwname+'_history'+".traj")
        else:
            latest.write(gpwname+'_history'+".traj")
    dyn=BFGS(atoms=atoms,trajectory=gpwname+'.traj',
            logfile = gpwname+'.log',maxstep=maxstep)
    # if history exists, read in hessian
    if _check_file_exists(gpwname+'_history'+".traj"):
        dyn.replay_trajectory(gpwname+'_history'+".traj")
    # optimize
    dyn.run(fmax=fmax)
    atoms.calc.write(gpwname+'.gpw')