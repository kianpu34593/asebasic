import os
from ase.eos import EquationOfState
from ase.optimize import BFGS
from gpaw import GPAW
from ase.io.trajectory import Trajectory
from ase.io import read,write
from fractions import Fraction
import numpy as np
from ase.dft.bee import BEEFEnsemble
from ase.parallel import parprint,world,barrier


def optimize_bulk(atoms, #how to write type for ase atom object
                bulk_path: str,
                name_extension: str,
                eos_step: float,
                fmax: float,
                maxstep: float,
                ):
    """
    Optimize for the bulk structure

    Parameters
    ----------

    atoms: Atoms object
        The Atoms object to relax.
        
    bulk_path:
        Path to the bulk directory of the materials computing
    
    name_extension:
        Name extension for the calculation. Possible values are grid spacing (h) and kpts density (kdens)
    
    eos_step: 
        Equation of state step size for lattice optimization. 

    maxstep:
        Maxstep for BFGS solver.
    
    fmax:
        Maximum force for BFGS solver.
    """
    cell=atoms.get_cell()
    name=atoms.get_chemical_formula(mode='hill')
    volume=atoms.get_volume()
    
    volume_lst=[]
    energy_lst=[]
    if atoms.calc.parameters['spinpol'] in [None, False]:
        atoms.set_initial_magnetic_moments()
    for scale in np.linspace(1-2*eos_step,1+2*eos_step,5):
        atoms.set_cell(cell*scale,scale_atoms=True)
        path_to_txt=os.path.join(bulk_path, 'eos_fit', f"{name}_{str(np.round(scale,decimals=2))}-{name_extension}.txt")
        atoms.calc.set(txt=path_to_txt)
        energy_lst.append(atoms.get_potential_energy())
        volume_lst.append(atoms.get_volume())
    #plot curve
    eos=EquationOfState(volume_lst,energy_lst,eos='birchmurnaghan')
    v0=eos.fit()[0]
    eos.plot(os.path.join(bulk_path,'eos_fit.png'))
    scale=(v0/volume)**Fraction('1/3')
    atoms.set_cell(cell*scale,scale_atoms=True)

    path_to_txt=os.path.join(bulk_path,f"{name}-{name_extension}.txt")
    path_to_traj=os.path.join(bulk_path,f"{name}-{name_extension}.traj")
    path_to_log=os.path.join(bulk_path,f"{name}-{name_extension}.log")
    atoms.calc.set(txt=path_to_txt)
    
    dyn=BFGS(atoms=atoms,trajectory=path_to_traj, logfile=path_to_log,maxstep=maxstep)
    dyn.run(fmax=fmax)
    atoms.calc.write(os.path.join(bulk_path,f"{name}-{name_extension}_finish.gpw"))
    atoms.write(os.path.join(bulk_path,f"{name}-{name_extension}_finish.traj"))
    ## TO-DO: add ensemble energies to file

def restart_calculation_check(path:str,
                        name: str,
                        hist_name: str):
    """
    Check if .traj history exist and extend it with current .traj

    Paramters
    ---------
    path(REQUIRED): 
        Path to the .traj directory

    name(REQUIRED):
        Name of the current .traj file
    
    hist_name(REQUIRED):
        Name of the history .traj file
    """
    def _check_file_exists(filename):
        """Check if file exists and is not empty"""
        if os.path.isfile(filename):
            return os.path.getsize(filename) > 0
        else:
            return False
    path_to_slab_traj = os.path.join(path,f"{name}.traj")
    path_to_hist_slab_traj = os.path.join(path,f"{hist_name}.traj")
    # check if it is a restart
    barrier()
    if _check_file_exists(path_to_slab_traj):
        latest = read(path_to_slab_traj, index=":")
        # check if already restarted previously and extend history if needed
        if not (_check_file_exists(path_to_slab_traj)):
            barrier()
            write(path_to_hist_slab_traj,latest)
        else:
            hist = read(path_to_hist_slab_traj, index=":")
            hist.extend(latest)
            write(path_to_hist_slab_traj,hist)
    return _check_file_exists(path_to_hist_slab_traj)

def relax_slab(atoms, 
                slab_dir, 
                name_extension,
                restart_calculation,
                maxstep,
                fmax, 
                ):
    """
    Optimize for the slab structure

    Parameters
    ----------

    atoms: Atoms object
        The Atoms object to relax.
        
    slab_dir:
        Path to the slab directory of the materials computing
    
    name_extension:
        Name extension for the calculation. Possible values are grid spacing (h) and kpts density (kdens)
    
    restart_calculation: 
        Whether to restart with the previous calculaiton. 

    maxstep:
        Maxstep for BFGS solver.
    
    fmax:
        Maximum force for BFGS solver.
    """
    
    name=atoms.get_chemical_formula(mode='hill')
    # slab_name=os.path.join(slab_path, name)
    # slab_hist_name=f"{slab_name}_history"
    if atoms.calc.parameters['spinpol'] in [None, False]:
        atoms.set_initial_magnetic_moments()
    atoms.calc.__dict__['observers']=[]
    atoms.calc.attach(atoms.calc.write, 10, os.path.join(slab_dir, f"{name}-{name_extension}_interm.gpw"))

    hist_exist=restart_calculation_check(slab_dir, name, f"{name}-{name_extension}_history")


    path_to_txt=os.path.join(slab_dir,f"{name}-{name_extension}.txt")
    path_to_traj=os.path.join(slab_dir,f"{name}-{name_extension}.traj")
    path_to_log=os.path.join(slab_dir,f"{name}-{name_extension}.log")
    atoms.calc.set(txt=path_to_txt)
    dyn=BFGS(atoms=atoms,
            trajectory=path_to_traj,
            logfile = path_to_log,
            maxstep=maxstep)
    

    # if history exists, read in hessian
    if hist_exist and restart_calculation:
        dyn.replay_trajectory(os.path.join(slab_dir, f"{name}-{name_extension}_history.traj"))

    # optimize
    dyn.run(fmax=fmax)
    atoms.calc.write(os.path.join(slab_dir, f"{name}-{name_extension}_finish.gpw"))
    atoms.write(os.path.join(slab_dir, f"{name}-{name_extension}_finish.traj"))


# class OccasionalMagCalc():
#     def __init__(self,atoms):
#         self.iter=0
#         self.atoms=atoms
    
#     def magmom_calc(self):
#         self.atoms.get_magnetic_moment()

