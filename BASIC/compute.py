## this is a wrapper function file for bulk, adsorption and surface calculation ##
import os

from typing import Tuple
from typing import Union
from typing import Dict
from typing import Any

from gpaw import restart

from ase.io import read
from ase.db import connect
from ase.parallel import barrier,world
from ase.constraints import FixAtoms

import BASIC.optimize as opt
import BASIC.message as msg
import BASIC.utils as ut 

# class surface_single_compute:
#     pass


def slab_compute(element:str,
                calculator_setting,
                computation_setting: Dict[str, Any],
                converge_parameter: Tuple[str, Union[float,str,int]],
                restart_calculation: bool,
                target_dir: str = None):

    defaultkwargs = {'solver_maxstep': 0.05, 'solver_fmax':0.03}
    optimizer_setting = {**defaultkwargs, **kwargs}

    miller_plane = computation_setting['miller_plane']
    shift = computation_setting['shift']
    order = computation_setting['order_val']
    full_name = '_'.join([miller_plane, shift, order])
    layer = str(computation_setting['layer'])
    fix_layer = str(computation_setting['fix_layer'])
    fix_mode = computation_setting['fix_mode']
    surface_energy_calculation_mode = computation_setting['surface_energy_calculation_mode']
    
    if surface_energy_calculation_mode == 'linear_fit' and 'fitted_bulk_energy' not in computation_setting.keys():
        raise RuntimeError ('Fitted bulk energy not specified. linear_fit is not supported.')
    # generate report
    if converge_parameter[0] == "single_compute":
        
        target_dir=os.path.join('results',element,'surf',full_name)
        converge_parameter_dir = os.path.join(target_dir,'single_compute')
        compute_dir=os.path.join(converge_parameter_dir, layer)
        if world.rank==0 and not os.path.isdir(layer_dir):
            os.makedirs(layer_dir,exist_ok=True)
        report_path=os.path.join(converge_parameter_dir, layer+'_results_report.txt')

        msg.initialize_report(report_path,calculator_setting.parameters,compute_mode=converge_parameter[0])

    else:
        compute_dir = os.path.join(target_dir,'convergence_test', 'calculator_parameter', layer, converge_parameter[0])
        report_path = os.path.join(target_dir,'convergence_test','calculator_parameter',layer, f"{converge_parameter[0]}_report.txt")


    # layer=self.ascend_all_cif_files_full_path[iters].split('/')[-1].split('.')[0]
    # location=self.target_sub_dir+layer+'x1x1'
    if restart_calculation and os.path.isfile(os.path.join(compute_dir,f"*-{str(converge_parameter[1])}_interm.gpw")):
        slab, calculator_setting = restart(os.path.join(compute_dir,f"*-{str(converge_parameter[1])}_interm.gpw"))
        msg.write_message_in_report(report_path,'restart_calculation == True and .gpw file exist. Restart from previous calculation')
    else:
        slab=read(os.path.join(target_dir, 'input_slab', f"{layer}.traj"))
        slab.pbc = [1,1,0]
        if 'vacuum' in computation_setting.keys():
            vacuum_size = computation_setting['vacuum']
        else:
            vacuum_size = 10
        slab.center(vacuum=vacuum_size,axis=2)
        slab=ut.fix_layer(slab,fix_layer,fix_mode)
        slab.set_calculator(calculator_setting)
    opt.relax_slab(slab,compute_dir,fmax=optimizer_setting['solver_fmax'],maxstep=optimizer_setting['solver_maxstep'])

    #TO-DO finalize single_compute

    return slab

def bulk_compute(
            element: str,
            calculator_setting,
            converge_parameter: Tuple[str, Union[float,str,int]],
            target_dir: str = None,
            **kwargs,
            # eos_step: float = 0.05,
            # solver_maxstep: float = 0.05,
            # solver_fmax: float = 0.03,
            # restart_calculation: bool=True, 
            ):
    """
    Bulk crystal structure computation without convergence.
   
    Parameters
    ----------

    element (REQUIRED):
        Chemical symbols and the materials project id of the bulk structures to be computed. E.g. Cu_mp-30

    calculator_setting (REQUIRED):
        Dictionary of calculator setting from ASE interface (description needs improvement).
    
    converge_parameter:
        Parameter to converge for. For convergence test, available options are ('grid_spacing', 0.16) or ('kdensity', 3.5); 
        for single compute, available options is ('single_compute', '')
        
    eos_step: 
        Equation of state step size for lattice optimization. 
        If not specified, default is 0.05.

    solver_maxstep:
        Maxstep for BFGS solver.
        If not specified, default is 0.05.
    
    solver_fmax:
        Maximum force for BFGS solver.
        If not specified, default is 0.03.
    """
    defaultkwargs = {'eos_step': 0.05,'solver_maxstep': 0.05, 'solver_fmax':0.03}
    optimizer_setting = {**defaultkwargs, **kwargs}
    # generate report
    if converge_parameter[0] == "single_compute":
        target_dir=os.path.join('results',element,'bulk',converge_parameter[0])
        if world.rank==0 and not os.path.isdir(target_dir):
            os.makedirs(target_dir,exist_ok=True)
        report_path=os.path.join(target_dir, 'results_report.txt')
        msg.initialize_report(report_path,calculator_setting.parameters,compute_mode=converge_parameter[0])
    
    eos_fit_dir=os.path.join(target_dir,'eos_fit')
    if world.rank==0 and not os.path.isdir(eos_fit_dir):
        os.makedirs(eos_fit_dir,exist_ok=True)
    barrier()

    # lattice optimization and relax
    traj_file_path=os.path.join('orig_cif_data',element,'input.traj')
    atoms=read(traj_file_path)
    atoms.set_calculator(calculator_setting)
    opt.optimize_bulk(atoms, bulk_path = target_dir, name_extension = str(converge_parameter[1]), eos_step = optimizer_setting['eos_step'], fmax = optimizer_setting['solver_fmax'], maxstep = optimizer_setting['solver_maxstep'])

    #finalize #TO-DO need some rethink on this
    if converge_parameter[0] == 'single_compute':
        db_final=connect('final_database'+'/'+'bulk_single.db')
        id=db_final.reserve(full_name=element)
        if id is None:
            id=db_final.get(full_name=element).id
            db_final.update(id=id,atoms=atoms)
        else:
            db_final.write(atoms,id=id,full_name=element)

        msg.write_message_in_report(report_path, message='single_compute complete!')
    
    return atoms



# class adsorption_compute:
#     pass

