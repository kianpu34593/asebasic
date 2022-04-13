## this is a wrapper function file for bulk, adsorption and surface calculation ##
import os
from typing import List
from typing import Dict
from typing import Any
from typing import Tuple
from typing import Union

from ase.parallel import paropen, parprint, world
from ase.io import read
from ase.db import connect

import BASIC.optimizer as opt

from glob import glob

# class surface_single_compute:
#     pass




def bulk_single_compute(
            element: str,
            calculator_setting,
            converge_parameter: Tuple(str, Union[float,str,int]),
            **kwargs
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
    defaultKwargs = {'eos_step': 0.05,'solver_maxstep': 0.05, 'solver_fmax':0.03}
    optimizer_setting = {**defaultKwargs, **kwargs}
    # generate report
    if converge_parameter[0] == "single_compute":
        target_dir=os.path.join('results',element,'bulk',converge_parameter[0])
        report_path=os.path.join(target_dir, 'results_report.txt')
        initialize_report(report_path,calculator_setting.parameters,compute_type=converge_parameter[0])

    # lattice optimization and relax
    traj_file_path=os.path.join('orig_cif_data',element,'input.traj')
    atoms=read(traj_file_path)
    atoms.set_calculator(calculator_setting)
    opt.optimize_bulk(atoms, bulk_path = target_dir, name_extension = str(converge_parameter[1]), eos_step = optimizer_setting['eos_step'], fmax = optimizer_setting['solver_fmax'], maxstep = optimizer_setting['solver_maxstep'])

    #finalize #TO-DO need some rethink on this
    if converge_parameter[0] == 'single_compute':
        db_final=connect('final_database'+'/'+'bulk.db')
        id=db_final.reserve(name=element)
        if id is None:
            id=db_final.get(name=element).id
            db_final.update(id=id,atoms=atoms,name=element,converge_test=False,converge_parameter='')
        else:
            db_final.write(atoms,id=id,name=element,converge_test=False,converge_parameter='')

        write_message_in_report(report_path, message='single_compute complete!')
    
    return atoms



class adsorption_compute:
    pass


def initialize_report(report_path: str,
                    calculator_parameters: Dict[str, Any],
                    compute_type: str,
                    relative_tolerance: float = 0, #eV/atom
                    ):
    """
    Generate report and print initial parameters for computation.

    Parameters
    ----------

    report_path (REQUIRED):
        Path to the report.
    
    calculator_parameters (REQUIRED):
        Paramters of the calculator.

    compute_type (REQUIRED):
        The type of the computation. Two types are available: `single_compute` or `convergence_test`.
    
    relative_tolerance:
        Tolerance for convergence test. Need to specify when `compute_type = 'convergence_test'`.

    magnetic_moments (REQUIRED):
        Dictionary with magentic moments to be set for the chemical species with `spinpol=True` (spin polarized computation).
    """
    if world.rank==0 and os.path.isfile(report_path):
        os.remove(report_path)
    f = paropen(report_path,'a')
    parprint('Computation Type: '+str(compute_type), file=f)
    parprint('Calculator Parameters:', file=f)
    parprint('\t'+'xc: '+calculator_parameters['xc'],file=f)
    parprint('\t'+'h: '+str(calculator_parameters['h']),file=f)
    parprint('\t'+'kpts: '+str(calculator_parameters['kpts']),file=f)
    parprint('\t'+'sw: '+str(calculator_parameters['occupations']),file=f)
    parprint('\t'+'spin polarized: '+str(calculator_parameters['spinpol']),file=f)
    if compute_type == 'convergence_test':
        parprint('\t'+'convergence tolerance: '+str(relative_tolerance)+'eV/atom',file=f)
    parprint(' \n',file=f)
    f.close()

def write_message_in_report(report_path: str,
                            message: str):
    """
    A helper function to write a message in the report.

    Parameters
    ----------

    report_path (required):
        Path to the report.

    message:
        The message to be print out in the report.
    """

    f = paropen(report_path,'a')
    parprint('Computation message: '+message,file=f)
    parprint(' \n',file=f)
    f.close()

def final_report(report_path: str,
                calculator_parameters: Dict[str, Any],
                ):
    """
    Print final parameters of the computation.

    Parameters
    ----------

    report_path (REQUIRED):
        Path to the report.
    
    calculator_parameters (REQUIRED):
        Paramters of the calculator.
    """
    f = paropen(report_path,'a')
    parprint('Calculator Parameters:', file=f)
    parprint('\t'+'xc: '+calculator_parameters['xc'],file=f)
    parprint('\t'+'h: '+str(calculator_parameters['h']),file=f)
    parprint('\t'+'kpts: '+str(calculator_parameters['kpts']),file=f)
    parprint('\t'+'sw: '+str(calculator_parameters['occupations']),file=f)
    parprint('\t'+'spin polarized: '+str(calculator_parameters['spinpol']),file=f)
    parprint(' ',file=f)
    f.close()

