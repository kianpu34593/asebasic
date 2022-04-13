import os
from typing import List
from typing import Dict
from typing import Any

import numpy as np

from ase.parallel import paropen, parprint, world

def initialize_report(report_path: str,
                    calculator_parameters: Dict[str, Any],
                    compute_mode: str,
                    **kwargs,
                    #relative_tolerance: float = 0, #eV/atom
                    ):
    """
    Generate report and print initial parameters for computation.

    Parameters
    ----------

    report_path (REQUIRED):
        Path to the report.
    
    calculator_parameters (REQUIRED):
        Paramters of the calculator.

    compute_mode (REQUIRED):
        The type of the computation. Two types are available: `single_compute` or `convergence_test`.
    
    relative_tolerance:
        Tolerance for convergence test. Need to specify when `compute_mode = 'convergence_test'`.
    """
    defaultkwargs = {'relative_tolerance': 0.015} #eV/atom
    tolerance_dict = {**defaultkwargs, **kwargs}

    if world.rank==0 and os.path.isfile(report_path):
        os.remove(report_path)
    f = paropen(report_path,'a')
    parprint('Computation Type: '+str(compute_mode), file=f)
    parprint('Calculator Parameters:', file=f)
    parprint('\t'+'xc: '+calculator_parameters['xc'],file=f)
    parprint('\t'+'h: '+str(calculator_parameters['h']),file=f)
    parprint('\t'+'kpts: '+str(calculator_parameters['kpts']),file=f)
    parprint('\t'+'sw: '+str(calculator_parameters['occupations']),file=f)
    parprint('\t'+'spin polarized: '+str(calculator_parameters['spinpol']),file=f)
    if compute_mode == 'convergence_test':
        parprint('\t'+'convergence tolerance: '+str(tolerance_dict['relative_tolerance'])+'eV/atom',file=f)
    parprint(' \n',file=f)
    f.close()

def write_message_in_report(report_path: str,
                            message: str):
    """
    A helper function to write a message in the report.

    Parameters
    ----------

    report_path (REQUIRED):
        Path to the report.

    message (REQUIRED):
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
        Parameters of the calculator.
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

def convergence_update_report(parameter: str,
                            parameter_converge_dict: Dict[str, List],
                            report_path: str,
                            energy_difference_array: np.ndarray,
                            ):
    """
    Print out convergence update report.

    Parameters
    ----------
    
    parameter (REQUIRED):
        Parameter to do convergence test.

    parameter_converge_dict (REQUIRED):
        Parameter convergence dictionary.

    report_path (REQUIRED):
        Path to the report.
    
    energy_difference_array (REQUIRED):
        Convergence energy difference arrray.
    """
    f = paropen(report_path,'a')
    parprint('Optimizing parameter: '+parameter,file=f)
    parameter_val_str_1st='1st: '+str(parameter_converge_dict['parameter_converge_lst'][-3])
    parameter_val_str_2nd=' 2nd: '+str(parameter_converge_dict['parameter_converge_lst'][-2])
    parameter_val_str_3rd=' 3rd: '+str(parameter_converge_dict['parameter_converge_lst'][-1])
    parameter_val_str=parameter_val_str_1st+parameter_val_str_2nd+parameter_val_str_3rd
    parprint('\t'+parameter_val_str,file=f)
    divider_str='-'
    parprint('\t'+divider_str*len(parameter_val_str),file=f)
    substrat_str='| '+'2nd-1st'+' | '+'3rd-2nd'+' | '+'3rd-1st'+' |'
    parprint('\t'+substrat_str,file=f)
    energies_str='\t'+'| '
    for i in range(3):
        energies_str+=str(energy_difference_array[i])+'  '+'|'+' '
    energies_str+='eV/atom'
    parprint(energies_str,file=f)
    parprint(' ',file=f)
    f.close()

def restart_report(parameter: str,
                calculator_parameters: Dict[str, Any],
                report_path: str,
                ):
    """
    Print out restart report.

    Parameters
    ----------
    
    parameter (REQUIRED):
        Parameter to do convergence test.

    calculator_parameters (REQUIRED):
        Parameters of the calculator.

    report_path (REQUIRED):
        Path to the report.
    """
    f = paropen(report_path,'a')
    parprint('Restarting '+parameter+' convergence test...',file=f)
    parameter_i=parameter.split('_')[0]
    parprint('\t'+'Last computation:'+'\t'+parameter+'='+str(calculator_parameters[parameter_i]),file=f)
    parprint(' ',file=f)
    f.close()