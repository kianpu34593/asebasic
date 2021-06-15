import os
from ase.parallel import paropen, parprint, world
from ase.db import connect
from glob import glob
import numpy as np
from gpaw import *
def bulk_calc_conv(element,gpaw_calc,
                    rela_tol=15*10**(-3), #eV/atom
                    init_magmom=0,
                    temp_print=True,
                    solver_step=0.05,
                    solver_fmax=0.05,
                    restart = False):
    # generate report
    target_dir='results/'+element+'/'+'bulk'+'/'
    rep_location=(target_dir+'results_report.txt')
    calc_dict=gpaw_calc.__dict__['parameters']
    if world.rank==0 and os.path.isfile(rep_location):
        os.remove(rep_location)
    write_report(calc_dict,rep_location,init_magmom,rela_tol)

    # convergence test 

    ## h size 
    # db_h=connect(target_dir+'h_converge.db')
    # iters=len(db_h)
    ### restart 
    if restart:
        gpw_files_dir=glob(target_dir+'results_h/'+'*.gpw')
        gpw_files_name=[name.split('/')[-1] for name in gpw_files_dir]
        if len(gpw_files_name) < 3:
            updated_gpw=np.sort(gpw_files_name)[0]
            atoms, calc = restart(updated_gpw)
        else:
            


    ## kpts
    ### restart

    ## sw
    ### restart
    
    

    
    return 'yes'

def write_report(calc_dict,location,init_magmom,rela_tol):
    f = paropen(location,'a')
    parprint('Initial Parameters:', file=f)
    parprint('\t'+'xc: '+calc_dict['xc'],file=f)
    parprint('\t'+'h: '+str(calc_dict['h']),file=f)
    parprint('\t'+'kpts: '+str(calc_dict['kpts']),file=f)
    parprint('\t'+'sw: '+str(calc_dict['occupations']),file=f)
    parprint('\t'+'spin polarized: '+str(calc_dict['spinpol']),file=f)
    if calc_dict['spinpol']:
        parprint('\t'+'magmom: '+str(init_magmom),file=f)
    parprint('\t'+'rela_tol: '+str(rela_tol)+'eV',file=f)
    f.close()
# class bulk_calc_conv(element,
#                     gpaw_calc,
#                     rela_tol=15*10**(-3), #eV/atom
#                     init_magmom=0,
#                     temp_print=True,
#                     solver_step=0.05,
#                     solver_fmax=0.05):
#     def __init__(self,element,gpaw_calc):
#         self.rep_location=('results/'+element+'/'+'bulk'+'/'+'results_report.txt')
#         self.calc_dict=gpaw_calc.__dict__['parameters']
#     def generate_report(self,)

