import os
from ase.parallel import paropen, parprint, world
from ase.db import connect
from ase.io import read
from glob import glob
import numpy as np
from numpy.lib.function_base import _diff_dispatcher
from gpaw import *

# def bulk_calc_conv(element,gpaw_calc,
#                     rela_tol=15*10**(-3), #eV/atom
#                     init_magmom=0,
#                     temp_print=True,
#                     solver_step=0.05,
#                     solver_fmax=0.05,
#                     restart = False):
#     # generate report
#     target_dir='results/'+element+'/'+'bulk'+'/'
#     rep_location=(target_dir+'results_report.txt')
#     calc_dict=gpaw_calc.__dict__['parameters']
#     if world.rank==0 and os.path.isfile(rep_location):
#         os.remove(rep_location)
#     write_report(calc_dict,rep_location,init_magmom,rela_tol)

#     # convergence test 

#     ## h size 
#     # db_h=connect(target_dir+'h_converge.db')
#     # iters=len(db_h)
#     ### restart 
#     if restart:
#         gpw_files_dir=glob(target_dir+'results_h/'+'*.gpw')
#         gpw_files_name=[name.split('/')[-1] for name in gpw_files_dir]
#         if len(gpw_files_name) < 3:
#             updated_gpw=np.sort(gpw_files_name)[0]
#             atoms, calc = restart(updated_gpw)
#         else:
            


#     ## kpts
#     ### restart

#     ## sw
#     ### restart
    
    

    
#     return 'yes'

def bulk_builder(element):
    location='orig_cif_data'+'/'+element+'.cif'
    atoms=read(location)
    return atoms

class bulk_calc_conv:
    def __init__(self,element,gpaw_calc,rela_tol,init_magmom,solver_step,solver_fmax,restart_calc):

        # generate report
        self.target_dir='results/'+element+'/'+'bulk'+'/'
        self.rep_location=(self.target_dir+'results_report.txt')
        self.gpaw_calc=gpaw_calc
        self.calc_dict=gpaw_calc.__dict__['parameters']
        self.rela_tol = rela_tol
        self.init_magmom = init_magmom
        self.initialize_report()
        self.element=element

        # convergence test 

        ## h size 
        param='h'
        ### restart 
        if restart_calc:
            gpw_files_dir=glob(self.target_dir+'results_h/'+'*.gpw')
            gpw_files_name=[name.split('/')[-1] for name in gpw_files_dir]
            h_ls=[float(i.split('-')[-1][:-4]) for i in gpw_files_name]
            descend_order=np.argsort(h_ls)[::-1]
            self.descend_gpw_files_dir=[gpw_files_dir[i] for i in descend_order]
            if len(self.descend_gpw_files_dir) < 3:
                self.restart_calculation(param)
            else: 
                self.convergence_update(param)
                diff_primary=max(self.energies_diff_mat[-1,0],self.energies_diff_mat[-1,2])
                diff_second=self.energies_diff_mat[-1,0]
            iters=len(h_ls)
        else:
            h_ls=[]
        self.convergence_loop(param,h_ls,diff_primary,diff_second,iters)
        
        

                
    def convergence_loop(self,param,param_ls,diff_p=100,diff_s=100,iters=0):
        atoms=bulk_builder(self.element)
        if self.calc_dict['spinpol']:
            atoms.set_initial_magnetic_moments(self.init_magmom*np.ones(len(atoms)))
        if len(param_ls)>0:
            param_val=param_ls[-1]
            if param == 'h':
                self.gpaw_calc.__dict__['parameters'][param]=np.round(param_val-0.02,decimals=2)
            if param == 'kpts':
                self.gpaw_calc.__dict__['parameters'][param]=np.round(param_val-0.02,decimals=2) ## TO-DO
        atoms.set_calculator(self.gpaw_calc)
        opt.optimize_bulk()

    def convergence_update_report(self,param,param_ls):
        f = paropen(self.rep_location,'a')
        for i in self.energies_diff_mat.shape[0]:
            parprint('Optimizing parameter: '+param,file=f)
            parprint('\t'+'1st: '+str(param_ls[i])+' 2nd: '+str(param_ls[i+1])+' 3rd: '+str(param_ls[i+2]),file=f)
            parprint('\t'+'2nd-1st'+'\t'+'3rd-2nd'+'\t'+'3rd-1st',file=f)
            parprint('\t'+str(self.energies_diff_mat[i,:])+'eV')
        f.close()

    def convergence_update(self,param):
        energies=[]
        param_ls=[]
        for gpw_dir in self.descend_gpw_files_dir:
            atoms, calc = restart(gpw_dir)
            param_ls.append(calc.__dict__['parameters'][param])
            energies.append(atoms.get_potential_energy())
        energies_mat = np.zeros(((len(energies)-3)+1,2))
        for i in range((len(energies)-3)+1):
            energies_mat[i,:]=energies[i:i+2]
        energies_mat_rep = (np.concatenate((energies_mat,energies_mat),axis=1))[:,1:3]
        self.energies_diff_mat=np.round(np.abs(energies_mat-energies_mat_rep),decimals=4)
        self.convergence_update_report(param,param_ls)

    def restart_calculation(self,param):
        updated_gpw=self.descend_gpw_files_dir[-1]
        atoms, calc = restart(updated_gpw)
        f = paropen(self.rep_location,'a')
        parprint('Restarting '+param+' size calculation...',file=f)
        parprint('\t'+'h: '+str(calc.__dict__['parameters'][param]),file=f)
        f.close()


    def initialize_report(self):
        if world.rank==0 and os.path.isfile(self.rep_location):
            os.remove(self.rep_location)
        f = paropen(self.rep_location,'a')
        parprint('Initial Parameters:', file=f)
        parprint('\t'+'xc: '+self.calc_dict['xc'],file=f)
        parprint('\t'+'h: '+str(self.calc_dict['h']),file=f)
        parprint('\t'+'kpts: '+str(self.calc_dict['kpts']),file=f)
        parprint('\t'+'sw: '+str(self.calc_dict['occupations']),file=f)
        parprint('\t'+'spin polarized: '+str(self.calc_dict['spinpol']),file=f)
        if self.calc_dict['spinpol']:
            parprint('\t'+'magmom: '+str(self.init_magmom),file=f)
        parprint('\t'+'rela_tol: '+str(self.rela_tol)+'eV',file=f)
        parprint('\n',file=f)
        f.close()
    


