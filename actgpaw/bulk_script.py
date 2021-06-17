import os
from re import S
from ase.parallel import paropen, parprint, world
from ase.db import connect
from ase.io import read
from glob import glob
import numpy as np
from numpy.lib.function_base import _diff_dispatcher
from gpaw import *
import actgpaw.optimizer as opt
import sys
from ase.calculators.calculator import kptdensity2monkhorstpack as kdens2mp
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
        self.solver_step=solver_step
        self.solver_fmax=solver_fmax

        # convergence test 

        ## h size 
        param='h'
        ### restart 
        if restart_calc and len(glob(self.target_dir+'results_'+param+'/'+'*.gpw'))>0:
            descend_param_ls,descend_gpw_files_dir=self.gather_gpw_file(param)
            if len(descend_gpw_files_dir) < 3:
                self.restart_report(param,descend_gpw_files_dir[-1])
            else: 
                for i in range((len(descend_param_ls)-3)+1):
                    print(i)
                    self.convergence_update(param,i,descend_gpw_files_dir)
                    diff_primary=max(self.energies_diff_mat[0],self.energies_diff_mat[2])
                    diff_second=self.energies_diff_mat[1]
            self.gpaw_calc.__dict__['parameters'][param]=np.round(descend_param_ls[-1]-0.02,decimals=2)
            self.calc_dict=self.gpaw_calc.__dict__['parameters']
            print('exit')
            sys.exit()
        else:
            descend_param_ls=[]
            diff_primary=100
            diff_second=100
        ### convergence loop
        iters=len(descend_param_ls)
        self.convergence_loop(param,iters,diff_primary,diff_second)
       
        ## kpts size 
        param='kdens'
        ### jump the first calculation
        descend_gpw_files_dir=self.gather_gpw_file('h')[1]
        atoms, calc = restart(descend_gpw_files_dir[-3])
        self.gpaw_calc=calc
        self.calc_dict=self.gpaw_calc.__dict__['parameters']
        param_val=self.calc_dict['kpts']['density']
        opt.optimize_bulk(atoms,
                    step=self.solver_step,fmax=self.solver_fmax,
                    location=self.target_dir+'results_'+param,
                    extname=param_val)

        ### restart 
        if restart_calc and len(glob(self.target_dir+'results_'+param+'/'+'*.gpw'))>0:
            descend_param_ls,descend_gpw_files_dir=self.gather_gpw_file(param)
            if len(descend_gpw_files_dir) < 3:
                self.restart_report(param,descend_gpw_files_dir[0])
            else: 
                for i in range((len(descend_param_ls)-3)+1):
                    self.convergence_update(param,i,descend_gpw_files_dir)
                    diff_primary=max(self.energies_diff_mat[0],self.energies_diff_mat[2])
                    diff_second=self.energies_diff_mat[1]
                # atoms,calc=restart(descend_gpw_files_dir[0])
                atoms=bulk_builder(self.element)
                kpts=kdens2mp(atoms,kptdensity=descend_param_ls[0])
                new_kpts=kpts.copy()
                new_kdens=descend_param_ls[0].copy()
                while np.mean(kpts)==np.mean(new_kpts):
                    new_kdens+=0.1
                    new_kpts=kdens2mp(atoms,kptdensity=new_kdens)
                new_kdens_dict={'density':new_kdens,'even':True}
            self.gpaw_calc.__dict__['parameters']['kpts']=new_kdens_dict
            self.calc_dict=self.gpaw_calc.__dict__['parameters']
        else:
            descend_param_ls=[]
            diff_primary=100
            diff_second=100
        ### convergence loop
        iters=len(descend_param_ls)
        self.convergence_loop(param,iters,diff_primary,diff_second)

        #finalize
        descend_gpw_files_dir=self.gather_gpw_file(param)[1]
        final_atoms, calc = restart(descend_gpw_files_dir[-3])
        self.gpaw_calc=calc
        self.calc_dict=self.gpaw_calc.__dict__['parameters']
        if self.calc_dict['spinpol']:
            self.final_magmom=final_atoms.get_magnetic_moments()
        db_final=connect('final_database'+'/'+'bulk.db')
        id=db_final.reserve(name=element)
        if id is None:
            id=db_final.get(name=element).id
            db_final.update(id=id,atoms=final_atoms,name=element,
                            restart_dir=descend_gpw_files_dir[-3])
        else:
            db_final.write(final_atoms,id=id,name=element,
                            restart_dir=descend_gpw_files_dir[-3])
        self.final_report()

        


    def convergence_loop(self,param,iters,diff_p,diff_s):
        while (diff_p>self.rela_tol or diff_s>self.rela_tol) and iters <= 6:
            atoms=bulk_builder(self.element)
            if self.calc_dict['spinpol']:
                atoms.set_initial_magnetic_moments(self.init_magmom*np.ones(len(atoms)))
            atoms.set_calculator(self.gpaw_calc)
            if param == 'h':
                param_val=self.calc_dict[param]
            elif param == 'kdens':
                param_val=self.calc_dict['kpts']['density']
            opt.optimize_bulk(atoms,
                                step=self.solver_step,fmax=self.solver_fmax,
                                location=self.target_dir+'results_'+param,
                                extname=param_val)
            #convergence update
            descend_param_ls,descend_gpw_files_dir=self.gather_gpw_file(param)
            if iters>2:
                iter=iters-3
                self.convergence_update(param,iter,descend_gpw_files_dir)
                diff_p=max(self.energies_diff_mat[0],self.energies_diff_mat[2])
                diff_s=self.energies_diff_mat[1]
            #update param
            if param == 'h':
                self.gpaw_calc.__dict__['parameters'][param]=np.round(param_val-0.02,decimals=2)
            elif param == 'kdens':
                atoms=bulk_builder(self.element)
                kpts=kdens2mp(atoms,kptdensity=descend_param_ls[0])
                new_kpts=kpts.copy()
                new_kdens=descend_param_ls[0].copy()
                while np.mean(kpts)==np.mean(new_kpts):
                    new_kdens+=0.1
                    new_kpts=kdens2mp(atoms,kptdensity=new_kdens)
                new_kdens_dict={'density':new_kdens,'even':True}
                self.gpaw_calc.__dict__['parameters']['kpts']=new_kdens_dict
            self.calc_dict=self.gpaw_calc.__dict__['parameters']
            iters=len(descend_param_ls)
        #check iteration
        self.check_convergence(diff_p,diff_s,iters,param)
    
    def check_convergence(self,diff_p,diff_s,iters,param):
        if iters>=6:
            if diff_p>self.rela_tol or diff_s>self.rela_tol:
                f=paropen(self.rep_location,'a')
                parprint("WARNING: Max iterations reached! "+param+" convergence test failed.",file=f)
                parprint("Computation Suspended!",file=f)
                parprint(' ',file=f)
                f.close()
                sys.exit()
        else:
            f=paropen(self.rep_location,'a')
            parprint(param+" convergence test success!",file=f)
            parprint('\n',file=f)
            f.close() 


    def gather_gpw_file(self,param):
        gpw_files_dir=glob(self.target_dir+'results_'+param+'/'+'*.gpw')
        gpw_files_name=[name.split('/')[-1] for name in gpw_files_dir]
        param_ls=[float(i.split('-')[-1][:-4]) for i in gpw_files_name]
        descend_order=np.argsort(param_ls)[::-1]
        descend_gpw_files_dir=[gpw_files_dir[i] for i in descend_order]
        descend_param_ls=np.sort(param_ls)
        return descend_param_ls,descend_gpw_files_dir

    def convergence_update(self,param,iter,gpw_files_dir):
        energies=[]
        param_ls=[]
        if param == 'kdens':
            gpw_files_dir=gpw_files_dir[::-1]
        for i in range(iter,iter+3,1):

            atoms, calc = restart(gpw_files_dir[i])
            if param == 'kdens':
                kdens=calc.__dict__['parameters']['kpts']['density']
                param_ls.append(kdens)
            elif param == 'h':
                param_ls.append(calc.__dict__['parameters'][param])
            energies.append(atoms.get_potential_energy()/len(atoms)) #eV/atom
        energies_mat = np.array(energies)
        energies_mat_rep = np.array((energies+energies)[1:4])
        self.energies_diff_mat=np.round(np.abs(energies_mat-energies_mat_rep),decimals=4)
        print(energies_mat)
        print(energies_mat_rep)
        self.convergence_update_report(param,param_ls)

    def convergence_update_report(self,param,param_ls):
        f = paropen(self.rep_location,'a')
        parprint('Optimizing parameter: '+param,file=f)
        param_val_str='1st: '+str(param_ls[0])+' 2nd: '+str(param_ls[1])+' 3rd: '+str(param_ls[2])
        parprint('\t'+param_val_str,file=f)
        divider_str='-'
        parprint('\t'+divider_str*len(param_val_str),file=f)
        substrat_str='2nd-1st'+' | '+'3rd-2nd'+' | '+'3rd-1st'
        parprint('\t'+substrat_str,file=f)
        energies_str='\t'
        for i in range(3):
            energies_str+=str(self.energies_diff_mat[i])+'   '+'|'+' '
        energies_str+='\t'+'eV/atom'
        parprint(energies_str,file=f)
        parprint(' ',file=f)
        f.close()

    def restart_report(self,param,updated_gpw):
        calc = restart(updated_gpw)[1]
        f = paropen(self.rep_location,'a')
        parprint('Restarting '+param+' convergence test...',file=f)
        parprint('\t'+'Last computation:'+'\t'+param+'='+str(calc.__dict__['parameters'][param]),file=f)
        parprint(' ',file=f)
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
        parprint('\t'+'convergence tolerance: '+str(self.rela_tol)+'eV/atom',file=f)
        parprint(' ',file=f)
        f.close()
    
    def final_report(self):
        f = paropen(self.rep_location,'a')
        parprint('Final Parameters:', file=f)
        parprint('\t'+'xc: '+self.calc_dict['xc'],file=f)
        parprint('\t'+'h: '+str(self.calc_dict['h']),file=f)
        parprint('\t'+'kpts: '+str(self.calc_dict['kpts']),file=f)
        parprint('\t'+'sw: '+str(self.calc_dict['occupations']),file=f)
        parprint('\t'+'spin polarized: '+str(self.calc_dict['spinpol']),file=f)
        if self.calc_dict['spinpol']:
            parprint('\t'+'magmom: '+str(self.final_magmom),file=f)
        parprint(' ',file=f)
        f.close()
    


