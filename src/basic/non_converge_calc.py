from copy import Error
import os
from typing import Type
from ase.parallel import paropen, parprint, world
from ase.db import connect
from ase.io import read
from glob import glob
import numpy as np
from gpaw import restart
import BASIC.optimizer as opt
import sys
from ase.constraints import FixAtoms,FixedLine
import pandas as pd
from BASIC.utils import detect_cluster

def pbc_checker(slab):
    anlges_arg=[angle != 90.0000 for angle in np.round(slab.cell.angles(),decimals=4)[:2]]
    if np.any(anlges_arg):
        slab.pbc=[1,1,1]
    else:
        slab.pbc=[1,1,0]

# def detect_cluster(slab,tol=0.1):
#     n=len(slab)
#     dist_matrix=np.zeros((n, n))
#     slab_c=np.sort(slab.get_positions()[:,2])
#     for i, j in itertools.combinations(list(range(n)), 2):
#         if i != j:
#             cdist = np.abs(slab_c[i] - slab_c[j])
#             dist_matrix[i, j] = cdist
#             dist_matrix[j, i] = cdist
#     condensed_m = squareform(dist_matrix)
#     z = linkage(condensed_m)
#     clusters = fcluster(z, tol, criterion="distance")
#     return slab_c,list(clusters)

def apply_magmom_opt_slab(opt_slab_magmom,ads_slab,adatom=1):
    if adatom == 1:
        magmom_ls=np.append(opt_slab_magmom,0)
    elif adatom == 2:
        magmom_ls=np.append(opt_slab_magmom,0)
        magmom_ls=np.append(magmom_ls,0)
    ads_slab.set_initial_magnetic_moments(magmom_ls)
    return ads_slab

def apply_magmom_manual(ads_slab,magmom_slab,magmom_ads,adatom=1):
    #what if it is alloy?
    magmom_ls=np.ones(len(ads_slab)-1)*magmom_slab
    if adatom == 1:
        magmom_ls=np.append(magmom_ls,magmom_ads)
    ads_slab.set_initial_magnetic_moments(magmom_ls)
    return ads_slab

def get_clean_slab(element, 
                    miller_index,
                    report_location,
                    target_dir,
                    size,
                    fix_layer,
                    solver_fmax,
                    solver_maxstep,
                    gpaw_calc):
    f = paropen(report_location,'a')
    parprint('Start clean slab calculation: ', file=f)
    if size != '1x1':
        clean_slab_gpw_path=target_dir+'/clean_slab/slab.gpw'
        
        if os.path.isfile(clean_slab_gpw_path):
            opt_slab, pre_calc = restart(clean_slab_gpw_path)
            pre_kpts=list(pre_calc.__dict__['parameters']['kpts'])
            set_kpts=list(gpaw_calc.__dict__['parameters']['kpts'])
            if pre_kpts == set_kpts:
                parprint('\t'+size+' clean slab is pre-calculated with kpts matched.',file=f)
            else:
                parprint('\t'+size+' clean slab pre-calculated has different kpts. Clean slab needs to re-calculate.', file=f)
                parprint('\t'+'Calculating '+size+' clean slab...',file=f)
                clean_slab=read(target_dir+'/clean_slab/input.traj')
                opt_slab=clean_slab_calculator(clean_slab,fix_layer,gpaw_calc,target_dir,solver_fmax,solver_maxstep)
        else:
            parprint('\t'+size+' clean slab is not pre-calculated.',file=f)
            parprint('\t'+'Calculating '+size+' clean slab...',file=f)
            interm_gpw=target_dir+'/clean_slab/slab_interm.gpw'
            if os.path.isfile(interm_gpw):
                clean_slab, gpaw_calc=restart(interm_gpw)
            else:
                clean_slab=read(target_dir+'/clean_slab/input.traj')
            opt_slab=clean_slab_calculator(clean_slab,fix_layer,gpaw_calc,target_dir,solver_fmax,solver_maxstep)
    else:
        parprint('\tslab size is 1x1. Clean slab calculation is skipped.', file=f)
        opt_slab=connect('final_database'+'/'+'surf.db').get_atoms(simple_name=element+'_'+miller_index)  
    parprint(' ',file=f)
    f.close()
    return opt_slab.get_potential_energy(), opt_slab.get_magnetic_moments()

def clean_slab_calculator(clean_slab,
                        fix_layer,
                        gpaw_calc,
                        target_dir,
                        solver_fmax,
                        solver_maxstep,
                        fix_option='bottom'):
    pbc_checker(clean_slab)
    calc_dict=gpaw_calc.__dict__['parameters']
    if calc_dict['spinpol']:
        clean_slab.set_initial_magnetic_moments([0]*len(clean_slab))
    slab_c_coord,cluster=detect_cluster(clean_slab)
    if fix_option == 'bottom':
        unique_cluster_index=sorted(set(cluster), key=cluster.index)[fix_layer-1]
        max_height_fix=max(slab_c_coord[cluster==unique_cluster_index])
        fix_mask=clean_slab.positions[:,2]<(max_height_fix+0.05) #add 0.05 Ang to make sure all bottom fixed
        fixed_atom_constrain=FixAtoms(mask=fix_mask)
        clean_slab.set_constraint(fixed_atom_constrain)
    clean_slab.set_calculator(gpaw_calc)
    opt.relax(clean_slab,target_dir+'/clean_slab',fmax=solver_fmax,maxstep=solver_maxstep)
    return clean_slab

def adsorption_energy_calculator(traj_file,
                                report_location,
                                opt_slab_energy,
                                adatom_pot_energy,
                                opt_slab_magmom,
                                gpaw_calc,
                                solver_fmax,
                                solver_maxstep,
                                magmom_option,
                                magmom_slab,
                                magmom_ads,
                                calc_type,
                                fix_layer,
                                fix_option = 'bottom'):
    
    interm_gpw='/'.join(traj_file.split('/')[:-1]+['slab_interm.gpw'])
    if os.path.isfile(interm_gpw):
        ads_slab, gpaw_calc=restart(interm_gpw)
    else:
        ads_slab=read(traj_file)
        pbc_checker(ads_slab)
        calc_dict=gpaw_calc.__dict__['parameters']
        if calc_dict['spinpol']:
            if magmom_option=='use_opt_slab':
                ads_slab=apply_magmom_opt_slab(opt_slab_magmom,ads_slab)
            elif magmom_option=='use_manual':
                ads_slab=apply_magmom_manual(ads_slab,magmom_slab,magmom_ads,adatom=1)
        fixed_line_constrain=FixedLine(a=-1,direction=[0,0,1])
        slab_c_coord,cluster=detect_cluster(ads_slab)
        if fix_option == 'bottom':
            unique_cluster_index=sorted(set(cluster), key=cluster.index)[fix_layer-1]
            max_height_fix=max(slab_c_coord[cluster==unique_cluster_index])
            fix_mask=ads_slab.positions[:,2]<(max_height_fix+0.05) #add 0.05 Ang to make sure all bottom fixed
        if calc_type == 'grid':
            fixed_atom_constrain=FixAtoms(mask=fix_mask)
            ads_slab.set_constraint([fixed_atom_constrain,fixed_line_constrain])
        elif calc_type == 'normal' and fix_option == 'bottom':
            fixed_atom_constrain=FixAtoms(mask=fix_mask)
            ads_slab.set_constraint(fixed_atom_constrain)
        ads_slab.set_calculator(gpaw_calc)
    location='/'.join(traj_file.split('/')[:-1])
    f=paropen(report_location,'a')
    parprint('Calculating '+('/'.join(location.split('/')[-2:]))+' adsorption site...',file=f)
    f.close()
    opt.relax(ads_slab,location,fmax=solver_fmax,maxstep=solver_maxstep)
    init_ads_site=traj_file.split('/')[-2]
    E_slab_ads=ads_slab.get_potential_energy()
    opt_slab_energy=opt_slab_energy
    adsorption_energy=E_slab_ads-(opt_slab_energy+adatom_pot_energy)
    final_ads_site=list(np.round(ads_slab.get_positions()[-1][:2],decimals=3))
    final_ads_site_str='_'.join([str(i) for i in final_ads_site])
    return init_ads_site, adsorption_energy, final_ads_site_str

def skip_ads_calculated(report_location,
                        all_gpw_files,
                        init_adsorbates_site_lst,
                        adsorption_energy_lst,
                        final_adsorbates_site_lst,
                        opt_slab_energy,
                        adatom_pot_energy):
    f = paropen(report_location,'a')
    parprint('Restarting...',file=f)
    for gpw_file in all_gpw_files:
        location='/'.join(gpw_file.split('/')[:-1])
        parprint('Skipping '+('/'.join(location.split('/')[-2:]))+' adsorption site...',file=f)
        atoms=restart(gpw_file)[0]
        init_adsorbates_site_lst.append(gpw_file.split('/')[-2])
        E_slab_ads=atoms.get_potential_energy()
        adsorption_energy=E_slab_ads-(opt_slab_energy+adatom_pot_energy)
        adsorption_energy_lst.append(adsorption_energy)
        final_ads_site=list(np.round(atoms.get_positions()[-1][:2],decimals=3))
        final_ads_site_str='_'.join([str(i) for i in final_ads_site])
        final_adsorbates_site_lst.append(final_ads_site_str)
    parprint(' ',file=f)
    f.close()
    return init_adsorbates_site_lst,adsorption_energy_lst,final_adsorbates_site_lst

def initialize_report(report_location,gpaw_calc):
    calc_dict=gpaw_calc.__dict__['parameters']
    if world.rank==0 and os.path.isfile(report_location):
        os.remove(report_location)
    f = paropen(report_location,'a')
    parprint('Initial Parameters:', file=f)
    parprint('\t'+'xc: '+calc_dict['xc'],file=f)
    parprint('\t'+'h: '+str(calc_dict['h']),file=f)
    parprint('\t'+'kpts: '+str(calc_dict['kpts']),file=f)
    parprint('\t'+'sw: '+str(calc_dict['occupations']),file=f)
    parprint('\t'+'spin polarized: '+str(calc_dict['spinpol']),file=f)
    if calc_dict['spinpol']:
        parprint('\t'+'magmom: initialize magnetic moment from slab calculation.',file=f)
    parprint(' ',file=f)
    f.close()

class ads_auto_select:
    def __init__(self,
                element,
                miller_index_tight,
                gpaw_calc,
                ads,
                adatom_pot_energy,
                solver_fmax,
                solver_max_step,
                restart_calc,
                magmom_slab,
                magmom_ads,
                magmom_option,
                size=(1,1), #xy size
                fix_layer=2,
                fix_option='bottom'):
        #initalize variable
        size_xy=str(size[0])+'x'+str(size[1])
        target_dir='results/'+element+'/'+'ads/'+size_xy+'/'+miller_index_tight
        report_location=target_dir+'_'+str(ads)+'_autocat_results_report.txt' 
        all_ads_file_loc=target_dir+'/'+'adsorbates/'+str(ads)+'/'
        ## TO-DO: need to figure out how to calculate adsorption energy for larger system
        # self.gpaw_calc=gpaw_calc
        # self.calc_dict=self.gpaw_calc.__dict__['parameters']
        # self.ads=ads
        # self.all_ads_file_loc=self.target_dir+'/'+'adsorbates/'+str(self.ads)+'/'
        # self.adatom_pot_energy=adatom_pot_energy
        ##generate report
        initialize_report(report_location, gpaw_calc)

        ##compute clean slab energy
        opt_slab_energy, opt_slab_magmom=get_clean_slab(element, miller_index_tight,
                                                    report_location, target_dir,size_xy,
                                                    fix_layer,solver_fmax,solver_max_step,
                                                    gpaw_calc)
        #opt_slab=self.get_clean_slab()

        ##start adsorption calculation
        adsorption_energy_dict={}
        init_adsorbates_site_lst=[]
        final_adsorbates_site_lst=[]
        adsorption_energy_lst=[]
        all_bridge_traj_files=glob(all_ads_file_loc+'bridge/*/input.traj')
        all_ontop_traj_files=glob(all_ads_file_loc+'ontop/*/input.traj')
        all_hollow_traj_files=glob(all_ads_file_loc+'hollow/*/input.traj')
        all_traj_files=all_bridge_traj_files+all_ontop_traj_files+all_hollow_traj_files
        all_bridge_gpw_files=glob(all_ads_file_loc+'bridge/*/slab.gpw')
        all_ontop_gpw_files=glob(all_ads_file_loc+'ontop/*/slab.gpw')
        all_hollow_gpw_files=glob(all_ads_file_loc+'hollow/*/slab.gpw')
        all_gpw_files=all_bridge_gpw_files+all_ontop_gpw_files+all_hollow_gpw_files

        ## restart 
        if restart_calc==True and len(all_gpw_files)>=1:
            init_adsorbates_site_lst,adsorption_energy_lst,final_adsorbates_site_lst=skip_ads_calculated(report_location,
                                                                                                all_gpw_files,
                                                                                                init_adsorbates_site_lst,
                                                                                                adsorption_energy_lst,
                                                                                                final_adsorbates_site_lst,
                                                                                                opt_slab_energy,
                                                                                                adatom_pot_energy)
            all_gpw_files_ads_site=['/'.join(i.split('/')[:-1]) for i in all_gpw_files]
            all_traj_files=[i for i in all_traj_files if '/'.join(i.split('/')[:-1]) not in all_gpw_files_ads_site]

        for traj_file in all_traj_files:
            #init_adsobates_site, adsorption_energy, final_adsorbates_site=self.adsorption_energy_calculator(traj_file,opt_slab)
            output_lst=adsorption_energy_calculator(traj_file,report_location,
                                                    opt_slab_energy,adatom_pot_energy,
                                                    opt_slab_magmom,gpaw_calc,
                                                    solver_fmax,solver_max_step,
                                                    magmom_option,magmom_slab,magmom_ads,
                                                    calc_type='normal',
                                                    fix_layer=fix_layer,fix_option = fix_option,
                                                    )
            init_adsorbates_site_lst.append(output_lst[0])
            adsorption_energy_lst.append(output_lst[1])
            final_adsorbates_site_lst.append(output_lst[2])
        
        adsorption_energy_dict['init_sites[x_y](Ang)']=init_adsorbates_site_lst
        adsorption_energy_dict['final_sites[x_y](Ang)']=final_adsorbates_site_lst
        adsorption_energy_dict['adsorption_energy(eV)']=adsorption_energy_lst
        ads_df=pd.DataFrame(adsorption_energy_dict)
        # ads_df.set_index('init_adsorbates_sites[x_y](Ang)',inplace=True)
        ads_df.sort_values(by=['adsorption_energy(eV)'],inplace=True)
        pd.set_option("display.max_rows", None, "display.max_columns", None)
        f=paropen(report_location,'a')
        parprint(ads_df,file=f)
        parprint('',file=f)
        f.close()
        min_adsorbates_site=ads_df.iloc[[0]]['init_sites[x_y](Ang)'].to_list()[0]
        lowest_ads_energy_slab=read(glob(all_ads_file_loc+'*/'+min_adsorbates_site+'/slab.traj')[0])

        #finalize
        final_slab_simple_name=element+'_'+miller_index_tight
        ads_db=connect('final_database/ads'+'_'+str(ads)+'_'+str(size_xy)+'.db')
        id=ads_db.reserve(name=final_slab_simple_name)
        if id is None:
            id=ads_db.get(name=final_slab_simple_name).id
            ads_db.update(id=id,atoms=lowest_ads_energy_slab,name=final_slab_simple_name,
                        ads_pot_e=float(ads_df.iloc[[0]]['adsorption_energy(eV)'].to_list()[0]))
        else:
            ads_db.write(lowest_ads_energy_slab,
                        id=id,
                        name=final_slab_simple_name,
                        ads_pot_e=float(ads_df.iloc[[0]]['adsorption_energy(eV)'].to_list()[0]))
        
        f=paropen(report_location,'a')
        parprint('Adsorption energy calculation complete.',file=f)
        parprint('Selected ads site is: ',file=f)
        parprint(min_adsorbates_site,file=f)
        f.close()

    # def get_clean_slab(self):
    #     f = paropen(self.report_location,'a')
    #     parprint('Start clean slab calculation: ', file=f)
    #     if self.size != '1x1':
    #         clean_slab_gpw_path=self.target_dir+'/clean_slab/slab.gpw'
    #         clean_slab=read(self.target_dir+'/clean_slab/input.traj')
    #         if os.path.isfile(clean_slab_gpw_path):
    #             opt_slab, pre_calc = restart(clean_slab_gpw_path)
    #             pre_kpts=pre_calc.__dict__['parameters']['kpts']
    #             set_kpts=self.calc_dict['kpts']
    #             if pre_kpts == set_kpts:
    #                 parprint('\t'+self.size+' clean slab is pre-calculated with kpts matched.',file=f)
    #             else:
    #                 parprint('\t'+self.size+' clean slab pre-calculated has different kpts. Clean slab needs to re-calculate.', file=f)
    #                 parprint('\t'+'Calculating '+self.size+' clean slab...',file=f)
    #                 opt_slab=self.clean_slab_calculator(clean_slab)
    #         else:
    #             parprint('\t'+self.size+' clean slab is not pre-calculated.',file=f)
    #             parprint('\t'+'Calculating '+self.size+' clean slab...',file=f)
    #             opt_slab=self.clean_slab_calculator(clean_slab)
    #     else:
    #         parprint('slab size is 1x1. Clean slab calculation is skipped.', file=f)
    #         opt_slab=connect('final_database'+'/'+'surf.db').get_atoms(simple_name=self.element+'_'+self.miller_index_tight)  
    #     f.close()
    #     return opt_slab

    # def clean_slab_calculator(self,clean_slab):
    #     pbc_checker(clean_slab)
    #     if self.calc_dict['spinpol']:
    #         clean_slab.set_initial_magnetic_moments([0]*len(clean_slab))
    #     slab_c_coord,cluster=detect_cluster(clean_slab)
    #     if self.fix_option == 'bottom':
    #         unique_cluster_index=sorted(set(cluster), key=cluster.index)[self.fix_layer-1]
    #         max_height_fix=max(slab_c_coord[cluster==unique_cluster_index])
    #         fix_mask=clean_slab.positions[:,2]<(max_height_fix+0.05) #add 0.05 Ang to make sure all bottom fixed
    #     else:
    #         raise RuntimeError('Only bottom fix option available now.')
    #     fixed_atom_constrain=FixAtoms(mask=fix_mask)
    #     clean_slab.set_constraint(fixed_atom_constrain)
    #     clean_slab.set_calculator(self.gpaw_calc)
    #     opt.relax(clean_slab,self.target_dir+'/clean_slab',fmax=self.solver_fmax,maxstep=self.solver_max_step)
    #     return clean_slab

    # def adsorption_energy_calculator(self,traj_file,opt_slab):
    #     ads_slab=read(traj_file)
    #     pbc_checker(ads_slab)
    #     if self.calc_dict['spinpol']:
    #         ads_slab=apply_magmom_opt_slab(opt_slab,ads_slab)
    #     slab_c_coord,cluster=detect_cluster(ads_slab)
    #     if self.fix_option == 'bottom':
    #         unique_cluster_index=sorted(set(cluster), key=cluster.index)[self.fix_layer-1]
    #         max_height_fix=max(slab_c_coord[cluster==unique_cluster_index])
    #         fix_mask=ads_slab.positions[:,2]<(max_height_fix+0.05) #add 0.05 Ang to make sure all bottom fixed
    #     else:
    #         raise RuntimeError('Only bottom fix option available now.')
    #     fixed_atom_constrain=FixAtoms(mask=fix_mask)
    #     ads_slab.set_constraint(fixed_atom_constrain)
    #     ads_slab.set_calculator(self.gpaw_calc)
    #     location='/'.join(traj_file.split('/')[:-1])
    #     f=paropen(self.report_location,'a')
    #     parprint('Calculating '+('/'.join(location.split('/')[-2:]))+' adsorption site...',file=f)
    #     f.close()
    #     opt.relax(ads_slab,location,fmax=self.solver_fmax,maxstep=self.solver_max_step)
    #     init_ads_site=traj_file.split('/')[-2]
    #     E_slab_ads=ads_slab.get_potential_energy()
    #     opt_slab_energy=opt_slab.get_potential_energy()*int(self.size[0])*int(self.size[2])
    #     adsorption_energy=E_slab_ads-(opt_slab_energy+self.adatom_pot_energy)
    #     final_ads_site=list(np.round(ads_slab.get_positions()[-1][:2],decimals=3))
    #     final_ads_site_str='_'.join([str(i) for i in final_ads_site])
    #     return init_ads_site, adsorption_energy, final_ads_site_str

    # def apply_magmom_opt_slab(self,opt_slab,ads_slab):
    #     slab_formula=ads_slab.get_chemical_symbols()
    #     magmom=opt_slab.get_magnetic_moments()
    #     magmom_ls=np.append(magmom,np.mean(magmom))
    #     magmom_ls[slab_formula.index(self.ads)]=0
    #     ads_slab.set_initial_magnetic_moments(magmom_ls)

    # def initialize_report(self,report_location,gpaw_calc):
    #     calc_dict=gpaw_calc.__dict__['parameters']
    #     if world.rank==0 and os.path.isfile(report_location):
    #         os.remove(report_location)
    #     f = paropen(report_location,'a')
    #     parprint('Initial Parameters:', file=f)
    #     parprint('\t'+'xc: '+calc_dict['xc'],file=f)
    #     parprint('\t'+'h: '+str(calc_dict['h']),file=f)
    #     parprint('\t'+'kpts: '+str(calc_dict['kpts']),file=f)
    #     parprint('\t'+'sw: '+str(calc_dict['occupations']),file=f)
    #     parprint('\t'+'spin polarized: '+str(calc_dict['spinpol']),file=f)
    #     if calc_dict['spinpol']:
    #         parprint('\t'+'magmom: initialize magnetic moment from slab calculation.',file=f)
    #     parprint(' ',file=f)
    #     f.close()


class ads_grid_calc:
    def __init__(self,
                element,
                miller_index_tight,
                gpaw_calc,
                ads,
                adatom_pot_energy,
                solver_fmax,
                solver_max_step,
                restart_calc,
                size,
                fix_layer=2,
                fix_option='bottom'):
        #initalize  variables
        size_xy=str(size[0])+'x'+str(size[1])
        target_dir='results/'+element+'/'+'ads/'+size_xy+'/'+miller_index_tight
        report_location=target_dir+'_grid_results_report.txt' 
        all_ads_file_loc=target_dir+'/'+'adsorbates/'+str(ads)+'/'
        ## TO-DO: need to figure out how to calculate adsorption energy for larger system
        # self.gpaw_calc=gpaw_calc
        # self.calc_dict=self.gpaw_calc.__dict__['parameters']
        # self.ads=ads
        #self.all_ads_file_loc=self.target_dir+'/'+'adsorbates/'+str(self.ads)+'/'
        #self.adatom_pot_energy=adatom_pot_energy

        ##generate report
        initialize_report(report_location,gpaw_calc)

        
        ##compute clean slab energy
        opt_slab_energy, opt_slab_magmom=get_clean_slab(element, miller_index_tight,
                                            report_location, target_dir, size_xy,
                                            fix_layer,solver_fmax,solver_max_step,
                                            gpaw_calc)
        ##start adsorption calculation
        adsorption_energy_dict={}
        init_adsorbates_site_lst=[]
        adsorption_energy_lst=[]
        final_adsorbates_site_lst=[]
        all_traj_files=glob(all_ads_file_loc+'grid/*/input.traj')
        all_gpw_files=glob(all_ads_file_loc+'grid/*/slab.gpw')

        ## restart
        if restart_calc==True and len(all_gpw_files)>=1:
            init_adsorbates_site_lst,adsorption_energy_lst=skip_ads_calculated(report_location,
                                                                            all_gpw_files,
                                                                            init_adsorbates_site_lst,
                                                                            adsorption_energy_lst,
                                                                            final_adsorbates_site_lst,
                                                                            opt_slab_energy,
                                                                            adatom_pot_energy)[0:2]
            all_gpw_files_ads_site=['/'.join(i.split('/')[:-1]) for i in all_gpw_files]
            all_traj_files=[i for i in all_traj_files if '/'.join(i.split('/')[:-1]) not in all_gpw_files_ads_site]


        for traj_file in all_traj_files:
            output_lst=adsorption_energy_calculator(traj_file,report_location,
                                                    opt_slab_energy,adatom_pot_energy,
                                                    opt_slab_magmom,gpaw_calc,
                                                    solver_fmax,solver_max_step,
                                                    calc_type='grid',
                                                    fix_layer=fix_layer,fix_option = 'bottom',
                                                    )
            init_adsorbates_site_lst.append(output_lst[0])
            adsorption_energy_lst.append(output_lst[1])
        
        adsorption_energy_dict['init_sites[x_y](Ang)']=init_adsorbates_site_lst
        adsorption_energy_dict['adsorption_energy(eV)']=adsorption_energy_lst
        ads_df=pd.DataFrame(adsorption_energy_dict)
        #ads_df.set_index('init_adsorbates_sites[x_y](Ang)',inplace=True)
        ads_df.sort_values(by=['adsorption_energy(eV)'],inplace=True)
        ads_df.to_csv(target_dir+'_ads_grid.csv')
        pd.set_option("display.max_rows", None, "display.max_columns", None)
        f=paropen(report_location,'a')
        parprint(ads_df,file=f)
        parprint('',file=f)
        parprint('Grid adsorption energy calculation complete.',file=f)
        f.close()

    # def get_clean_slab(self):
    #     f = paropen(self.report_location,'a')
    #     parprint('Start clean slab calculation: ', file=f)
    #     if self.size != '1x1':
    #         clean_slab_gpw_path=self.target_dir+'/clean_slab/slab.gpw'
    #         clean_slab=read(self.target_dir+'/clean_slab/input.traj')
    #         if os.path.isfile(clean_slab_gpw_path):
    #             opt_slab, pre_calc = restart(clean_slab_gpw_path)
    #             pre_kpts=pre_calc.__dict__['parameters']['kpts']
    #             set_kpts=self.calc_dict['kpts']
    #             if pre_kpts == set_kpts:
    #                 parprint('\t'+self.size+' clean slab is pre-calculated with kpts matched.',file=f)
    #             else:
    #                 parprint('\t'+self.size+' clean slab pre-calculated has different kpts. Clean slab needs to re-calculate.', file=f)
    #                 parprint('\t'+'Calculating '+self.size+' clean slab...',file=f)
    #                 opt_slab=self.clean_slab_calculator(clean_slab)
    #         else:
    #             parprint('\t'+self.size+' clean slab is not pre-calculated.',file=f)
    #             parprint('\t'+'Calculating '+self.size+' clean slab...',file=f)
    #             opt_slab=self.clean_slab_calculator(clean_slab)
    #     else:
    #         parprint('slab size is 1x1. Clean slab calculation is skipped.', file=f)
    #         opt_slab=connect('final_database'+'/'+'surf.db').get_atoms(simple_name=self.element+'_'+self.miller_index_tight)  
    #     f.close()
    #     return opt_slab

    # def clean_slab_calculator(self,clean_slab):
    #     pbc_checker(clean_slab)
    #     if self.calc_dict['spinpol']:
    #         clean_slab.set_initial_magnetic_moments([0]*len(clean_slab))
    #     slab_c_coord,cluster=detect_cluster(clean_slab)
    #     if self.fix_option == 'bottom':
    #         unique_cluster_index=sorted(set(cluster), key=cluster.index)[self.fix_layer-1]
    #         max_height_fix=max(slab_c_coord[cluster==unique_cluster_index])
    #         fix_mask=clean_slab.positions[:,2]<(max_height_fix+0.05) #add 0.05 Ang to make sure all bottom fixed
    #     else:
    #         raise RuntimeError('Only bottom fix option available now.')
    #     fixed_atom_constrain=FixAtoms(mask=fix_mask)
    #     clean_slab.set_constraint(fixed_atom_constrain)
    #     clean_slab.set_calculator(self.gpaw_calc)
    #     opt.relax(clean_slab,self.target_dir+'/clean_slab',fmax=self.solver_fmax,maxstep=self.solver_max_step)
    #     return clean_slab

    # def adsorption_energy_calculator(self,traj_file,opt_slab):
    #     ads_slab=read(traj_file)
    #     pbc_checker(ads_slab)
    #     if self.calc_dict['spinpol']:
    #         ads_slab=apply_magmom_opt_slab(opt_slab,ads_slab)
    #     fixed_line_constrain=FixedLine(a=-1,direction=[0,0,1])
    #     slab_c_coord,cluster=detect_cluster(ads_slab)
    #     if self.fix_option == 'bottom':
    #         unique_cluster_index=sorted(set(cluster), key=cluster.index)[self.fix_layer-1]
    #         max_height_fix=max(slab_c_coord[cluster==unique_cluster_index])
    #         fix_mask=ads_slab.positions[:,2]<(max_height_fix+0.05) #add 0.05 Ang to make sure all bottom fixed
    #     else:
    #         raise RuntimeError('Only bottom fix option available now.')
    #     fixed_atom_constrain=FixAtoms(mask=fix_mask)
    #     ads_slab.set_constraint([fixed_atom_constrain,fixed_line_constrain])
    #     ads_slab.set_calculator(self.gpaw_calc)
    #     location='/'.join(traj_file.split('/')[:-1])
    #     f=paropen(self.report_location,'a')
    #     parprint('Calculating '+('/'.join(location.split('/')[-2:]))+' adsorption site...',file=f)
    #     f.close()
    #     opt.relax(ads_slab,location,fmax=self.solver_fmax,maxstep=self.solver_max_step)
    #     init_ads_site=traj_file.split('/')[-2]
    #     adsorption_energy=ads_slab.get_potential_energy()-(opt_slab.get_potential_energy()+self.adatom_pot_energy)
    #     return init_ads_site, adsorption_energy

    # def apply_magmom_opt_slab(self,opt_slab,ads_slab):
    #     slab_formula=ads_slab.get_chemical_symbols()
    #     magmom=opt_slab.get_magnetic_moments()
    #     magmom_ls=np.append(magmom,np.mean(magmom))
    #     magmom_ls[slab_formula.index(self.ads)]=0
    #     ads_slab.set_initial_magnetic_moments(magmom_ls)

    # def initialize_report(self):
    #     if world.rank==0 and os.path.isfile(self.report_location):
    #         os.remove(self.report_location)
    #     f = paropen(self.report_location,'a')
    #     parprint('Initial Parameters:', file=f)
    #     parprint('\t'+'xc: '+self.calc_dict['xc'],file=f)
    #     parprint('\t'+'h: '+str(self.calc_dict['h']),file=f)
    #     parprint('\t'+'kpts: '+str(self.calc_dict['kpts']),file=f)
    #     parprint('\t'+'sw: '+str(self.calc_dict['occupations']),file=f)
    #     parprint('\t'+'spin polarized: '+str(self.calc_dict['spinpol']),file=f)
    #     if self.calc_dict['spinpol']:
    #         parprint('\t'+'magmom: initial magnetic moment from slab calculation.',file=f)
    #     parprint(' ',file=f)
    #     f.close()

class ads_lowest_ads_site_calc:
    def __init__(self,
                element,
                miller_index_tight,
                gpaw_calc,
                ads,
                adatom_pot_energy,
                solver_fmax,
                solver_max_step,
                restart_calc,
                magmom_slab,
                magmom_ads,
                magmom_option,
                size, #xy size
                fix_layer=2,
                fix_option='bottom'):
        #initalize
        ##globlalize variable
        size_xy=str(size[0])+'x'+str(size[1])
        target_dir='results/'+element+'/'+'ads/'+size_xy+'/'+miller_index_tight
        report_location=target_dir+'_lowest_ads_results_report.txt' 
        all_ads_file_loc=target_dir+'/'+'adsorbates/'+str(ads)+'/'
        ##generate report
        initialize_report(report_location, gpaw_calc)

        ##compute clean slab energy
        opt_slab_energy, opt_slab_magmom=get_clean_slab(element, miller_index_tight,
                                    report_location, target_dir, size_xy,
                                    fix_layer,solver_fmax,solver_max_step,
                                    gpaw_calc)

        ##start adsorption calculation
        adsorption_energy_dict={}
        init_adsorbates_site_lst=[]
        final_adsorbates_site_lst=[]
        adsorption_energy_lst=[]
        all_traj_files=glob(all_ads_file_loc+'lowest_ads_site/*/input.traj')
        all_gpw_files=glob(all_ads_file_loc+'lowest_ads_site/*/slab.gpw')


        if restart_calc==True and len(all_gpw_files)>=1:
            init_adsorbates_site_lst,adsorption_energy_lst=skip_ads_calculated(report_location,
                                                                            all_gpw_files,
                                                                            init_adsorbates_site_lst,
                                                                            adsorption_energy_lst,
                                                                            final_adsorbates_site_lst,
                                                                            opt_slab_energy,
                                                                            adatom_pot_energy)[0:2]
            all_gpw_files_ads_site=['/'.join(i.split('/')[:-1]) for i in all_gpw_files]
            all_traj_files=[i for i in all_traj_files if '/'.join(i.split('/')[:-1]) not in all_gpw_files_ads_site]

        for traj_file in all_traj_files:
            output_lst=adsorption_energy_calculator(traj_file,report_location,
                                                    opt_slab_energy,adatom_pot_energy,
                                                    opt_slab_magmom,gpaw_calc,
                                                    solver_fmax,solver_max_step,
                                                    magmom_option,magmom_slab,magmom_ads,
                                                    calc_type='normal',
                                                    fix_layer=fix_layer,fix_option = fix_option,
                                                    )
            init_adsorbates_site_lst.append(output_lst[0])
            adsorption_energy_lst.append(output_lst[1])
            final_adsorbates_site_lst.append(output_lst[2])
        
        adsorption_energy_dict['init_sites[x_y](Ang)']=init_adsorbates_site_lst
        adsorption_energy_dict['final_sites[x_y](Ang)']=final_adsorbates_site_lst
        adsorption_energy_dict['adsorption_energy(eV)']=adsorption_energy_lst
        ads_df=pd.DataFrame(adsorption_energy_dict)
        # ads_df.set_index('init_adsorbates_sites[x_y](Ang)',inplace=True)
        ads_df.sort_values(by=['adsorption_energy(eV)'],inplace=True)
        pd.set_option("display.max_rows", None, "display.max_columns", None)
        f=paropen(report_location,'a')
        parprint(ads_df,file=f)
        parprint('',file=f)
        f.close()
        min_adsorbates_site=ads_df.iloc[[0]]['init_sites[x_y](Ang)'].to_list()[0]
        lowest_ads_energy_slab=read(glob(all_ads_file_loc+'*/'+min_adsorbates_site+'/slab.traj')[0])

        #finalize
        final_slab_simple_name=element+'_'+miller_index_tight
        ads_db=connect('final_database/ads_'+str(ads)+'_'+size_xy+'.db')
        id=ads_db.reserve(name=final_slab_simple_name)
        if id is None:
            id=ads_db.get(name=final_slab_simple_name).id
            ads_db.update(id=id,atoms=lowest_ads_energy_slab,name=final_slab_simple_name,
                        ads_pot_e=float(ads_df.iloc[[0]]['adsorption_energy(eV)'].to_list()[0]))
        else:
            ads_db.write(lowest_ads_energy_slab,
                        id=id,
                        name=final_slab_simple_name,
                        ads_pot_e=float(ads_df.iloc[[0]]['adsorption_energy(eV)'].to_list()[0]))
        
        f=paropen(report_location,'a')
        parprint('Adsorption energy calculation complete.',file=f)
        parprint('Selected ads site is: ',file=f)
        parprint(min_adsorbates_site,file=f)
        f.close()

    # def get_clean_slab(self):
    #     f = paropen(self.report_location,'a')
    #     parprint('Start clean slab calculation: ', file=f)
    #     if self.size != '1x1':
    #         clean_slab_gpw_path=self.target_dir+'/clean_slab/slab.gpw'
    #         clean_slab=read(self.target_dir+'/clean_slab/input.traj')
    #         if os.path.isfile(clean_slab_gpw_path):
    #             opt_slab, pre_calc = restart(clean_slab_gpw_path)
    #             pre_kpts=pre_calc.__dict__['parameters']['kpts']
    #             set_kpts=self.calc_dict['kpts']
    #             if pre_kpts == set_kpts:
    #                 parprint('\t'+self.size+' clean slab is pre-calculated with kpts matched.',file=f)
    #             else:
    #                 parprint('\t'+self.size+' clean slab pre-calculated has different kpts. Clean slab needs to re-calculate.', file=f)
    #                 parprint('\t'+'Calculating '+self.size+' clean slab...',file=f)
    #                 opt_slab=self.clean_slab_calculator(clean_slab)
    #         else:
    #             parprint('\t'+self.size+' clean slab is not pre-calculated.',file=f)
    #             parprint('\t'+'Calculating '+self.size+' clean slab...',file=f)
    #             opt_slab=self.clean_slab_calculator(clean_slab)
    #     else:
    #         parprint('slab size is 1x1. Clean slab calculation is skipped.', file=f)
    #         opt_slab=connect('final_database'+'/'+'surf.db').get_atoms(simple_name=self.element+'_'+self.miller_index_tight)  
    #     parprint(' ',file=f)
    #     f.close()
    #     return opt_slab

    # def clean_slab_calculator(self,clean_slab):
    #     pbc_checker(clean_slab)
    #     if self.calc_dict['spinpol']:
    #         clean_slab.set_initial_magnetic_moments([0]*len(clean_slab))
    #     slab_c_coord,cluster=detect_cluster(clean_slab)
    #     if self.fix_option == 'bottom':
    #         unique_cluster_index=sorted(set(cluster), key=cluster.index)[self.fix_layer-1]
    #         max_height_fix=max(slab_c_coord[cluster==unique_cluster_index])
    #         fix_mask=clean_slab.positions[:,2]<(max_height_fix+0.05) #add 0.05 Ang to make sure all bottom fixed
    #     else:
    #         raise RuntimeError('Only bottom fix option available now.')
    #     fixed_atom_constrain=FixAtoms(mask=fix_mask)
    #     clean_slab.set_constraint(fixed_atom_constrain)
    #     clean_slab.set_calculator(self.gpaw_calc)
    #     opt.relax(clean_slab,self.target_dir+'/clean_slab',fmax=self.solver_fmax,maxstep=self.solver_max_step)
    #     return clean_slab

    # def adsorption_energy_calculator(self,traj_file,opt_slab):
    #     ads_slab=read(traj_file)
    #     pbc_checker(ads_slab)
    #     if self.calc_dict['spinpol']:
    #         ads_slab=apply_magmom_opt_slab(opt_slab,ads_slab)
    #     slab_c_coord,cluster=detect_cluster(ads_slab)
    #     if self.fix_option == 'bottom':
    #         unique_cluster_index=sorted(set(cluster), key=cluster.index)[self.fix_layer-1]
    #         max_height_fix=max(slab_c_coord[cluster==unique_cluster_index])
    #         fix_mask=ads_slab.positions[:,2]<(max_height_fix+0.05) #add 0.05 Ang to make sure all bottom fixed
    #     else:
    #         raise RuntimeError('Only bottom fix option available now.')
    #     fixed_atom_constrain=FixAtoms(mask=fix_mask)
    #     ads_slab.set_constraint(fixed_atom_constrain)
    #     ads_slab.set_calculator(self.gpaw_calc)
    #     location='/'.join(traj_file.split('/')[:-1])
    #     f=paropen(self.report_location,'a')
    #     parprint('\tCalculating '+('/'.join(location.split('/')[-2:]))+' adsorption site...',file=f)
    #     f.close()
    #     opt.relax(ads_slab,location,fmax=self.solver_fmax,maxstep=self.solver_max_step)
    #     init_ads_site=traj_file.split('/')[-2]
    #     E_slab_ads=ads_slab.get_potential_energy()
    #     opt_slab_energy=opt_slab.get_potential_energy()
    #     adsorption_energy=E_slab_ads-(opt_slab_energy+self.adatom_pot_energy)
    #     final_ads_site=list(np.round(ads_slab.get_positions()[-1][:2],decimals=3))
    #     final_ads_site_str='_'.join([str(i) for i in final_ads_site])
    #     return init_ads_site, adsorption_energy, final_ads_site_str

    # def initialize_report(self):
    #     if world.rank==0 and os.path.isfile(self.report_location):
    #         os.remove(self.report_location)
    #     f = paropen(self.report_location,'a')
    #     parprint('Initial Parameters:', file=f)
    #     parprint('\t'+'xc: '+self.calc_dict['xc'],file=f)
    #     parprint('\t'+'h: '+str(self.calc_dict['h']),file=f)
    #     parprint('\t'+'kpts: '+str(self.calc_dict['kpts']),file=f)
    #     parprint('\t'+'sw: '+str(self.calc_dict['occupations']),file=f)
    #     parprint('\t'+'spin polarized: '+str(self.calc_dict['spinpol']),file=f)
    #     if self.calc_dict['spinpol']:
    #         parprint('\t'+'magmom: initial magnetic moment from slab calculation.',file=f)
    #     parprint(' ',file=f)
    #     f.close()


class ads_NN_interact_calc:
    def __init__(self,
                element,
                miller_index_tight,
                gpaw_calc,
                ads,
                solver_fmax,
                solver_max_step,
                restart_calc,
                size, #xy size
                sub_dir,
                fix_layer=2,
                fix_option='bottom'):
        #initalize
        ##globlalize variable
        size_xy=str(size[0])+'x'+str(size[1])
        target_dir='results/'+element+'/'+'ads/'+size_xy+'/'+miller_index_tight
        #report_location=target_dir+'_lowest_ads_results_report.txt' 
        all_ads_file_loc=target_dir+'/'+'adsorbates/'+str(ads)+'/'


        ##start adsorption calculation
        # adsorption_energy_dict={}
        # init_adsorbates_site_lst=[]
        # final_adsorbates_site_lst=[]
        # adsorption_energy_lst=[]
        all_traj_files=glob(all_ads_file_loc+sub_dir+'/*/input.traj')
        all_gpw_files=glob(all_ads_file_loc+sub_dir+'/*/slab.gpw')


        if restart_calc==True and len(all_gpw_files)>=1:
            all_gpw_files_ads_site=['/'.join(i.split('/')[:-1]) for i in all_gpw_files]
            all_traj_files=[i for i in all_traj_files if '/'.join(i.split('/')[:-1]) not in all_gpw_files_ads_site]

        for traj_file in all_traj_files:
            interm_gpw='/'.join(traj_file.split('/')[:-1]+['slab_interm.gpw'])
            if os.path.isfile(interm_gpw):
                ads_slab, gpaw_calc=restart(interm_gpw)
            else:
                ads_slab=read(traj_file)
                pbc_checker(ads_slab)
                calc_dict=gpaw_calc.__dict__['parameters']
                if calc_dict['spinpol']:
                    raise RuntimeError('spin polarization calculation not supported.')
                slab_c_coord,cluster=detect_cluster(ads_slab)
                if fix_option == 'bottom':
                    unique_cluster_index=sorted(set(cluster), key=cluster.index)[fix_layer-1]
                    max_height_fix=max(slab_c_coord[cluster==unique_cluster_index])
                    fix_mask=ads_slab.positions[:,2]<(max_height_fix+0.05) #add 0.05 Ang to make sure all bottom fixed
                else:
                    raise RuntimeError('Only bottom fix option available now.')
                fixed_atom_constrain=FixAtoms(mask=fix_mask)
                ads_slab.set_constraint(fixed_atom_constrain)
                ads_slab.set_calculator(gpaw_calc)
            location='/'.join(traj_file.split('/')[:-1])
            opt.relax(ads_slab,location,fmax=solver_fmax,maxstep=solver_max_step)


class ads_custom_ads_site_calc:
    def __init__(self,
                element,
                miller_index_tight,
                gpaw_calc,
                ads,
                adatom_pot_energy,
                solver_fmax,
                solver_max_step,
                restart_calc,
                size, #xy size
                fix_layer=2,
                fix_option='bottom'):
        #initalize
        ##globlalize variable
        size_xy=str(size[0])+'x'+str(size[1])
        target_dir='results/'+element+'/'+'ads/'+size_xy+'/'+miller_index_tight
        report_location=target_dir+'_custom_ads_results_report.txt' 
        all_ads_file_loc=target_dir+'/'+'adsorbates/'+str(ads)+'/'
        ##generate report
        initialize_report(report_location, gpaw_calc)

        ##compute clean slab energy
        opt_slab_energy, opt_slab_magmom=get_clean_slab(element, miller_index_tight,
                                    report_location, target_dir, size_xy,
                                    fix_layer,solver_fmax,solver_max_step,
                                    gpaw_calc)

        ##start adsorption calculation
        adsorption_energy_dict={}
        init_adsorbates_site_lst=[]
        final_adsorbates_site_lst=[]
        adsorption_energy_lst=[]
        all_traj_files=glob(all_ads_file_loc+'custom/*/input.traj')
        all_gpw_files=glob(all_ads_file_loc+'custom/*/slab.gpw')


        if restart_calc==True and len(all_gpw_files)>=1:
            init_adsorbates_site_lst,adsorption_energy_lst=skip_ads_calculated(report_location,
                                                                            all_gpw_files,
                                                                            init_adsorbates_site_lst,
                                                                            adsorption_energy_lst,
                                                                            final_adsorbates_site_lst,
                                                                            opt_slab_energy,
                                                                            adatom_pot_energy)[0:2]
            all_gpw_files_ads_site=['/'.join(i.split('/')[:-1]) for i in all_gpw_files]
            all_traj_files=[i for i in all_traj_files if '/'.join(i.split('/')[:-1]) not in all_gpw_files_ads_site]

        for traj_file in all_traj_files:
            output_lst=adsorption_energy_calculator(traj_file,report_location,
                                                    opt_slab_energy,adatom_pot_energy,
                                                    opt_slab_magmom,gpaw_calc,
                                                    solver_fmax,solver_max_step,
                                                    calc_type='normal',
                                                    fix_layer=fix_layer,fix_option = 'bottom',
                                                    )
            init_adsorbates_site_lst.append(output_lst[0])
            adsorption_energy_lst.append(output_lst[1])
            final_adsorbates_site_lst.append(output_lst[2])
        
        adsorption_energy_dict['init_sites[x_y](Ang)']=init_adsorbates_site_lst
        adsorption_energy_dict['final_sites[x_y](Ang)']=final_adsorbates_site_lst
        adsorption_energy_dict['adsorption_energy(eV)']=adsorption_energy_lst
        ads_df=pd.DataFrame(adsorption_energy_dict)
        # ads_df.set_index('init_adsorbates_sites[x_y](Ang)',inplace=True)
        ads_df.sort_values(by=['adsorption_energy(eV)'],inplace=True)
        pd.set_option("display.max_rows", None, "display.max_columns", None)
        f=paropen(report_location,'a')
        parprint(ads_df,file=f)
        parprint('',file=f)
        f.close()
        min_adsorbates_site=ads_df.iloc[[0]]['init_sites[x_y](Ang)'].to_list()[0]
        #lowest_ads_energy_slab=read(glob(all_ads_file_loc+'*/'+min_adsorbates_site+'/slab.traj')[0])

        #finalize
        # final_slab_simple_name=element+'_'+miller_index_tight
        # ads_db=connect('final_database/ads_'+size_xy+'.db')
        # id=ads_db.reserve(name=final_slab_simple_name)
        # if id is None:
        #     id=ads_db.get(name=final_slab_simple_name).id
        #     ads_db.update(id=id,atoms=lowest_ads_energy_slab,name=final_slab_simple_name,
        #                 ads_pot_e=float(ads_df.iloc[[0]]['adsorption_energy(eV)'].to_list()[0]))
        # else:
        #     ads_db.write(lowest_ads_energy_slab,
        #                 id=id,
        #                 name=final_slab_simple_name,
        #                 ads_pot_e=float(ads_df.iloc[[0]]['adsorption_energy(eV)'].to_list()[0]))
        
        f=paropen(report_location,'a')
        parprint('Adsorption energy calculation complete.',file=f)
        parprint('Selected ads site is: ',file=f)
        parprint(min_adsorbates_site,file=f)
        f.close()
