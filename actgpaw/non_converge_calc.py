import os
from typing import Type
from ase.parallel import paropen, parprint, world
from ase.db import connect
from ase.io import read
from glob import glob
import numpy as np
from gpaw import restart
import actgpaw.optimizer as opt
import sys
from ase.calculators.calculator import kptdensity2monkhorstpack as kdens2mp
import itertools
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import fcluster, linkage
from ase.constraints import FixAtoms,FixedLine
import pandas as pd

def pbc_checker(slab):
    anlges_arg=[angle != 90.0000 for angle in np.round(slab.cell.angles(),decimals=4)[:2]]
    if np.any(anlges_arg):
        slab.pbc=[1,1,1]
    else:
        slab.pbc=[1,1,0]

def detect_cluster(slab,tol=0.1):
    n=len(slab)
    dist_matrix=np.zeros((n, n))
    slab_c=np.sort(slab.get_positions()[:,2])
    for i, j in itertools.combinations(list(range(n)), 2):
        if i != j:
            cdist = np.abs(slab_c[i] - slab_c[j])
            dist_matrix[i, j] = cdist
            dist_matrix[j, i] = cdist
    condensed_m = squareform(dist_matrix)
    z = linkage(condensed_m)
    clusters = fcluster(z, tol, criterion="distance")
    return slab_c,list(clusters)

class ads_auto_select:
    def __init__(self,
                element,
                miller_index,
                gpaw_calc,
                ads,
                ads_pot_energy,
                solver_fmax,
                solver_max_step,
                size,
                restart_calc,
                fix_layer=2,
                fix_option='bottom'):
        #initalize
        ##globlalize variable
        self.element=element
        self.solver_max_step=solver_max_step
        self.solver_fmax=solver_fmax
        self.fix_option=fix_option
        self.fix_layer=fix_layer
        self.miller_index_tight=miller_index
        self.miller_index_loose=tuple(map(int,miller_index)) #tuple
        self.gpaw_calc=gpaw_calc
        self.target_dir='results/'+element+'/'+'ads/'
        self.report_location=self.target_dir+self.miller_index_tight+'_autocat_results_report.txt' 
        ## TO-DO: need to figure out how to calculate adsorption energy for larger system
        self.gpaw_calc=gpaw_calc
        self.calc_dict=self.gpaw_calc.__dict__['parameters']
        self.ads=ads
        self.all_ads_file_loc=self.target_dir+self.miller_index_tight+'/'+'adsorbates/'+str(self.ads)+'/'
        self.ads_pot_energy=ads_pot_energy
        ##generate report
        self.initialize_report()
        ## connect to opt_slab
        opt_slab=connect('final_database'+'/'+'surf.db').get_atoms(simple_name=self.element+'_'+self.miller_index_tight)
        
        ##start adsorption calculation
        adsorption_energy_dict={}
        init_adsorbates_site_lst=[]
        final_adsorbates_site_lst=[]
        adsorption_energy_lst=[]
        all_bridge_traj_files=glob(self.all_ads_file_loc+'bridge/*/input.traj')
        all_ontop_traj_files=glob(self.all_ads_file_loc+'ontop/*/input.traj')
        all_hollow_traj_files=glob(self.all_ads_file_loc+'hollow/*/input.traj')
        all_traj_files=all_bridge_traj_files+all_ontop_traj_files+all_hollow_traj_files
        all_bridge_gpw_files=glob(self.all_ads_file_loc+'bridge/*/slab.gpw')
        all_ontop_gpw_files=glob(self.all_ads_file_loc+'ontop/*/slab.gpw')
        all_hollow_gpw_files=glob(self.all_ads_file_loc+'hollow/*/slab.gpw')
        all_gpw_files=all_bridge_gpw_files+all_ontop_gpw_files+all_hollow_gpw_files

        ## restart 
        if restart_calc==True and len(all_gpw_files)>=1:
            f = paropen(self.report_location,'a')
            parprint('Restarting...',file=f)
            for gpw_file in all_gpw_files:
                location='/'.join(gpw_file.split('/')[:-1])
                parprint('Skipping '+('/'.join(location.split('/')[-2:]))+' adsorption site...',file=f)
                atoms=restart(gpw_file)[0]
                init_adsorbates_site_lst.append(gpw_file.split('/')[-2])
                adsorption_energy=atoms.get_potential_energy()-(opt_slab.get_potential_energy()+self.ads_pot_energy)
                adsorption_energy_lst.append(adsorption_energy)
                final_ads_site=list(np.round(atoms.get_positions()[-1][:2],decimals=3))
                final_ads_site_str='_'.join([str(i) for i in final_ads_site])
                final_adsorbates_site_lst.append(final_ads_site_str)
            parprint(' ',file=f)
            f.close()
            all_gpw_files_ads_site=['/'.join(i.split('/')[:-1]) for i in all_gpw_files]
            all_traj_files=[i for i in all_traj_files if '/'.join(i.split('/')[:-1]) not in all_gpw_files_ads_site]
        
        for traj_file in all_traj_files:
            init_adsobates_site, adsorption_energy, final_adsorbates_site=self.adsorption_energy_calculator(traj_file,opt_slab)
            init_adsorbates_site_lst.append(init_adsobates_site)
            adsorption_energy_lst.append(adsorption_energy)
            final_adsorbates_site_lst.append(final_adsorbates_site)
        
        adsorption_energy_dict['init_sites[x_y](Ang)']=init_adsorbates_site_lst
        adsorption_energy_dict['final_sites[x_y](Ang)']=final_adsorbates_site_lst
        adsorption_energy_dict['adsorption_energy(eV)']=adsorption_energy_lst
        ads_df=pd.DataFrame(adsorption_energy_dict)
        # ads_df.set_index('init_adsorbates_sites[x_y](Ang)',inplace=True)
        ads_df.sort_values(by=['adsorption_energy(eV)'],inplace=True)
        f=paropen(self.report_location,'a')
        parprint(' ',file=f)
        for i in len(ads_df):
            parprint(ads_df.iloc[i],file=f)
        f.close()
        parprint(ads_df.iloc[[0]]['init_sites[x_y](Ang)'])
        min_adsorbates_site=ads_df.iloc[[0]]['init_sites[x_y](Ang)'].values()
        parprint(min_adsorbates_site)
        lowest_ads_energy_slab=read(self.all_ads_file_loc+'*/'+min_adsorbates_site+'/slab.traj')
        
        #finalize
        final_slab_simple_name=self.element+'_'+self.miller_index_tight
        ads_db=connect('final_database/ads_'+str(size)+'.db')
        id=ads_db.reserve(name=final_slab_simple_name)
        if id is None:
            id=ads_db.get(name=final_slab_simple_name).id
            ads_db.update(id=id,atoms=lowest_ads_energy_slab,name=final_slab_simple_name,
                        ads_pot_e=float(ads_df[min_adsorbates_site].values()))
        else:
            ads_db.write(lowest_ads_energy_slab,id=id,name=final_slab_simple_name,ads_pot_e=float(ads_df[min_adsorbates_site].values()))
        
        f=paropen(self.report_location,'a')
        parprint('Adsorption energy calculation complete.',file=f)
        parprint('Selected ads site is: ',file=f)
        parprint(ads_df.loc[[min_adsorbates_site]],file=f)
        f.close()

    def adsorption_energy_calculator(self,traj_file,opt_slab):
        ads_slab=read(traj_file)
        pbc_checker(ads_slab)
        if self.calc_dict['spinpol']:
            self.apply_magmom(opt_slab,ads_slab)
        slab_c_coord,cluster=detect_cluster(ads_slab)
        if self.fix_option == 'bottom':
            unique_cluster_index=sorted(set(cluster), key=cluster.index)[self.fix_layer-1]
            max_height_fix=max(slab_c_coord[cluster==unique_cluster_index])
            fix_mask=ads_slab.positions[:,2]<(max_height_fix+0.05) #add 0.05 Ang to make sure all bottom fixed
        else:
            raise RuntimeError('Only bottom fix option available now.')
        fixed_atom_constrain=FixAtoms(mask=fix_mask)
        ads_slab.set_constraint(fixed_atom_constrain)
        ads_slab.set_calculator(self.gpaw_calc)
        location='/'.join(traj_file.split('/')[:-1])
        f=paropen(self.report_location,'a')
        parprint('Calculating '+('/'.join(location.split('/')[-2:]))+' adsorption site...',file=f)
        f.close()
        opt.relax(ads_slab,location,fmax=self.solver_fmax,maxstep=self.solver_max_step)
        init_ads_site=traj_file.split('/')[-2]
        adsorption_energy=ads_slab.get_potential_energy()-(opt_slab.get_potential_energy()+self.ads_pot_energy)
        final_ads_site=list(np.round(ads_slab.get_positions()[-1][:2],decimals=3))
        final_ads_site_str='_'.join([str(i) for i in final_ads_site])
        return init_ads_site, adsorption_energy, final_ads_site_str

    def apply_magmom(self,opt_slab,ads_slab):
        slab_formula=ads_slab.get_chemical_symbols()
        magmom=opt_slab.get_magnetic_moments()
        magmom_ls=np.append(magmom,np.mean(magmom))
        magmom_ls[slab_formula.index(self.ads)]=0
        ads_slab.set_initial_magnetic_moments(magmom_ls)

    def initialize_report(self):
        if world.rank==0 and os.path.isfile(self.report_location):
            os.remove(self.report_location)
        f = paropen(self.report_location,'a')
        parprint('Initial Parameters:', file=f)
        parprint('\t'+'xc: '+self.calc_dict['xc'],file=f)
        parprint('\t'+'h: '+str(self.calc_dict['h']),file=f)
        parprint('\t'+'kpts: '+str(self.calc_dict['kpts']),file=f)
        parprint('\t'+'sw: '+str(self.calc_dict['occupations']),file=f)
        parprint('\t'+'spin polarized: '+str(self.calc_dict['spinpol']),file=f)
        if self.calc_dict['spinpol']:
            parprint('\t'+'magmom: initial magnetic moment from slab calculation.',file=f)
        parprint(' ',file=f)
        f.close()


class ads_grid_calc:
    def __init__(self,
                element,
                miller_index,
                gpaw_calc,
                ads,
                ads_pot_energy,
                solver_fmax,
                solver_max_step,
                size,
                restart_calc,
                fix_layer=2,
                fix_option='bottom'):
        #initalize
        ##globlalize variable
        self.element=element
        self.solver_max_step=solver_max_step
        self.solver_fmax=solver_fmax
        self.fix_option=fix_option
        self.fix_layer=fix_layer
        self.miller_index_tight=miller_index
        self.miller_index_loose=tuple(map(int,miller_index)) #tuple
        self.gpaw_calc=gpaw_calc
        self.target_dir='results/'+element+'/'+'ads/'
        self.report_location=self.target_dir+self.miller_index_tight+'_grid_results_report.txt' 
        ## TO-DO: need to figure out how to calculate adsorption energy for larger system
        self.gpaw_calc=gpaw_calc
        self.calc_dict=self.gpaw_calc.__dict__['parameters']
        self.ads=ads
        self.all_ads_file_loc=self.target_dir+self.miller_index_tight+'/'+'adsorbates/'+str(self.ads)+'/'
        self.ads_pot_energy=ads_pot_energy
        ##generate report
        self.initialize_report()
        ## connect to opt_slab
        opt_slab=connect('final_database'+'/'+'surf.db').get_atoms(simple_name=self.element+'_'+self.miller_index_tight)
        
        ##start adsorption calculation
        adsorption_energy_dict={}
        init_adsorbates_site_lst=[]
        adsorption_energy_lst=[]
        all_traj_files=glob(self.all_ads_file_loc+'grid/*/input.traj')
        all_gpw_files=glob(self.all_ads_file_loc+'grid/*/slab.gpw')

        ## restart
        if restart_calc==True and len(all_gpw_files)>=1:
            f = paropen(self.report_location,'a')
            parprint('Restarting...',file=f)
            for gpw_file in all_gpw_files:
                location='/'.join(gpw_file.split('/')[:-1])
                parprint('Skipping '+('/'.join(location.split('/')[-2:]))+' adsorption site...',file=f)
                atoms=restart(gpw_file)[0]
                init_adsorbates_site_lst.append(gpw_file.split('/')[-2])
                adsorption_energy=atoms.get_potential_energy()-(opt_slab.get_potential_energy()+self.ads_pot_energy)
                adsorption_energy_lst.append(adsorption_energy)
            f.close()
        
            all_gpw_files_ads_site=['/'.join(i.split('/')[:-1]) for i in all_gpw_files]
            all_traj_files=[i for i in all_traj_files if '/'.join(i.split('/')[:-1]) not in all_gpw_files_ads_site]


        for traj_file in all_traj_files:
            init_adsobates_site, adsorption_energy=self.adsorption_energy_calculator(traj_file,opt_slab)
            init_adsorbates_site_lst.append(init_adsobates_site)
            adsorption_energy_lst.append(adsorption_energy)
        
        adsorption_energy_dict['init_adsorbates_sites[x_y](Ang)']=init_adsorbates_site_lst
        adsorption_energy_dict['adsorption_energy(eV)']=adsorption_energy_lst
        ads_df=pd.DataFrame(adsorption_energy_dict)
        ads_df.set_index('init_adsorbates_sites[x_y](Ang)',inplace=True)
        ads_df.sort_values(by=['adsorption_energy(eV)'],inplace=True)
        f=paropen(self.report_location,'a')
        parprint(ads_df,file=f)
        f.close()
        ads_df.to_csv(self.target_dir+self.miller_index_tight+'_ads_grid.csv')
        f=paropen(self.report_location,'a')
        parprint('Grid adsorption energy calculation complete.',file=f)
        f.close()

    def adsorption_energy_calculator(self,traj_file,opt_slab):
        ads_slab=read(traj_file)
        pbc_checker(ads_slab)
        if self.calc_dict['spinpol']:
            self.apply_magmom(opt_slab,ads_slab)
        fixed_line_constrain=FixedLine(a=-1,direction=[0,0,1])
        slab_c_coord,cluster=detect_cluster(ads_slab)
        if self.fix_option == 'bottom':
            unique_cluster_index=sorted(set(cluster), key=cluster.index)[self.fix_layer-1]
            max_height_fix=max(slab_c_coord[cluster==unique_cluster_index])
            fix_mask=ads_slab.positions[:,2]<(max_height_fix+0.05) #add 0.05 Ang to make sure all bottom fixed
        else:
            raise RuntimeError('Only bottom fix option available now.')
        fixed_atom_constrain=FixAtoms(mask=fix_mask)
        ads_slab.set_constraint([fixed_atom_constrain,fixed_line_constrain])
        ads_slab.set_calculator(self.gpaw_calc)
        location='/'.join(traj_file.split('/')[:-1])
        f=paropen(self.report_location,'a')
        parprint('Calculating '+('/'.join(location.split('/')[-2:]))+' adsorption site...',file=f)
        f.close()
        opt.relax(ads_slab,location,fmax=self.solver_fmax,maxstep=self.solver_max_step)
        init_ads_site=traj_file.split('/')[-2]
        adsorption_energy=ads_slab.get_potential_energy()-(opt_slab.get_potential_energy()+self.ads_pot_energy)
        return init_ads_site, adsorption_energy

    def apply_magmom(self,opt_slab,ads_slab):
        slab_formula=ads_slab.get_chemical_symbols()
        magmom=opt_slab.get_magnetic_moments()
        magmom_ls=np.append(magmom,np.mean(magmom))
        magmom_ls[slab_formula.index(self.ads)]=0
        ads_slab.set_initial_magnetic_moments(magmom_ls)

    def initialize_report(self):
        if world.rank==0 and os.path.isfile(self.report_location):
            os.remove(self.report_location)
        f = paropen(self.report_location,'a')
        parprint('Initial Parameters:', file=f)
        parprint('\t'+'xc: '+self.calc_dict['xc'],file=f)
        parprint('\t'+'h: '+str(self.calc_dict['h']),file=f)
        parprint('\t'+'kpts: '+str(self.calc_dict['kpts']),file=f)
        parprint('\t'+'sw: '+str(self.calc_dict['occupations']),file=f)
        parprint('\t'+'spin polarized: '+str(self.calc_dict['spinpol']),file=f)
        if self.calc_dict['spinpol']:
            parprint('\t'+'magmom: initial magnetic moment from slab calculation.',file=f)
        parprint(' \n',file=f)
        f.close()

##big TO-DO:
    #add adsorption energy calculation for larger surface 
    #both autocat and mapping