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
            for gpw_file in all_gpw_files:
                atoms=restart(gpw_file)[0]
                init_adsorbates_site_lst.append(gpw_file.split('/')[-2])
                adsorption_energy=atoms.get_potential_energy()-(opt_slab.get_potential_energy()+self.ads_pot_energy)
                adsorption_energy_lst.append(adsorption_energy)
                final_adsorbates_site='_'.join(list(np.round(atoms.get_positions()[-1][:2],decimals=3)))
                final_adsorbates_site_lst.append(final_adsorbates_site)
            all_gpw_files_ads_site=['/'.join(i.split('/')[:-1]) for i in all_gpw_files]
            all_traj_files=[i for i in all_traj_files if '/'.join(i.split('/')[:-1]) not in all_gpw_files_ads_site]
        
        for traj_file in all_traj_files:
            init_adsobates_site, adsorption_energy, final_adsorbates_site=self.adsorption_energy_calculator(traj_file,opt_slab)
            init_adsorbates_site_lst.append(init_adsobates_site)
            adsorption_energy_lst.append(adsorption_energy)
            final_adsorbates_site_lst.append(final_adsorbates_site)
        
        adsorption_energy_dict['init_adsorbates_sites[x_y](Ang)']=init_adsorbates_site_lst
        adsorption_energy_dict['adsorption_energy(eV)']=adsorption_energy_lst
        adsorption_energy_dict['final_adsorbates_sites[x_y](Ang)']=final_adsorbates_site_lst
        ads_df=pd.DataFrame(adsorption_energy_dict)
        ads_df.set_index('adsorbates_sites[x_y](Ang)',inplace=True)
        ads_df.sort_values(by=['adsorption_energy(eV)'],inplace=True)
        f=paropen(self.report_location,'a')
        parprint(ads_df,file=f)
        f.close()
        min_adsorbates_site=ads_df.idxmin()['adsorbates_sites[x_y](Ang)']
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
        ads_slab.set_calculator(self.gpaw_calc)
        location='/'.join(traj_file.split('/')[:-1])
        f=paropen(self.report_location,'a')
        parprint('Calculating '+traj_file.split('/')[-3]+' '+traj_file.split('/')[-2]+' adsorption site...',file=f)
        f.close()
        opt.relax(ads_slab,location,fmax=self.solver_fmax,maxstep=self.solver_max_step)
        init_ads_site=traj_file.split('/')[-2]
        adsorption_energy=ads_slab.get_potential_energy()-(opt_slab.get_potential_energy()+self.ads_pot_energy)
        final_ads_site='_'.join(list(np.round(ads_slab.get_positions()[-1][:2],decimals=3)))
        return init_ads_site, adsorption_energy, final_ads_site

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
            for gpw_file in all_gpw_files:
                atoms=restart(gpw_file)[0]
                init_adsorbates_site_lst.append(gpw_file.split('/')[-2])
                adsorption_energy=atoms.get_potential_energy()-(opt_slab.get_potential_energy()+self.ads_pot_energy)
                adsorption_energy_lst.append(adsorption_energy)
            all_gpw_files_ads_site=['/'.join(i.split('/')[:-1]) for i in all_gpw_files]
            all_traj_files=[i for i in all_traj_files if '/'.join(i.split('/')[:-1]) not in all_gpw_files_ads_site]
        
        for traj_file in all_traj_files:
            init_adsobates_site, adsorption_energy=self.adsorption_energy_calculator(traj_file,opt_slab)
            init_adsorbates_site_lst.append(init_adsobates_site)
            adsorption_energy_lst.append(adsorption_energy)
        
        adsorption_energy_dict['init_adsorbates_sites[x_y](Ang)']=init_adsorbates_site_lst
        adsorption_energy_dict['adsorption_energy(eV)']=adsorption_energy_lst
        ads_df=pd.DataFrame(adsorption_energy_dict)
        ads_df.set_index('adsorbates_sites[x_y](Ang)',inplace=True)
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
        mask=FixedLine(a=len(ads_slab)-1,direction=[0,0,1])
        ads_slab.set_constraint(mask)
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


class surf_calc:
    def __init__(self,
                element: str,
                miller_index: str,
                pymatgen_layer: int,
                shift: float,
                gpaw_calc,
                rela_tol: float=0.015,
                #restart_calc: bool=False,
                fix_layer: int=2,
                vacuum: int=10,
                solver_fmax: float=0.01,
                solver_max_step: float=0.05,
                surf_energy_calc_mode: str='regular', #because only calculate one 
                fix_option: str='bottom'):
        #intialize
        ##globalize variables
        self.element=element
        self.solver_max_step=solver_max_step
        self.solver_fmax=solver_fmax
        self.surf_energy_calc_mode=surf_energy_calc_mode
        self.vacuum=vacuum
        self.fix_option=fix_option
        self.fix_layer=fix_layer
        self.miller_index_tight=miller_index
        self.miller_index_loose=tuple(map(int,miller_index)) #tuple
        self.shift=shift
        self.gpaw_calc=gpaw_calc
        self.final_slab_name=self.element+'_'+self.miller_index_tight+'_'+str(self.shift)
        self.raw_slab_dir='results/'+element+'/'+'raw_surf/'
        self.target_dir='results/'+element+'/'+'surf/'
        self.target_sub_dir=self.target_dir+self.miller_index_tight+'_'+str(self.shift)+'/'
        self.report_location=(self.target_dir+self.miller_index_tight+'_'+str(self.shift)+'_results_report.txt')
        self.rela_tol = rela_tol

        ##connect to optimize bulk database to get gpw_dir and bulk potential_energy
        db_bulk=connect('final_database/bulk.db')
        kdensity=db_bulk.get(name=self.element).kdensity
        self.bulk_potential_energy=(db_bulk.get_atoms(name=self.element).get_potential_energy())/len(db_bulk.get_atoms(name=element))
        
        ##read the smallest slab to get the kpoints
        self.ascend_all_cif_files_full_path=self.sort_raw_slab()
        raw_slab_smallest=read(self.ascend_all_cif_files_full_path[0])
        raw_slab_smallest.pbc=[1,1,0]
        kpts=kdens2mp(raw_slab_smallest,kptdensity=kdensity,even=True)
        self.gpaw_calc.__dict__['parameters']['kpts']=kpts
        self.calc_dict=self.gpaw_calc.__dict__['parameters']

        ##generate report
        if self.calc_dict['spinpol']:
            self.init_magmom=np.mean(db_bulk.get_atoms(name=element).get_magnetic_moments())
        self.initialize_report()


        # convergence test 

        ## number of layers
        ### restart 
        if restart_calc and len(glob(self.target_sub_dir+'*/*.gpw'))>2:
            ascend_layer_ls,ascend_gpw_files_dir=self.gather_gpw_file()
            if len(ascend_gpw_files_dir) > 2:
                for i in range((len(ascend_layer_ls)-3)+1):
                    self.convergence_update(i,ascend_gpw_files_dir)
                    diff_primary=max(self.surf_energies_diff_arr[0],self.surf_energies_diff_arr[2])
                    diff_second=self.surf_energies_diff_arr[1]
        else:
            #os.remove(self.target_dir+self.miller_index_tight+'_'+str(self.shift)+'/'+)
            ascend_layer_ls=[]
            diff_primary=100
            diff_second=100
        iters=len(ascend_layer_ls)
        self.convergence_loop(iters,diff_primary,diff_second)

        #finalize
        ascend_gpw_files_dir=self.gather_gpw_file()[1]
        ## calculate the surface energy
        if self.surf_energy_calc_mode == 'regular':
            final_atoms,self.gpaw_calc=restart(ascend_gpw_files_dir[-3])
            slab_energy=[final_atoms.get_potential_energy()]
            surface_area=[2*final_atoms.cell[0][0]*final_atoms.cell[1][1]]
            num_of_atoms=[len(final_atoms)]
            surf_energy=np.round(self.surface_energy_calculator(np.array(slab_energy),np.array(surface_area),np.array(num_of_atoms))[0],decimals=4)
            self.calc_dict=self.gpaw_calc.__dict__['parameters']
        elif self.surf_energy_calc_mode == 'linear-fit':
            slab_energy_lst=[]
            surface_area_total_lst=[]
            num_of_atoms_lst=[]
            for gpw_file_dir in ascend_gpw_files_dir[-3:]:
                interm_atoms=restart(gpw_file_dir)[0]
                slab_energy_lst.append(interm_atoms.get_potential_energy())
                surface_area_total_lst.append(2*interm_atoms.cell[0][0]*interm_atoms.cell[1][1])
                num_of_atoms_lst.append(len(interm_atoms))
            surf_energy=np.round(self.surface_energy_calculator(np.array(slab_energy_lst),np.array(surface_area_total_lst),np.array(num_of_atoms_lst))[0],decimals=4)
            final_atoms,self.gpaw_calc=restart(ascend_gpw_files_dir[-3])
            self.calc_dict=self.gpaw_calc.__dict__['parameters']
        else:
            raise RuntimeError(self.surf_energy_calc_mode+'mode not avilable. Available modes are regular, linear-fit.')
        
        ##save to database
        db_slab_interm=connect(self.target_dir+'all_miller_indices_all_shift'+'.db')
        id=db_slab_interm.reserve(name=self.final_slab_name)
        if id is None:
            id=db_slab_interm.get(name=self.final_slab_name).id
            db_slab_interm.update(id=id,atoms=final_atoms,name=self.final_slab_name,
                                    surf_energy=surf_energy,
                                    kpts=str(','.join(map(str, self.calc_dict['kpts']))))
        else:
            db_slab_interm.write(final_atoms,id=id,name=self.final_slab_name,
                                    surf_energy=surf_energy,
                                    kpts=str(','.join(map(str, self.calc_dict['kpts']))))
        f = paropen(self.report_location,'a')
        parprint('Surface energy calculation complete.', file=f)
        f.close()

    def convergence_loop(self,iters,diff_p,diff_s):
        while (diff_p>self.rela_tol or diff_s>self.rela_tol) and iters <= 6:
            slab=read(self.ascend_all_cif_files_full_path[iters])
            pbc_checker(slab)
            slab.center(vacuum=self.vacuum,axis=2)
            if self.calc_dict['spinpol']:
                slab.set_initial_magnetic_moments(self.init_magmom*np.ones(len(slab)))
            slab_c_coord,cluster=detect_cluster(slab)
            if self.fix_option == 'bottom':
                unique_cluster_index=sorted(set(cluster), key=cluster.index)[self.fix_layer-1]
                max_height_fix=max(slab_c_coord[cluster==unique_cluster_index])
                fix_mask=slab.positions[:,2]<(max_height_fix+0.05) #add 0.05 Ang to make sure all bottom fixed
            else:
                raise RuntimeError('Only bottom fix option available now.')
            slab.set_constraint(FixAtoms(mask=fix_mask))
            slab.set_calculator(self.gpaw_calc)
            layer=self.ascend_all_cif_files_full_path[iters].split('/')[-1].split('_')[-1].split('-')[0]
            location=self.target_sub_dir+layer+'x1x1'
            opt.relax(slab,location,fmax=self.solver_fmax,maxstep=self.solver_max_step)
            ascend_layer_ls,ascend_gpw_files_dir=self.gather_gpw_file()
            iters=len(ascend_layer_ls)
            if iters>2:
                iter=iters-3
                self.convergence_update(iter,ascend_gpw_files_dir)
                diff_p=max(self.surf_energies_diff_arr[0],self.surf_energies_diff_arr[2])
                diff_s=self.surf_energies_diff_arr[1]
        self.check_convergence(diff_p,diff_s,iters)
    
    def check_convergence(self,diff_p,diff_s,iters):
        if iters>=6:
            if diff_p>self.rela_tol or diff_s>self.rela_tol:
                f=paropen(self.report_location,'a')
                parprint("WARNING: Max iterations reached! layer convergence test failed.",file=f)
                parprint("Computation Suspended!",file=f)
                parprint(' ',file=f)
                f.close()
                sys.exit()
        else:
            f=paropen(self.report_location,'a')
            parprint("layer convergence test success!",file=f)
            parprint("="*44,file=f)
            parprint('\n',file=f)
            f.close() 

    def convergence_update(self,iter,gpw_files_dir):
        slab_energy_lst=[]
        num_of_atoms_lst=[]
        surface_area_total_lst=[]
        pymatgen_layer_ls=[]
        for i in range(iter,iter+3,1):
            atoms=restart(gpw_files_dir[i])[0]
            slab_energy_lst.append(atoms.get_potential_energy())
            surface_area_total_lst.append(2*atoms.cell[0][0]*atoms.cell[1][1])
            num_of_atoms_lst.append(len(atoms))
            pymatgen_layer_ls.append(int(gpw_files_dir[i].split('/')[-2].split('x')[0]))
        surf_energy_lst=self.surface_energy_calculator(np.array(slab_energy_lst),np.array(surface_area_total_lst),np.array(num_of_atoms_lst))
        surf_energy_arr=np.array(surf_energy_lst)
        surf_energy_arr_rep= np.array((surf_energy_lst+surf_energy_lst)[1:4])
        self.surf_energies_diff_arr=np.round(np.abs(surf_energy_arr-surf_energy_arr_rep),decimals=4)
        self.convergence_update_report(pymatgen_layer_ls)


    def convergence_update_report(self,layer_ls):
        f = paropen(self.report_location,'a')
        parprint('Optimizing parameter: '+'layers',file=f)
        param_val_str='1st: '+str(layer_ls[0])+' 2nd: '+str(layer_ls[1])+' 3rd: '+str(layer_ls[2])
        parprint('\t'+param_val_str,file=f)
        divider_str='-'
        parprint('\t'+divider_str*len(param_val_str),file=f)
        substrat_str='| '+'2nd-1st'+' | '+'3rd-2nd'+' | '+'3rd-1st'+' |'
        parprint('\t'+substrat_str,file=f)
        energies_str='\t'+'| '
        for i in range(3):
            energies_str+=str(self.surf_energies_diff_arr[i])+'  '+'|'+' '
        energies_str+='eV/Ang^2'
        parprint(energies_str,file=f)
        parprint(' ',file=f)
        f.close()

    def surface_energy_calculator(self,slab_energy,surface_area_total,num_of_atoms):
        if self.surf_energy_calc_mode=='regular':
            surf_energy_lst=(1/surface_area_total)*(slab_energy-num_of_atoms*self.bulk_potential_energy)
            # for slab_energy,surface_area,num_of_atoms in zip(slab_energies,surface_area_total_lst,num_of_atoms_lst):
            #     surf_energy=(1/surface_area)*(slab_energy-num_of_atoms*self.bulk_potential_energy)
            #     surf_energy_lst.append(surf_energy)
        elif self.surf_energy_calc_mode=='linear-fit': ## TO-DO: need to think about how to fit to all slab energies, right now this is localize fitting
            assert type(num_of_atoms)==type(np.array([1,2,3])), 'In linear-fit mode, the type of num_of_atoms variable must be numpy.ndarray'
            assert type(slab_energy)==type(np.array([1,2,3])), 'In linear-fit mode, the type of slab_energy variable must be numpy.ndarray'
            assert len(num_of_atoms)==3, 'In linear-fit mode, the size of num_of_atoms variable must be 3'
            assert len(slab_energy)==3, 'In linear-fit mode, the size of slab_energy variable must be 3'
            self.fitted_bulk_potential_energy=np.round(np.polyfit(num_of_atoms,slab_energy,1)[0],decimals=5)
            surf_energy_lst=(1/surface_area_total)*(slab_energy-num_of_atoms*self.fitted_bulk_potential_energy)
            # for slab_energy,surface_area,num_of_atoms in zip(slab_energies,surface_area_total_lst,num_of_atoms_lst):
            #     surf_energy=(1/surface_area)*(slab_energy-num_of_atoms*self.fitted_bulk_potential_energy)
            #     surf_energy_lst.append(surf_energy)
        else:
            raise RuntimeError(self.surf_energy_calc_mode+'mode not avilable. Available modes are regular, linear-fit.')
        return list(surf_energy_lst)

    def gather_gpw_file(self):
        #need to make sure there are no gpw files from previous run
        gpw_files_dir=glob(self.target_sub_dir+'*/*.gpw')
        gpw_slab_size=[gpw_file.split('/')[-2] for gpw_file in gpw_files_dir]
        slab_layers=[int(i.split('x')[0]) for i in gpw_slab_size]
        ascend_order=np.argsort(slab_layers)
        ascend_gpw_files_dir=[gpw_files_dir[i] for i in ascend_order]
        ascend_param_ls=np.sort(slab_layers)
        return ascend_param_ls,ascend_gpw_files_dir
        
    def sort_raw_slab(self):
        all_cif_files_full_path=glob(self.raw_slab_dir+str(self.miller_index_loose)+'_*'+'-'+str(self.shift)+'.cif')

        cif_files_name=[cif_file.split('/')[-1] for cif_file in all_cif_files_full_path]
        layers_and_shift=[name.split('_')[1] for name in cif_files_name]
        layers=[int(name.split('-')[0]) for name in layers_and_shift]
        ascend_order=np.argsort(layers)
        ascend_all_cif_files_full_path=[all_cif_files_full_path[i] for i in ascend_order]
        return ascend_all_cif_files_full_path

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
            parprint('\t'+'magmom: '+str(self.init_magmom),file=f)
        parprint('\t'+'convergence tolerance: '+str(self.rela_tol)+'eV/Ang^2',file=f)
        parprint('\t'+'surface energy calculation mode: '+str(self.surf_energy_calc_mode),file=f)
        parprint('\t'+'fixed layers: '+str(self.fix_layer),file=f)
        parprint('\t'+'fixed option: '+str(self.fix_option),file=f)
        parprint(' \n',file=f)
        f.close()