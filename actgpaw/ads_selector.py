from gpaw import GPAW,Mixer,MixerDif,Davidson
import glob
import actgpaw.crystal.optimizer as opt
import numpy as np 
from ase.parallel import parprint,world,paropen
import os
from ase.db import connect
import sys
from ase.calculators.calculator import kptdensity2monkhorstpack as kdens2mp
from ase.io import read,write

def ads_auto_select(element,
                struc,
                gpaw_calc,
                ads,
                ads_pot_e,
                temp_print=True,
                size='1x1'):

    #convert str ind to tuple
    m_ind=tuple(map(int,struc))

    #set up the workspace
    code_dir=os.getcwd() #get the current working dir
    #struc_dir=element+'/'+'ads'+'/'+struc #get the structure dir
    #
    
    #create report
    rep_location=code_dir+'/'+element+'/'+'ads'+'/'+struc+'_results_report.txt' 
    if world.rank==0 and os.path.isfile(rep_location):
        os.remove(rep_location)
    #connect to the surface database to get the parameters for calculation
    opt_slab=connect('final_database'+'/'+'surf.db').get_atoms(name=element+'('+struc+')')
    calc_dict=gpaw_calc.__dict__['parameters']
    if calc_dict['spinpol']:
        magmom=opt_slab.get_magnetic_moments()

    with paropen(rep_location,'a') as f:
        parprint('Initial Parameters:',file=f)
        parprint('\t'+'Materials: '+element,file=f)
        parprint('\t'+'Miller Index: '+str(m_ind),file=f)
        parprint('\t'+'Adsorbate: '+str(ads),file=f)
        parprint('\t'+'xc: '+calc_dict['xc'],file=f)
        parprint('\t'+'h: '+str(calc_dict['h']),file=f)
        parprint('\t'+'kpts: '+str(calc_dict['kpts']),file=f)
        parprint('\t'+'sw: '+str(calc_dict['occupations']),file=f)
        parprint('\t'+'spin polarized: '+str(calc_dict['spinpol']),file=f)
        if calc_dict['spinpol']:
            parprint('\t'+'init_magmom: '+str(magmom),file=f)
    f.close()
    
    ads_file_loc=code_dir+'/'+element+'/'+'ads'+'/'+struc
    fils=glob.glob(ads_file_loc+'/'+'adsorbates/Li/**/**/*.traj',recursive=False)
    ads_dict={}
    with paropen(rep_location,'a') as f:
        parprint('Ads Site(Ang)\t\t\tAds Energy(eV)',file=f)
    f.close()
    for file_loc in fils:
        ads_slab=read(file_loc)
        # kpts=kdens2mp(ads_slab,kptdensity=k_density,even=True)
        slab_length=ads_slab.cell.lengths()
        slab_long_short_ratio=max(slab_length)/min(slab_length)
        if calc_dict['spinpol']:
            slab_formula=ads_slab.get_chemical_symbols()
            magmom_ls=np.append(magmom,np.mean(magmom))
            magmom_ls[slab_formula.index(ads)]=0
            ads_slab.set_initial_magnetic_moments(magmom_ls)
        if slab_long_short_ratio > 15:  
            with paropen(rep_location,'a') as f:
                parprint('WARNING: slab long-short side ratio is'+str(slab_long_short_ratio),file=f)
                parprint('Consider change the mixer setting, if not converged.',fiile=f)
            f.close()                           
        ads_slab.set_calculator(gpaw_calc)
        location='/'.join(file_loc.split('/')[:-1])
        opt.surf_relax(ads_slab, location, fmax=0.01, maxstep=0.04, replay_traj=None)
        ads_dict[location]=ads_slab.get_potential_energy()-(opt_slab.get_potential_energy()+ads_pot_e)
        if temp_print:
            with paropen(rep_location,'a') as f:
                parprint(str(file_loc.split('/')[-2])+'\t\t\t'+str(np.round(ads_dict[location],decimals=5)),file=f)
            f.close()
    ads_dict_sorted=sorted(ads_dict,key=ads_dict.get)
    lowest_ads_e_slab=read(ads_dict_sorted[0]+'/slab.traj')
    
    ads_db=connect('final_database/ads'+str(size)+'.db')
    id=ads_db.reserve(name=element+'('+struc+')')
    if id is None:
        id=ads_db.get(name=element+'('+struc+')').id
        ads_db.update(id=id,atoms=lowest_ads_e_slab,name=element+'('+struc+')',clean_slab_pot_e=opt_slab.get_potential_energy(),ads_pot_e=ads_dict[ads_dict_sorted[0]])
    else:
        ads_db.write(lowest_ads_e_slab,id=id,name=element+'('+struc+')',clean_slab_pot_e=opt_slab.get_potential_energy(),ads_pot_e=ads_dict[ads_dict_sorted[0]])
    with paropen(rep_location,'a') as f:
        parprint('Computation Complete. Selected ads site is: '+ads_dict_sorted[0].split('/')[-1],file=f)