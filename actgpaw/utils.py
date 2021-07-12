import os
import sys
from pymatgen.io.cif import CifWriter
from pymatgen.ext.matproj import MPRester
from autocat import adsorption
from ase.db import connect
from pymatgen.core.surface import SlabGenerator, generate_all_slabs
from pymatgen.io.ase import AseAtomsAdaptor
from collections import Counter 
from itertools import chain 
from ase.io import read,write
import numpy as np
import pandas as pd
from typing import List
from glob import glob
import warnings
import itertools

from scipy.cluster.hierarchy import ClusterNode

def pause():
    input('Press <ENTER> to continue...')

def create_big_dir():
    current_dir=os.getcwd()
    os.chdir(current_dir)
    #create the orig_cif_data and final_database dir
    if os.path.isdir('orig_cif_data'):
        print("WARNING: orig_cif_data/ directory already exists!")
        pause()
    else:
        os.makedirs('orig_cif_data',exist_ok=True)
    if os.path.isdir('final_database'):
        print("WARNING: final_database/ directory already exists!")
        pause()
    else:
        os.makedirs('final_database',exist_ok=True)
    if os.path.isdir('results'):
        print("WARNING: results/ directory already exists!")
        pause()
    else:
        os.makedirs('results',exist_ok=True)

def create_element_dir(element,
                miller_index=None,
                shift_lst: List[float]=None,
                options=['bulk','surf'],
                optimized_parameters=['h','kdens']):
    current_dir=os.getcwd()
    os.chdir(current_dir)
    element='results/'+element

    #create the element dir
    if os.path.isdir(element):
        print("WARNING: {}/ directory already exists!".format(element))
        pause()
    else:
        os.makedirs(element,exist_ok=True)
    
    #create the bulk dir
    if 'bulk' in options:
        if os.path.isdir(element+'/'+'bulk'):
            print("WARNING: {}/bulk/ directory already exists!".format(element))
            pause()
        else:
            os.makedirs(element+'/'+'bulk',exist_ok=True)
        for par in optimized_parameters:
            create_bulk_sub_dir(element,par)
        print("{}/bulk/ directory created!".format(element))

    #create the surf dir
    if 'surf' in options:
        if os.path.isdir(element+'/'+'surf'):
            print("WARNING: {}/surf/ directory already exists!".format(element))
            pause()
        else:
            os.makedirs(element+'/'+'surf',exist_ok=True)
        for shift in shift_lst:
            create_surf_sub_dir(element,miller_index,shift)
            # create_surf_vac_dir(element,struc,init_vac)
        print('{}/surf/ directories created!'.format(element))

def create_surf_sub_dir(element,miller_index_input,shift):
    miller_index=''.join(miller_index_input.split(','))
    miller_index_loose=tuple(map(int,miller_index_input.split(',')))
    raw_surf_dir=element+'/'+'raw_surf'
    if not os.path.isdir(raw_surf_dir):
        raise RuntimeError(raw_surf_dir+' does not exist.')
    else:
        raw_cif_path=element+'/'+'raw_surf/'+str(miller_index_loose)+'_*'+'-'+str(shift)+'.cif'
        raw_cif_files=glob(raw_cif_path)
        assert len(raw_cif_files)==6, 'The size of raw_cif_files is not 6.'
        cif_files_name=[cif_file.split('/')[-1] for cif_file in raw_cif_files]
        layers_and_shift=[name.split('_')[1] for name in cif_files_name]
        layers=[int(name.split('-')[0]) for name in layers_and_shift]
        sub_dir=element+'/'+'surf'+'/'+miller_index+'_'+str(shift)
        if os.path.isdir(element+'/'+'surf'+'/'+sub_dir):
            print('WARNING: '+sub_dir+'/ directory already exists!')
            pause()
        else:
            os.makedirs(sub_dir,exist_ok=True)
        for layer in layers:
            sub_sub_dir=sub_dir+'/'+str(layer)+'x1x1'
            os.makedirs(sub_sub_dir,exist_ok=True)
        
def create_bulk_sub_dir(element,par):
    sub_dir=element+'/'+'bulk'+'/'+'results'+'_'+par
    if os.path.isdir(sub_dir):
        print('WARNING: '+sub_dir+'/ directory already exists!')
        pause()
    else:
        os.makedirs(sub_dir,exist_ok=True)
    sub_sub_dir=element+'/'+'bulk'+'/'+'results'+'_'+par+'/'+'eos_fit'
    os.makedirs(sub_sub_dir,exist_ok=True)

def create_ads_and_dir_sym_sites(element, 
                        surf_struc,
                        ads_atom=['Li'],
                        ads_site=['ontop','hollow','bridge']):
    current_dir=os.getcwd()
    surf_db_path='final_database/surf.db'
    os.chdir(current_dir)
    if not os.path.isfile(surf_db_path):
        sys.exit("ERROR: surf database has not been established!")
    else:
        surf_db=connect(surf_db_path)
    for struc in surf_struc:
        os.makedirs('results/'+element+'/'+'ads',exist_ok=True) 
        os.makedirs('results/'+element+'/'+'ads'+'/'+'autocat',exist_ok=True)
        os.makedirs('results/'+element+'/'+'ads'+'/'+'autocat'+'/'+struc,exist_ok=True)
        slab = surf_db.get_atoms(simple_name=element+'_'+struc)
        sub_dir='results/'+element+'/'+'ads'+'/'+'autocat'+'/'+struc
        os.chdir(current_dir+'/'+sub_dir)
        adsorption.generate_rxn_structures(slab,ads=ads_atom,site_type=ads_site,write_to_disk=True)
        os.chdir(current_dir)

def create_ads_and_dir_grid_sites(element,
                                surf_struc,
                                ads_atom=['Li']):
    current_dir=os.getcwd()
    surf_db_path='final_database/surf.db'
    os.chdir(current_dir)
    if not os.path.isfile(surf_db_path):
        sys.exit("ERROR: surf database has not been established!")
    else:
        surf_db=connect(surf_db_path)
    for struc in surf_struc:
        os.makedirs('results/'+element+'/'+'ads',exist_ok=True) 
        os.makedirs('results/'+element+'/'+'ads'+'/'+'grid',exist_ok=True)
        os.makedirs('results/'+element+'/'+'ads'+'/'+'grid'+'/'+struc,exist_ok=True)
        slab = surf_db.get_atoms(simple_name=element+'_'+struc)
        slab_cell_x = slab.cell[0][0]
        slab_cell_y = slab.cell[1][1]
        ads_sites=[]
        for i, j in itertools.product(list(range(1,int(slab_cell_x//0.5+1))), list(range(1,int(slab_cell_y//0.5+1)))):
            ads_sites.append((i*0.5,j*0.5))
        sites_dict={'custom':ads_sites}
        sub_dir='results/'+element+'/'+'ads'+'/'+'grid'+'/'+struc
        os.chdir(current_dir+'/'+sub_dir)
        adsorption.generate_rxn_structures(slab,ads=ads_atom,all_sym_sites=False,sites=sites_dict,write_to_disk=True)
        os.chdir(current_dir)

def cif_grabber(API_key,pretty_formula):
    #currently will grab the cif of the lowest formation_energy_per_atom
    mpr=MPRester(str(API_key))
    doc_ls=mpr.query({'pretty_formula':pretty_formula},['material_id','formation_energy_per_atom'])
    form_dict={}
    for doc in doc_ls:
        form_dict[doc['material_id']]=doc['formation_energy_per_atom']
    id_sorted = sorted(form_dict,key=form_dict.get)
    lowest_matID=id_sorted[0]
    struc=mpr.get_structure_by_material_id(lowest_matID,final=True,conventional_unit_cell=True)
    Cif_temp=CifWriter(struc)
    name='orig_cif_data'+'/'+pretty_formula+'_'+lowest_matID
    Cif_temp.write_file('{}.cif'.format(name))

def sym_all_slab(element,max_ind,layers=10,vacuum_layer=10,symmetric=False):
    bulk_ase=connect('final_database/bulk.db').get_atoms(name=element)
    bulk_pym=AseAtomsAdaptor.get_structure(bulk_ase)
    slabgenall=generate_all_slabs(bulk_pym,max_ind,layers,vacuum_layer,
                                center_slab=True,symmetrize=symmetric,in_unit_planes=True)
    print('Miller Index'+'\t'+'Num of Different Shift(s)'+'\t'+'Shifts')
    slab_M=[]
    slabgenall_sym=[]
    for slab in slabgenall:
        if slab.is_symmetric():
            slab_M.append([slab.miller_index])
            slabgenall_sym.append(slab)
    slab_M_unique = Counter(chain(*slab_M))
    for key in list(slab_M_unique.keys()):
            print(str(key)+'\t'+str(slab_M_unique[key])+'\t\t\t\t'+str([np.round(slab.shift,decimals=4) for slab in slabgenall_sym if slab.miller_index==key]))

def surf_creator(element,ind,layers,vacuum_layer,unit,shift_to_save,save=False,orthogonalize=False,symmetric=False):
    bulk_ase=connect('final_database/bulk.db').get_atoms(name=element)
    bulk_pym=AseAtomsAdaptor.get_structure(bulk_ase)
    slabgen = SlabGenerator(bulk_pym, ind, layers, vacuum_layer,
                            center_slab=True,in_unit_planes=unit)
    slabs_all=slabgen.get_slabs(symmetrize=symmetric)
    slabs_symmetric=[slabs_all[i] for i, slab in enumerate(slabs_all) if slab.is_symmetric()]
    if len(slabs_symmetric) == 0:
        raise RuntimeError('No symmetric slab found!')
    else:
        shift_ls=[]
        slab_ase_ls=[]
        angle_ls=[]
        for slab in slabs_symmetric:
            #temp save for analysis
            os.makedirs('results/'+element+'/raw_surf',exist_ok=True)
            surf_location='results/'+element+'/raw_surf/'+str(ind)+'_temp'+'.cif'
            CifWriter(slab).write_file(surf_location)
            slab_ase=read(surf_location)
            angles=np.round(slab_ase.cell.angles(),decimals=4)
            anlges_arg=[angle != 90.0000 for angle in angles[:2]]
            if orthogonalize==True and np.any(anlges_arg):
                L=slab_ase.cell.lengths()[2]
                slab_ase.cell[2]=[0,0,L]
                slab_ase.wrap()
            slab_ase_ls.append(slab_ase)
            angle_ls.append(np.round(slab_ase.cell.angles(),decimals=4))
            shift_ls.append(np.round(slab.shift,decimals=4))
        slabs_info_dict={'shift':shift_ls,'angles':angle_ls}
        slabs_info_df=pd.DataFrame(slabs_info_dict)
        print(slabs_info_df)
        if os.path.isfile(surf_location):
            os.remove(surf_location)
    if save:
        slab_order_save=[i for i,slab in enumerate(slabs_symmetric) if np.round(slab.shift,decimals=4)==shift_to_save]
        if len(slab_order_save)==0:
            raise RuntimeError('No slab to save!')
        elif len(slab_order_save)>1:
            warnings.warn('More than one slabs to save! Current code only saves the first one!')
        surf_saver(element,slab_ase_ls[slab_order_save[0]],ind,layers,shift_ls[slab_order_save[0]])

def surf_saver(element,slab_to_save,ind,layers,shift):
    rep_location='results/'+element+'/raw_surf'
    os.makedirs(rep_location,exist_ok=True)
    surf_location='results/'+element+'/raw_surf/'+str(ind)+'_'+str(layers)+'-'+str(shift)+'.cif'
    if os.path.isfile(surf_location):
        raise RuntimeError(surf_location+' already exists!')
    else:
        slab_to_save.write(surf_location,format='cif')
        print('Raw surface saving complete!')



