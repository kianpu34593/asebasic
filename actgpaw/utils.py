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
from matplotlib import pyplot as plt
from ase.io import read,write
import numpy as np
import pandas as pd

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

def create_element_dir(element,options=['bulk','surf','ads'],
                surf_struc=['100','110','111'],
                optimized_parameters=['h','kdens'],
                starting_layer=3,
                ads_atom=['Li'],
                ads_site=['ontop','hollow','bridge'],
                interval=2):
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
        for struc in surf_struc:
            create_surf_sub_dir(element,struc,starting_layer,interval)
            # create_surf_vac_dir(element,struc,init_vac)
        print('{} surf directories created!'.format(element))

def create_ads_and_dir(element, 
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
        if os.path.isdir(element+'/'+'ads'):
            print("WARNING: {}/ads/ directory already exists!".format(element))
            pause()
        else:
            os.makedirs(element+'/'+'ads',exist_ok=True) 
        if os.path.isdir(element+'/'+'ads'+'/'+struc):
            print('WARNING: '+element+'/'+'ads'+'/'+'{}/ directory already exists!'.format(struc))
            pause()
        else:
            os.makedirs(element+'/'+'ads'+'/'+struc,exist_ok=True)
        surf = surf_db.get_atoms(name=element+'('+struc+')')
        sub_dir=element+'/'+'ads'+'/'+struc
        os.chdir(current_dir+'/'+sub_dir)
        adsorption.generate_rxn_structures(surf,ads=ads_atom,site_type=ads_site,write_to_disk=True)
        os.chdir(current_dir)


def create_surf_sub_dir(element,struc,starting_layer,interval):
    sub_dir=element+'/'+'surf'+'/'+struc
    if os.path.isdir(element+'/'+'surf'+'/'+struc):
        print('WARNING: '+sub_dir+'/ directory already exists!')
        pause()
    else:
        os.makedirs(sub_dir,exist_ok=True)
    for layer in range(starting_layer,starting_layer+6*interval,interval):
        sub_sub_dir=element+'/'+'surf'+'/'+struc+'/'+str(layer)+'x1x1'
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

def sym_all_slab(element,max_ind,layers=10,vacuum_layer=10):
    bulk_ase=connect('final_database/bulk.db').get_atoms(name=element)
    bulk_pym=AseAtomsAdaptor.get_structure(bulk_ase)
    slabgenall=generate_all_slabs(bulk_pym,max_ind,layers,vacuum_layer,
                                center_slab=True,symmetrize=False,in_unit_planes=True)
    print('Miller Index'+'\t'+'Num of Different Shift(s)')
    slab_M=[]
    for slab in slabgenall:
        if slab.is_symmetric():
            slab_M.append([slab.miller_index])
    slab_M_unique = Counter(chain(*slab_M))
    for key in list(slab_M_unique.keys()):
        print(str(key)+'\t'+str(slab_M_unique[key]))

def surf_creator(element,ind,layers,vacuum_layer,unit,shift_to_save,save=False):
    bulk_ase=connect('final_database/bulk.db').get_atoms(name=element)
    bulk_pym=AseAtomsAdaptor.get_structure(bulk_ase)
    slabgen = SlabGenerator(bulk_pym, ind, layers, vacuum_layer,
                            center_slab=True,in_unit_planes=unit)
    slabs_all=slabgen.get_slabs()
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
            angles=np.round(slab_ase.cell.cellpar()[3:],decimals=4)
            if angles[2] != 90.0000:
                L=slab_ase.cell.lengths()[2]
                slab_ase.cell[2]=[0,0,L]
                slab_ase.wrap()
            slab_ase_ls.append(slab_ase)
            angle_ls.append(np.round(slab_ase.cell.cellpar()[3:],decimals=4))
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
            raise RuntimeError('More than one slabs to save! Current code only saves the first one!')
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



