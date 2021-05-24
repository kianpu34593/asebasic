import os
import sys
from pymatgen.io.cif import CifWriter
from pymatgen.ext.matproj import Element, MPRester
from autocat import adsorption
from ase.db import connect
from ase.build import add_vacuum
import copy as cp
def pause():
    input('Press <ENTER> to continue...')

def create_big_dir():
    current_dir=os.getcwd()
    os.chdir(current_dir)
    #create the orig_cif_data and final_database dir
    if os.path.isdir('orig_cif_data'):
        print("WARNING: orig_cif_data directory already exists!")
        pause()
    else:
        os.makedirs('orig_cif_data',exist_ok=True)
    if os.path.isdir('final_database'):
        print("WARNING: final_database directory already exists!")
        pause()
    else:
        os.makedirs('final_database',exist_ok=True)

def create_element_dir(element,options=['bulk','surf','ads'],
                surf_struc=['100','110','111'],
                optimized_parameters=['h','k','sw'],
                starting_layer=3,
                # init_vac=5,
                ads_atom=['Li'],
                ads_site=['ontop','hollow','bridge'],
                interval=2):
    current_dir=os.getcwd()
    os.chdir(current_dir)

    #create the element dir
    if os.path.isdir(element):
        print("WARNING: {} directory already exists!".format(element))
        pause()
    else:
        os.makedirs(element,exist_ok=True)
    
    #create the bulk dir
    if 'bulk' in options:
        if os.path.isdir(element+'/'+'bulk'):
            print("WARNING: ./{}/bulk directory already exists!".format(element))
            pause()
        else:
            os.makedirs(element+'/'+'bulk',exist_ok=True)
        for par in optimized_parameters:
            create_bulk_sub_dir(element,par)
        print("{} bulk directory created!".format(element))

    #create the surf dir
    if 'surf' in options:
        if os.path.isdir(element+'/'+'surf'):
            print("WARNING: ./{}/surf directory already exists!".format(element))
            pause()
        else:
            os.makedirs(element+'/'+'surf',exist_ok=True)
        for struc in surf_struc:
            create_surf_sub_dir(element,struc,starting_layer,interval)
            # create_surf_vac_dir(element,struc,init_vac)
        print('{} surf directories created!'.format(element))

    #create adsorption dir
    if 'ads' in options:
        for struc in surf_struc:
            slab_db_path=element+'/'+'surf'+'/'+struc+'/'+'layer_converge.db'
            if not os.path.isfile(slab_db_path):
                sys.exit("ERROR: slab database has not been established!")
            slab_db=connect(slab_db_path)
            if os.path.isdir(element+'/'+'ads'):
                print("WARNING: ./{}/ads directory already exists!".format(element))
                pause()
            else:
                os.makedirs(element+'/'+'ads',exist_ok=True) 
            if os.path.isdir(element+'/'+'ads'+'/'+struc):
                print('WARNING: '+'./'+element+'/'+'ads'+'/'+'{} directory already exists!'.format(struc))
                pause()
            else:
                os.makedirs(element+'/'+'ads'+'/'+struc,exist_ok=True)            
            create_ads_sub_dir(element,struc,current_dir,ads_site,ads_atom,slab_db)
            os.chdir(current_dir)
            # create_ads_vac_dir(element,struc,slab_clean_row,current_dir,slab_clean,ads_atom,ads_height)
            # os.chdir(current_dir)
        print('{} ads directories created!'.format(element))

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
            print("WARNING: ./{}/ads directory already exists!".format(element))
            pause()
        else:
            os.makedirs(element+'/'+'ads',exist_ok=True) 
        if os.path.isdir(element+'/'+'ads'+'/'+struc):
            print('WARNING: '+'./'+element+'/'+'ads'+'/'+'{} directory already exists!'.format(struc))
            pause()
        else:
            os.makedirs(element+'/'+'ads'+'/'+struc,exist_ok=True)
        surf = surf_db.get_atoms(name=element+'('+struc+')')
        sub_dir=element+'/'+'ads'+'/'+struc
        os.chdir(current_dir+'/'+sub_dir)
        adsorption.generate_rxn_structures(surf,ads=ads_atom,site_type=ads_site,write_to_disk=True)
        os.chdir(current_dir)

def create_ads_sub_dir(element,struc,current_dir,site,ads_atom,slab_db):
    db_size=len(slab_db)
    for i in range(db_size):
        os.chdir(current_dir)
        layer=slab_db.get(i+1).act_layer
        sub_dir=element+'/'+'ads'+'/'+struc+'/'+str(layer)+'x1x1'
        if os.path.isdir(sub_dir):
            print('WARNING: '+'./'+sub_dir+' directory already exists!')
            pause()
        else:
            os.makedirs(sub_dir,exist_ok=True)
        os.chdir(current_dir+'/'+sub_dir)
        #adsorption.gen_rxn_int_sym(slab_db.get_atoms(i+1), ads=[ads_atom],height={ads_atom:ads_height})
        adsorption.generate_rxn_structures(slab_db.get_atoms(i+1),ads=ads_atom,site_type=site,write_to_disk=True)


# def create_surf_vac_dir(element,struc,init_vac):
#     sub_dir=element+'/'+'surf'+'/'+struc+'/'+'layer_optimized'
#     if os.path.isdir(sub_dir):
#         print('WARNING: '+'./'+sub_dir+' directory already exists!')
#         pause()
#     else:
#         os.makedirs(sub_dir,exist_ok=True)
#     for vac in range(init_vac,init_vac+6):
#         sub_sub_dir=element+'/'+'surf'+'/'+struc+'/'+'layer_optimized'+'/'+'vacuum_'+str(vac)
#         if os.path.isdir(sub_sub_dir):
#             print('WARNING: '+'./'+sub_sub_dir+' directory already exists!')
#             pause()
#         else:
#             os.makedirs(sub_sub_dir,exist_ok=True)

def create_surf_sub_dir(element,struc,starting_layer,interval,order=0):
    sub_dir=element+'/'+'surf'+'/'+struc+'_'+str(order)
    if os.path.isdir(sub_dir):
        print('WARNING: '+'./'+sub_dir+' directory already exists!')
        pause()
    else:
        os.makedirs(sub_dir,exist_ok=True)
    for layer in range(starting_layer,starting_layer+6*interval,interval):
        sub_sub_dir=sub_dir+'/'+str(layer)+'x1x1'
        if os.path.isdir(sub_sub_dir):
            print('WARNING: '+'./'+sub_sub_dir+' directory already exists!')
            pause()
        else:
            os.makedirs(sub_sub_dir,exist_ok=True)
        
def create_bulk_sub_dir(element,par):
    sub_dir=element+'/'+'bulk'+'/'+'results'+'_'+par
    if os.path.isdir(sub_dir):
        print('WARNING: '+'./'+sub_dir+' directory already exists!')
        pause()
    else:
        os.makedirs(sub_dir,exist_ok=True)
    sub_sub_dir=element+'/'+'bulk'+'/'+'results'+'_'+par+'/'+'eos_fit'
    if os.path.isdir(sub_sub_dir):
        print('WARNING: '+'./'+sub_sub_dir+' directory already exists!')
    else:
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
