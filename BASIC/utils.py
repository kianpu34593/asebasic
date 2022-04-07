import os
import sys
from pymatgen.io.cif import CifWriter
from pymatgen.core.structure import Structure
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
from typing import List, Type
from glob import glob
import warnings
import itertools
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import fcluster, linkage
from ase import Atom
import matplotlib.pyplot as plt
from ase.visualize.plot import plot_atoms

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
                order_lst: List[int]=None,
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
        for shift,order in zip(shift_lst,order_lst):
            create_surf_sub_dir(element,miller_index,shift,order)
            # create_surf_vac_dir(element,struc,init_vac)
        print('{}/surf/ directories created!'.format(element))

def create_surf_sub_dir(element,miller_index_input,shift,order):
    miller_index=''.join(miller_index_input.split(','))
    miller_index_loose=tuple(map(int,miller_index_input.split(',')))
    raw_surf_dir=element+'/'+'raw_surf'
    if not os.path.isdir(raw_surf_dir):
        raise RuntimeError(raw_surf_dir+' does not exist.')
    else:
        raw_cif_path=element+'/'+'raw_surf/'+str(miller_index)+'/'+str(shift)+'/'+str(order)+'/'+'*.cif'
        raw_cif_files=glob(raw_cif_path)
        assert len(raw_cif_files)==6, 'The size of raw_cif_files is not 6.'
        #cif_files_name=[cif_file.split('/')[-1] for cif_file in raw_cif_files]
        layers=[cif_file.split('/')[-1].split('.')[0] for cif_file in raw_cif_files]
        #layers=[int(name.split('-')[0]) for name in layers_and_shift]
        sub_dir=element+'/'+'surf'+'/'+miller_index+'_'+str(shift)+'_'+str(order)
        if os.path.isdir(sub_dir):
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

def create_ads_and_dir(element, 
                        surf_struc,
                        ads_option,
                        offset,
                        ortho=False,
                        ads_atom=['Li'],
                        ads_site=['ontop','hollow','bridge'],
                        grid_size=[0.5,0.5],
                        slab_size=(1,1,1),
                        tuple_list=[()],
                        height_dict=None,
                        custom_position=[0,0],
                        ):
    current_dir=os.getcwd()
    surf_db_path='final_database/surf.db'

    if not os.path.isfile(surf_db_path):
        sys.exit("ERROR: surf database has not been established!")
    else:
        surf_db=connect(surf_db_path)
    
    primitive_ads_db_path='final_database/ads_1x1.db'
    if os.path.isfile(primitive_ads_db_path):
        ads1x1_db=connect(primitive_ads_db_path)

    for struc in surf_struc:
        primitive_slab = surf_db.get_atoms(simple_name=element+'_'+struc)
        sub_dir='results/'+element+'/'+'ads'+'/'+str(slab_size[0])+'x'+str(slab_size[1])+'/'+struc
        os.makedirs(sub_dir,exist_ok=True)
        
        big_slab=primitive_slab*slab_size
        if ads_option=='autocat':
            os.chdir(current_dir+'/'+sub_dir)
            adsorption.generate_rxn_structures(big_slab,ads=ads_atom,site_type=ads_site,write_to_disk=True,height=height_dict)
        elif ads_option=='grid':
            single_cell_x=primitive_slab.cell.cellpar()[0]
            single_cell_y=primitive_slab.cell.cellpar()[1]
            single_frac_x=1/(single_cell_x//grid_size[0])
            single_frac_y=1/(single_cell_y//grid_size[1])
            if ortho:
                single_cell_x_element=np.array([primitive_slab.cell[0][0],0])*single_frac_x
                single_cell_y_element=np.array([0,primitive_slab.cell[1][1]])*single_frac_y
            else:
                single_cell_x_element=primitive_slab.cell[0][0:2]*single_frac_x
                single_cell_y_element=primitive_slab.cell[1][0:2]*single_frac_y

            ads_sites=[]
            for i, j in itertools.product(list(range(int(single_cell_x//grid_size[0]))), list(range(int(single_cell_y//grid_size[1])))):
                
                single_ads_site=np.round(i*single_cell_x_element+j*single_cell_y_element,decimals=3)
                ###testing
                single_ads_site+=offset
                ####
                ads_sites.append((single_ads_site))
            sites_dict={'grid':ads_sites}
            os.chdir(current_dir+'/'+sub_dir)
            adsorption.generate_rxn_structures(big_slab,ads=ads_atom,all_sym_sites=False,sites=sites_dict,write_to_disk=True,height=height_dict)
        elif ads_option=='lowest_ads_site':
            primitive_ads_slab=ads1x1_db.get_atoms(name=element+'_'+struc)
            primitive_ads_slab.wrap()
            ads_xy_position=np.round(primitive_ads_slab.get_positions()[-1,:2],decimals=3)
            ads_height=primitive_ads_slab.get_positions()[-1,2]-np.max(primitive_ads_slab.get_positions()[:-1,2])
            height_dict={ads_atom[0]:np.round(ads_height,decimals=3)}
            site_dict={'lowest_ads_site':[tuple(ads_xy_position)]}
            os.chdir(current_dir+'/'+sub_dir)
            adsorption.generate_rxn_structures(big_slab,ads=ads_atom,all_sym_sites=False,sites=site_dict,write_to_disk=True,height=height_dict)
        elif ads_option=='custom':
            site_dict={'custom':[tuple(custom_position)]}
            os.chdir(current_dir+'/'+sub_dir)
            adsorption.generate_rxn_structures(big_slab,ads=ads_atom,all_sym_sites=False,sites=site_dict,write_to_disk=True,height=height_dict)
        elif ads_option=='nearest-neighbors':
            big_ads_slab_path = 'final_database/ads_'+str(slab_size[0])+'x'+str(slab_size[1])+'.db'
            big_ads_db = connect(big_ads_slab_path)
            big_ads_slab = big_ads_db.get_atoms(name=element+'_'+struc)
            single_cell_x=primitive_slab.cell.cellpar()[0]
            single_cell_y=primitive_slab.cell.cellpar()[1]
            single_frac_x=1/single_cell_x
            single_frac_y=1/single_cell_y
            single_cell_x_element=primitive_slab.cell[0][0:2]
            single_cell_y_element=primitive_slab.cell[1][0:2]
            ads_xy_position=np.round(big_ads_slab.get_positions()[-1,:2],decimals=3)
            nearest_position_list=[]
            # if single_cell_x != single_cell_y:
            #     if single_cell_x > single_cell_y:
            #         fst_nearst_position[1]+=single_cell_y
            #         snd_nearst_position[0]+=single_cell_x
            #     else:
            #         fst_nearst_position[0]+=single_cell_x
            #         snd_nearst_position[1]+=single_cell_y
            # else: 
            #     fst_nearst_position[1]+=single_cell_y
            #     snd_nearst_position[0]+=single_cell_x
            #     snd_nearst_position[1]+=single_cell_y
            if len(tuple_list)==0:
                raise ValueError('Positions tuple is empty.')
            for i in tuple_list:
                nearst_position=ads_xy_position+single_cell_x_element*i[0]+single_cell_y_element*i[1]
                nearest_position_list.append(nearst_position)
            site_dict={str(i[0])+'_'+str(i[1]):[tuple(j)] for i,j in zip(tuple_list,nearest_position_list)}
            #site_dict={'fst_near':[tuple(fst_nearst_position)],'snd_near':[tuple(snd_nearst_position)]}
            ads_height=big_ads_slab.get_positions()[-1,2]-np.max(big_ads_slab.get_positions()[:-1,2])
            #height_dict={ads_atom[0]:np.round(ads_height,decimals=3)}
            height_dict={ads_atom[0]:0}
            os.chdir(current_dir+'/'+sub_dir)
            adsorption.generate_rxn_structures(big_ads_slab,ads=ads_atom,all_sym_sites=False,sites=site_dict,write_to_disk=True,height=height_dict)
        elif ads_option=='no-adatom':
            os.chdir(current_dir+'/'+sub_dir)
            os.makedirs('clean_slab')
            big_slab.write('clean_slab/input.traj')
            print('clean slab written to ./clean_slab/input.traj')
        else:
            raise TypeError('Specify the ads_option. Availble options: autocat, grid, custom and 2-adatoms')
        os.chdir(current_dir)

def adsobates_plotter(element,
                    miller_index,
                    ads,
                    slab_size=(1,1,1),
                    option='autocat',#grid
                    ):
    current_dir=os.getcwd()
    surf_db_path='final_database/surf.db'
    if not os.path.isfile(surf_db_path):
        sys.exit("ERROR: surf database has not been established!")
    else:
        surf_db=connect(surf_db_path)
    for m_ind in miller_index:
        base_slab = surf_db.get_atoms(simple_name=element+'_'+m_ind)
        base_slab=base_slab*slab_size
        sub_dir='results/'+element+'/'+'ads'+'/'+str(slab_size[0])+'x'+str(slab_size[1])+'/'+m_ind+'/adsorbates/'
        

        if option == 'autocat':
            os.chdir(current_dir+'/'+sub_dir)
            bridges=glob(str(ads)+'/bridge/*/input.traj')
            ontop=glob(str(ads)+'/ontop/*/input.traj')
            hollow=glob(str(ads)+'/hollow/*/input.traj')
            all_files=bridges+ontop+hollow
        elif option == 'grid':
            os.chdir(current_dir+'/'+sub_dir)
            all_files=glob(str(ads)+'/grid/*/input.traj')
        elif option == 'lowest_ads_site':
            os.chdir(current_dir+'/'+sub_dir)
            all_files=glob(str(ads)+'/lowest_ads_site/*/input.traj')
        elif option == 'nearest-neighbors':
            big_ads_slab_path = 'final_database/ads_'+str(slab_size[0])+'x'+str(slab_size[1])+'.db'
            big_ads_db = connect(big_ads_slab_path)
            base_slab = big_ads_db.get_atoms(name=element+'_'+m_ind)
            os.chdir(current_dir+'/'+sub_dir)
            all_files=glob(str(ads)+'/1_0/*/input.traj')+glob(str(ads)+'/0_1/*/input.traj')+glob(str(ads)+'/1_1/*/input.traj')+glob(str(ads)+'/0.5_0.5/*/input.traj')
        else:
            raise TypeError('Specify the option. Availble options: autocat, grid, custom and 2-adatoms')

        for file in all_files:
            slab=read(file)
            positions=slab.get_positions()
            ads_atom_index=[-1]
            Li_position=positions[ads_atom_index,:][0]
            base_slab.append(Atom('O',position=Li_position))
        fig, axarr = plt.subplots(1, 3, figsize=(15, 5))
        plot_atoms(base_slab,axarr[0],rotation=('0x,0y,0z'))
        plot_atoms(base_slab,axarr[1],rotation=('270x,0y,0z'))
        plot_atoms(base_slab,axarr[2],rotation=('270x,90y,0z'))
        fig.savefig("ads_sites.png")
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

def sym_all_slab(element,max_ind,layers=5,vacuum_layer=10,symmetric=False):
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

def surf_creator(element,ind,layers,vacuum_layer,unit,order_to_save,save=False,orthogonalize=False,symmetric=False):
    bulk_ase=connect('final_database/bulk.db').get_atoms(name=element)
    bulk_pym=AseAtomsAdaptor.get_structure(bulk_ase)
    slabgen = SlabGenerator(bulk_pym, ind, layers, vacuum_layer,
                            center_slab=True,in_unit_planes=unit)

    slabs_all=slabgen.get_slabs(symmetrize=symmetric)

    slabs_symmetric=slabs_all
    
    #slabs_symmetric=[slabs_all[i] for i, slab in enumerate(slabs_all) if slab.is_symmetric()]
    if len(slabs_symmetric) == 0:
        raise RuntimeError('No symmetric slab found!')
    else:
        shift_ls=[]
        slab_ase_ls=[]
        angle_ls=[]
        num_different_layers_ls=[]
        num_atom_ls=[]
        composition_ls=[]
        shift_ls = [] 
       

        for i,slab in enumerate(slabs_symmetric):
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
                slab_ase.center()
            slab_ase_ls.append(slab_ase)
            angle_ls.append(np.round(slab_ase.cell.angles(),decimals=4))
            shift_ls.append(np.round(slab.shift,decimals=4))
            unique_cluster=np.unique(detect_cluster(slab_ase)[1])
            num_different_layers_ls.append(len(unique_cluster))
            num_atom_ls.append(len(slab_ase))
            composition_dict=dict(Counter(slab_ase.get_chemical_symbols()))
            total_num_atoms=len(slab_ase)
            composition_ls.append({key: np.round(values/total_num_atoms,decimals=4) for key, values in composition_dict.items()})
        if len(slabs_symmetric)==len(slabgen._calculate_possible_shifts()):
            shift_ls=np.round(slabgen._calculate_possible_shifts(),decimals=4)

        slabs_info_dict={'shift':shift_ls,'angles':angle_ls,'actual_layers':num_different_layers_ls,'num_of_atoms':num_atom_ls,'composition':composition_ls}
        slabs_info_df=pd.DataFrame(slabs_info_dict)
        print(slabs_info_df)
        if os.path.isfile(surf_location):
            os.remove(surf_location)
    if save:
        #slab_order_save=[i for i,slab in enumerate(slabs_symmetric) if np.round(slab.shift,decimals=4)==shift_to_save]
        if len(slab_ase_ls)==0:
            raise RuntimeError('No slab to save!')
        #elif len(slab_order_save)>1:
            #warnings.warn('More than one slabs to save! Current code only saves the first one!')

        if order_to_save>len(slab_ase_ls)-1:
            raise RuntimeError('order_to_save exceeds the number of slabs!')
        surf_saver(element,slab_ase_ls[order_to_save],ind,layers,shift_ls[order_to_save],order_to_save)

def surf_saver(element,slab_to_save,ind,layers,shift,order_to_save):
    tight_ind=''.join(list(map(lambda x:str(x),ind)))
    rep_location='results/'+element+'/raw_surf/'+str(tight_ind)+'/'+str(shift)+'/'+str(order_to_save)
    os.makedirs(rep_location,exist_ok=True)
    surf_location='results/'+element+'/raw_surf/'+str(tight_ind)+'/'+str(shift)+'/'+str(order_to_save)+'/'+str(layers)+'.cif'
    if os.path.isfile(surf_location):
        raise RuntimeError(surf_location+' already exists!')
    else:
        slab_to_save.write(surf_location,format='cif')
        print('Raw surface saving complete!')

def detect_cluster(slab,tol=0.3):
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


