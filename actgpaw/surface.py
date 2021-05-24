from pymatgen.core.surface import SlabGenerator, generate_all_slabs
from pymatgen.io.ase import AseAtomsAdaptor
from ase.db import connect
from ase.build import surface
from pymatgen.io.cif import CifWriter
import numpy as np
from collections import Counter 
from itertools import chain 
from pymatgen.analysis.adsorption import plot_slab
from matplotlib import pyplot as plt
from ase.visualize.plot import plot_atoms
import os
from ase.io import read

def sym_all_slab(element,max_ind,layers,vacuum_layer,opt=True):
    if opt == True:
        bulk_ase=connect('final_database/bulk.db').get_atoms(name=element)
        bulk_pym=AseAtomsAdaptor.get_structure(bulk_ase)
    else:
        bulk_ase=read('orig_cif_data/'+element+'.cif')
        bulk_pym=AseAtomsAdaptor.get_structure(bulk_ase)
    slabgenall=generate_all_slabs(bulk_pym,max_ind,layers,vacuum_layer,
                                lll_reduce=True,center_slab=True,
                                symmetrize=True,in_unit_planes=True)
    print('Miller Index'+'\t'+'Num of Different Shift(s)')
    slab_M=[]
    for slab in slabgenall:
        slab_M.append([slab.miller_index])
    slab_M_unique = Counter(chain(*slab_M))
    for key in list(slab_M_unique.keys()):
        print(str(key)+'\t'+str(slab_M_unique[key]))

def surf_creator(element,ind,layers,vacuum_layer,option='pymatgen',max_ind=1,unit=True,order=0,save=False,plot=True,opt=True):
    if opt == True:
        bulk_ase=connect('final_database/bulk.db').get_atoms(name=element)
        bulk_pym=AseAtomsAdaptor.get_structure(bulk_ase)
    else:
        bulk_pym=read('orig_cif_data/'+element+'.cif')
    if option=='pymatgen':
        slabgen = SlabGenerator(bulk_pym, ind, layers, vacuum_layer,
                            center_slab=True,lll_reduce=True,in_unit_planes=unit)
        #slabs=slabgen.get_slabs()
        #slabs_symmetric=[slab for slab in slabs if slab.is_symmetric()]
        slabs_symmetric=slabgen.get_slabs(symmetrize=True)
        if len(slabs_symmetric) == 0:
            print('No symmetric slab found!')
        else:
            print('No.'+'\t'+'Layers'+'\t'+'Angles'+'\t\t\t\tCell Length')
            if plot:
                fig=plt.figure(figsize=(8,8))
            layers_ls=[]
            for n,slab in enumerate(slabs_symmetric):
                #temp save for analysis
                os.makedirs(element+'/raw_surf',exist_ok=True)
                surf_location=element+'/raw_surf/'+str(ind)+'_temp'+'.cif'
                CifWriter(slab).write_file(surf_location)
                slab_ase=read(surf_location)
                #slab_ase=AseAtomsAdaptor.get_atoms(slab)
                angles=np.round(slab_ase.get_cell_lengths_and_angles()[3:],decimals=4)
                cell_length=np.round(slab_ase.get_cell_lengths_and_angles()[:3],decimals=4)
                layers=len(np.unique(np.round(slab_ase.positions[:,2],decimals=4)))
                print(str(n)+'\t'+str(layers)+'\t'+str(angles)+'\t'+str(cell_length))
                layers_ls.append(layers)
                if plot:
                    ax=fig.add_subplot(np.ceil(len(slabs_symmetric)/2),2,n+1)
                    plot_slab(slab,ax,adsorption_sites=False,decay=0.25,window=1)
                    ax.set_title('{}: No. {}'.format(slab.miller_index,n),{'fontsize':20})
                    ax.set_xticks([])
                    ax.set_yticks([])
            if os.path.isfile(surf_location):
                os.remove(surf_location)
        if save:
            slab_to_save=slabs_symmetric[order]
            surf_saver(element,slab_to_save,ind,layers_ls[order],order)
    elif option=='pymatgen_all':
        max_ind=max(ind)
        slabgenall=generate_all_slabs(bulk_pym,max_ind,layers,vacuum_layer,
                                lll_reduce=True,center_slab=True,
                                symmetrize=True,in_unit_planes=True)
        slab_RM=[]
        for slab in slabgenall:
            if slab.miller_index == ind:
                slab_RM.append(slab)
        print('No.'+'\t'+'Layers'+'\t'+'Angles'+'\t\t\t\tCell Length')
        if plot:
            fig=plt.figure(figsize=(8,8))
        layers_ls=[]
        surf_location='.'
        for n,slab in enumerate(slab_RM):
            #temp save for analysis
            os.makedirs(element+'/raw_surf',exist_ok=True)
            surf_location=element+'/raw_surf/'+str(ind)+'_temp'+'.cif'
            CifWriter(slab).write_file(surf_location)
            slab_ase=read(surf_location)
            #slab_ase=AseAtomsAdaptor.get_atoms(slab)
            angles=np.round(slab_ase.get_cell_lengths_and_angles()[3:],decimals=4)
            cell_length=np.round(slab_ase.get_cell_lengths_and_angles()[:3],decimals=4)
            layers=len(np.unique(np.round(slab_ase.positions[:,2],decimals=4)))
            print(str(n)+'\t'+str(layers)+'\t'+str(angles)+'\t'+str(cell_length))
            layers_ls.append(layers)
            if plot:
                ax=fig.add_subplot(np.ceil(len(slab_RM)/2),2,n+1)
                plot_slab(slab,ax,adsorption_sites=False,decay=0.25,window=1)
                ax.set_title('{}: No. {}'.format(slab.miller_index,n),{'fontsize':20})
                ax.set_xticks([])
                ax.set_yticks([])
        if os.path.isfile(surf_location):
            os.remove(surf_location)
        if save:
            slab_to_save=slab_RM[order]
            surf_saver(element,slab_RM[order],ind,layers_ls[order],order)
    elif option=='ase':
        slab_ase=surface(bulk_ase,ind,layers=layers,vacuum=vacuum_layer)
        print('No.'+'\t'+'Layers'+'\t'+'Angles'+'\t\t\t\tCell Length')
        angles=np.round(slab_ase.get_cell_lengths_and_angles()[3:],decimals=4)
        cell_length=np.round(slab_ase.get_cell_lengths_and_angles()[:3],decimals=4)
        print(str(0)+'\t'+str(len(np.unique(np.round(slab_ase.positions[:,2],decimals=4))))+'\t'+str(angles)+'\t'+str(cell_length))
        if plot:
            fig=plt.figure(figsize=(8,8))
            ax=fig.add_subplot(111)
            plot_atoms(slab_ase,ax=ax)
            ax.set_title('ASE created: {}'.format(str(ind)),{'fontsize':20})
            ax.set_xticks([])
            ax.set_yticks([])
        if save:
            slab_struc=AseAtomsAdaptor.get_structure(slab_ase)
            layers=len(np.unique(np.round(slab_ase.positions[:,2],decimals=4)))
            surf_saver(element,slab_struc,ind,layers)

def surf_saver(element,slab_to_save,ind,layers,order=0):
    rep_location=element+'/raw_surf'
    if os.path.isdir(rep_location):
        print('WARNING: '+rep_location+' already exists!')
    os.makedirs(rep_location,exist_ok=True)
    surf_location=element+'/raw_surf/'+str(ind)+'_'+str(layers)+'_'+str(order)+'.cif'
    if os.path.isfile(surf_location):
        print('WARNING: '+surf_location+' already exists!')
        print('Raw surface saving fail!')
    else:
        CifWriter(slab_to_save).write_file(surf_location)
        print('Raw surface saving complete!')
