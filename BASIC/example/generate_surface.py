from BASIC.utils import generate_all_slab
import sys

element=f'{sys.argv[1]}_mp-{sys.argv[2]}'
slab_facet_dict=generate_all_slab(element,
                calculated=False,
                max_ind=1)
facet_ls=[]
shift_ls=[]
order_ls=[]
for facet,shift_temp_ls in slab_facet_dict.items():
    order_temp_ls=[]
    order_count=0
    for shift in shift_temp_ls:
        facet_ls.append(facet)
        shift_ls.append(shift)
        order_temp_ls.append(order_count)
        order_count+=1
    order_ls+=order_temp_ls

print(facet_ls,shift_ls,order_ls)