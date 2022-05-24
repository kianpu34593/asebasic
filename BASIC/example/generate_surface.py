from BASIC.utils import generate_all_slab
import sys

element=f'{sys.argv[1]}_mp-{sys.argv[2]}'
slab_facet_dict=generate_all_slab(element,
                calculated=False,
                max_ind=1)
print(slab_facet_dict)