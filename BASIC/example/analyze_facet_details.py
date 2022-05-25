from BASIC.utils import generate_facet_details
import sys
element=f'{sys.argv[1]}_mp-{sys.argv[2]}'
ind=sys.argv[3]
layers=int(sys.argv[4])
generate_facet_details(element,ind,layers,calculated=False)