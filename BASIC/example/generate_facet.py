from BASIC.utils import surface_creator
import sys

element=f'{sys.argv[1]}_mp-{sys.argv[2]}'
facet=(sys.argv[3], float(sys.argv[4]), int(sys.argv[5]))
surface_creator(element,
                facet=facet,
                calculated=False,
                layers=int(sys.argv[6]),
                save_to_disk=True)