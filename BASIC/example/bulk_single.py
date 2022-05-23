from BASIC.compute import bulk_compute
from gpaw import GPAW

calc=GPAW(xc = 'PBE')
converge_tuple=('single_compute','')

bulk=bulk_compute(element='Li_mp-135',
            calculator_setting=calc,
            converge_parameter=converge_tuple)