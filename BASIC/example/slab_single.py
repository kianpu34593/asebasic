from compute import slab_compute
from gpaw import GPAW

calc=GPAW(xc = 'PBE')
converge_tuple=('single_compute','')
computation_setting={'miller_plane':100,
                    'shift':0,
                    'order':0,
                    'fix_layer':2,
                    'fix_mode':'bottom',
                    }

slab=slab_compute(element='Li_mp-135',
            calculator_setting=calc,
            converge_parameter=converge_tuple,
            computation_setting=computation_setting,
            restart_calculation=False,
            compute_surface_energy=True,
            )