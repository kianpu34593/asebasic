from ase.db import connect
from gpaw import GPAW, Mixer, MixerDif, Davidson, PoissonSolver
from actgpaw import surf_autoconv as surf_ac

# read the optimized conventional cell in the database
element = "Cu_mp-30"
element_bulk = connect("final_database/bulk.db").get(name=element)
h = element_bulk.h
xc = element_bulk.xc
sw = element_bulk.sw
spin = element_bulk.spin
## all settings but kpts (due to structure dependence)

# miller index of interest
struc = "111"

# set up the calculator
calc = GPAW(
    xc=xc,
    h=h,
    symmetry={"point_group": False},
    eigensolver=Davidson(3),
    mixer=Mixer(beta=0.05, nmaxold=5, weight=50),
    spinpol=spin,
    maxiter=500,
    occupations={"name": "fermi-dirac", "width": sw},
    poissonsolver={'dipolelayer': 'xy'},
)

# call surf_auto_conv module
surf_ac.surf_auto_conv(
    element, #input material
    struc, #miller index of interest
    calc, #calculator
    generator="import", #import the slab model 
    pbc_all=False, #periodic boundary condition true for all directions
    init_layer=4, #initial slab layers
    interval=2, #interval between layers
    fix_layer=2, #number of fixed layers
    fix_option='bottom', #constrain the bottom layers
    vac=10, #vacuum size (Ang)
    rela_tol=5, #convergence criteria (%)
    temp_print=True, #print out the convergence process
)
