from actgpaw import ads_selector
from gpaw import GPAW, MixerSum, Mixer, MixerDif, Davidson
from ase.db import connect

# specify the material and miller index of interest
element = "Cu_mp-30"
struc = "111"

# read the optimized conventional cell in the database
element_surf = connect("final_database/surf.db").get(name=element + "(" + struc + ")")
h = element_surf.h
xc = element_surf.xc
sw = element_surf.sw
spin = element_surf.spin
kpts = [int(i) for i in (element_surf.kpts).split(",")]

# set up the calculator
calc = GPAW(
    xc=xc,
    h=h,
    kpts=kpts,
    symmetry={"point_group": False},
    eigensolver=Davidson(3),
    mixer=Mixer(beta=0.05, nmaxold=5, weight=50),
    spinpol=spin,
    maxiter=333,
    occupations={"name": "fermi-dirac", "width": sw},
    poissonsolver={"dipolelayer": "xy"},
)

# call ads_selector module
ads_selector.ads_auto_select(
    element,
    struc,
    calc,
    ads="Li",  # specify adsorbate
    ads_pot_e=-1.89678,  # adsorbate energy
    size="1x1",  # specify the size of the supercell (xy)
    temp_print=True,  # print out the convergence process
)
