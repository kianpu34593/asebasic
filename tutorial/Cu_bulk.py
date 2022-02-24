from gpaw import GPAW,Mixer,Davidson
from ase.calculators.calculator import kptdensity2monkhorstpack as kdens2mp
from actgpaw import bulk_autoconv as bulk_ac

# input the material
element = 'Cu_mp-30' #the name should be the same as the cif file
element_atom = bulk_ac.bulk_builder(element) #cif --> ase atom

# initalize the calulator
## convert k density to kpts based on cell size
kpts = kdens2mp(element_atom) 
calc=GPAW(xc = 'PBE',
            h = 0.16,
            kpts = kpts,
            spinpol = False,
            maxiter = 333,
            mixer = Mixer(0.05,5,50),
            eigensolver = Davidson(3),
            occupations = {'name':'fermi-dirac','width':0.1})

#call bulk_auto_conv module
bulk_ac.bulk_auto_conv(element, #input material
                        calc, #initial calculator
                        rela_tol=10*10**(-3), #convergence criteria
                        temp_print=True, #print out the convergence process
                        ) 