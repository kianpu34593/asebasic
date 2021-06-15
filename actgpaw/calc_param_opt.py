import numpy as np
class calc_param:
    def __init__(self,gpaw_calc,element,solver_step=0.05,solver_fmax=0.01,init_magmom):
        self.calc_dict = gpaw_calc
        self.atoms = bulk_builder(element)
        self.atoms = self.atoms.set_initial_magnetic_moments(init_magmom*np.ones(len(atoms)))
        self.solver_step=solver_step
        self.solver_fmax=solver_fmax

    def h_size(self):
        atoms=self.atoms.set_calculator(self.calc_dict)
        opt.optimize_bulk(atoms,step=self.solver_step,fmax=self.solver_fmax,
                            location=)
        return atoms
    def k_points(self):
        atoms=self.atoms.set_calculator(self.calc_dict)
        opt.optimize_bulk(atoms,step=self.solver_step,fmax=self.solver_fmax,
                            location=)
        return atoms
    def smearing_width(self,option='guess'):
        if option = 'guess':
            try:
                opt.optimize_bulk(atoms,step=self.solver_step,fmax=self.solver_fmax,
                            location=)
                return atoms
            except:
                print("SCF computation did not converge for the GUESS option.")

        elif option = 'converge':
            opt.optimize_bulk(atoms,step=self.solver_step,fmax=self.solver_fmax,
                            location=)
            return atom
            


def bulk_builder