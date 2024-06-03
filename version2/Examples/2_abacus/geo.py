from ase.calculators.FCPelectrochem import FCP
from ase.calculators.abacus import Abacus, AbacusProfile
from ase.io import read
from ase.optimize import LBFGS

# abacus parameters
pseudo_dir = r'/home/liuyu/data/ABACUS/PP_ORB'
basis_dir = r'/home/liuyu/data/ABACUS/PP_ORB'
pp = {'C':'C_ONCV_PBE-1.0.upf',
      'N':'N_ONCV_PBE-1.0.upf',
      'Fe':'Fe_ONCV_PBE-1.0.upf',
      'O':'O_ONCV_PBE-1.0.upf'}
basis = {'C':'C_gga_8au_60Ry_2s2p1d.orb',
         'N':'N_gga_8au_60Ry_2s2p1d.orb',
         'Fe':'Fe_gga_9au_60Ry_4s2p2d1f.orb',
         'O':'O_gga_7au_60Ry_2s2p1d.orb'}

kpts = {'size': [3, 3, 1], 'gamma': True}
parameters = {
   'xc': 'pbe',
   'ecutwfc': 60,
   'nspin': 2,
   'smearing_method': 'gaussian',
   'smearing_sigma': 0.01,
   'basis_type': 'lcao',
   'ks_solver': 'genelpa',
   'mixing_type': 'pulay',
   'scf_thr': 1e-5,
   'scf_nmax': 500,
   'calculation': 'scf',
   'cal_force' : 1,
   'cal_stress' : 1
}

profile = AbacusProfile(argv=['mpirun', '-n', '24', 'abacus'])
cal_abacus = Abacus(profile=profile, pp=pp, pseudo_dir=pseudo_dir, basis=basis, basis_dir=basis_dir, kpts=kpts, **parameters)

# FCP parameters
cal_FCP=FCP(innercalc = cal_abacus,
            fcptxt = 'log-fcp.txt',
            U = 0.8,
            NELECT = 216.5,
            C = 1/80,    #1/k  capacitance per A^2
            FCPmethod = 'Newton-fitting',
            FCPconv = 0.01,
            NELECT0 = 218, 
            adaptive_lr = False,
            work_ref = 4.6,
            max_FCP_iter = 10000
            )

#________________________________________________________________

atoms=read(filename='STRU', format='abacus')
atoms.calc = cal_FCP
dyn=LBFGS(atoms, trajectory='fp.traj')
dyn.run(fmax=0.01)
