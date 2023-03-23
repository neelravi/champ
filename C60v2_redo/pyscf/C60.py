#!/usr/bin/env python

import numpy as np
import sys
import os 
#sys.path.append('/p/home/jusers/landinez1/juwels/Codes/pyscf/') # if local installation


import pyscf
from pyscf import gto, scf, dft, cc
from pyscf.tools import c60struct


os.environ['OMP_NUM_THREADS'] = "48"


#basis='cc-pVTZ'
basisset= {'C': gto.basis.parse("""
#C.cc-pVDZ.nwchem
C s
13.073594 0.0051583
6.541187 0.0603424
4.573411 -0.1978471
1.637494 -0.0810340
0.819297 0.2321726
0.409924 0.2914643
0.231300 0.4336405
0.102619 0.2131940
0.051344 0.0049848
C s
0.127852 1.000000
C p
9.934169 0.0209076
3.886955 0.0572698
1.871016 0.1122682
0.935757 0.2130082
0.468003 0.2835815
0.239473 0.3011207
0.117063 0.2016934
0.058547 0.0453575
0.029281 0.0029775
C p
0.149161 1.000000
C d
0.561160 1.000000
""")}

#passing coordinates to Borh
coords=c60struct.make60(1.46,1.38)*1.8897259886

mol = pyscf.M(atom=[('C', r) for r in coords],
              basis=basisset,
              max_memory=40000)

mol.ecp={'C': gto.basis.parse_ecp(
"""
#C.ccECP.nwchem
C nelec 2
C ul
1 14.43502 4.00000
3 8.39889 57.74008
2 7.38188 -25.81955
C S
2 7.76079 52.13345
""")}
mol.unit     = 'B'
mol.symmetry = False
mol.verbose = 5
#cell.build(cart=True)
mol.build(cart=False)
print(mol.atom_coords())

mf = scf.fast_newton(mol.RHF())
print('Initial E(tot) %.15g ' % mf.e_tot)
#mf = scf.RHF(mol)
scfdump= "RHF.dump"               
mf.chkfile = scfdump;
ehf = mf.kernel()
print("HF energy  = %.17g" % ehf)


print("\n Analysis by PySCF \n")

mf.analyze()


##print(dir(mf))


print("\n mo's coeff \n")
print(mf.mo_coeff)


dm = mf.make_rdm1()

#build h1e first
#pure kinetic part
t = mol.intor_symmetric('int1e_kin')
print("t shape", t.shape)
t = np.einsum('ij,ji', t, dm).real
print("Kinectic energy only", t)

#pseudopotential part/ pen
pp=0.0
pen1=0.0

if mol._pseudo:
    # Although mol._pseudo for GTH PP is only available in Cell, GTH PP
    # may exist if mol is converted from cell object.
    from pyscf.gto import pp_int
    pen1 = pp_int.get_gth_pp(mol)
else:
    pen1 = mol.intor_symmetric('int1e_nuc')
    
if len(mol._ecpbas) > 0:
    pp = mol.intor_symmetric('ECPscalar')

pen=pp+pen1

print("pp",pp)
print("pen1",pen1)
print("pen",pen)
    
pen = np.einsum('ij,ji', pen, dm).real
print("pen (energy)",pen)

h1e_t=t+pen
print('h1e_t',h1e_t)


# pure pyscf function for testing
h1e = mf.get_hcore()
h1e = np.einsum('ij,ji->', h1e, dm).real
print('h1e (energy)',h1e)

##Returns matrix Vhf = 2*J - K.  Vhf can be a list matrices, corresponding to the rdm
#pyscf function
vhf = mf.get_veff(mol, dm)
ee_coul = np.einsum('ij,ji->', vhf, dm).real * .5
print('ee_coul',ee_coul)

te_t=h1e_t+ee_coul
te=h1e+ee_coul

print("Total electron Energy test", te_t)
print("Total electron Energy pyscf", te)


#nuclear repulsion energy
p_nn = mf.energy_nuc()
print("p_nn",p_nn)


print("Total Energy test", te_t+p_nn)
print("Total Energy split pyscf", te+p_nn)
print("Total Energy pyscf", ehf)





#mycc = cc.RCCSD(mf)
#scfdump= "RCCSD.dump"               
#mycc.chkfile = scfdump;
#mycc.kernel()
#print("RCCSD energy =", mycc.e_tot)
