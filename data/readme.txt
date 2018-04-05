data: This contains example matrices in Matrix Market Format.
These matrices serve as input Hamiltonians to ExaSP2.
They represent different the following.
 1) polyethylene chains of 512 and 1024 molecules.
 2) nitro-methane 100 molecules, from last SCF iteration.
 3) nitro-methane 200 molecules, from last SCF iteration.
 4) trpcage protein solvated in water, 8352 atoms.
 5) alanine 19945 atoms.
 6) glutamic acid 210 atoms.

poly_chain.512.mtx - 512 molecules, 3072 atoms, 6144 orbitals,
number of non-zeroes = 98304, N = 6144, M = 260.

poly_chain.1024.mtx - 1024 molecules, 6144 atoms, 12,288 orbitals,
number of non-zeroes = 196608, N = 12288, M = 260.

nm.100.mtx - 100 molecules, 700 atoms, 1900 orbitals
number of non-zeroes = 730866, N = 1900, M = 1900

nm.200.mtx - 200 molecules, 1400 atoms, 3800 orbitals
number of non-zeroes = 1442828, N = 3800, M = 3800

trpcage_8K.mtx - 8352 atoms, 16863 orbitals
number of non-zeroes = 116690, N = 16863, M = 1500, BNDFIL = 0.661626 (OCC = 11157) 

alanine.20k.mtx - 19945 atoms, 41185 orbitals
number of non-zeroes = 374997, N = 41185, M = 3000, BNDFIL = 0.6498968 (OCC = 26766)

F_ort_matrix_0008.mtx - 210 atoms, 1204 orbitals
number of non-zeroes = 1388944, N = 1204, M = 1204, BNDFIL = 0.35299 (OCC = 425)

Note the M values are suggested for SP2 Basic. The other variants may
require larger values.
