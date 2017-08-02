data: This contains example matrices in Matrix Market Format.
These matrices serve as input Hamiltonians to ExaSP2.
They represent different the following.
 1) polyethylene chains of 512 and 1024 molecules,
 2) nitro-methane 100 molecules, from last SCF iteration.

poly_chain.512.mtx - 512 molecules, 3072 atoms, 6144 orbitals,
number of non-zeroes = 98304, N = 6144, M = 260.

poly_chain.1024.mtx - 1024 molecules, 6144 atoms, 12,288 orbitals,
number of non-zeroes = 196608, N = 12288, M = 260.

nm.100.mtx - 100 molecules, 700 atoms, 1900 orbitals
number of non-zeroes = 730866, N = 1900, M = 1900

Note the M values are suggested for SP2 Basic. The other variants may
require larger values.
