/// \file
/// Handle command line arguments.

#ifndef MYCOMMAND_H
#define MYCOMMAND_H

#include <stdio.h>

#include "mytype.h"

/// A structure to hold the value of every run-time parameter that can
/// be read from the command line.
typedef struct CommandSt
{
   char hmatName[1024]; //!< name of the dense H matrix file
   int N;               //!< number of rows in H matrix (N x N)
   int M;               //!< max number of non-zeroes in H matrix row
   int mtype;           //!< matrix type (1-dense, 2-ellpack, ...)
   int dout;            //!< if == 1, write out density matrix
   int gen;             //!< if == 1, generate sparse hamiltonian
   int minsp2iter;      //!< minimum number of sp2 iterations
   int maxsp2iter;      //!< maximum number of sp2 iterations
   int debug;           //!< if == 1, write out debug messages
   int nsteps;          //!< number of SP2 steps
   int osteps;          //!< number of occupation loop steps

   real_t nocc;         //!< number of occupied states
   real_t eps;          //!< threshold for sparse math 
   real_t idemTol;      //!< threshold for SP2 loop
   real_t bndfil;       //!< band fill
   real_t tscale;       //!< scaling factor for SP2_Fermi
   real_t mu;           //!< the chemical potential 
   real_t beta;         //!< 1/KBT
   real_t occLimit;     //!< occupation error limit
   real_t traceLimit;   //!< trace comparison limit
} Command;

/// Process command line arguments into an easy to handle structure.
Command parseCommandLine(int argc, 
                         char** argv);

#endif
