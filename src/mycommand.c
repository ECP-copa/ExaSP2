/// \file
/// Handle command line arguments.

#include "mycommand.h"

#include <string.h>
#include <stdlib.h>

#include "cmdLineParser.h"
#include "parallel.h"
#include "mytype.h"

/// \page pg_running_cosp2 Running ExaSp2
///
/// \section sec_command_line_options Command Line Options
///
/// ExaSp2 accepts a number of command line options to set the parameters
/// of the simulation. Every option has both a long form and a short
/// form. The long and short form of the arguments are entirely
/// interchangeable and may be mixed. 
///
/// Supported options are:
///
/// | Long  Form    | Short Form  | Default Value | Description
/// | :------------ | :---------: | :-----------: | :----------
/// | \--help       | -h          | N/A           | print this message
/// | \--hmatName   | -f          |               | H matrix file name in Matrix Market format
/// | \--N          | -n          | 1600          | number of rows
/// | \--M          | -m          | 1600 or N     | max non-zeroes per row
/// | \--mtype      | -y          | 2 (ellpack)   | matrix type
/// | \--minIter    | -w          | 25            | min sp2 iters
/// | \--maxIter    | -x          | 100           | max sp2 iters
/// | \--nocc       | -c          | 0.0           | number of occupied states
/// | \--bndfil     | -b          | 0.5           | bndfil
/// | \--beta       | -k          | 0.0           | beta=1/KBT
/// | \--eps        | -e          | 1.0E-05       | threshold for sparse math
/// | \--idemtol    | -i          | 1.0E-14       | threshold for SP2 loop
/// | \--gen        | -g          | 0             | generate H matrix if 1
/// | \--dout       | -o          | 0             | write out density matrix if 1
/// | \--dbg        | -d          | 0             | write debug messages if 1
///
/// Notes: 
/// 
/// 
/// \subsection ssec_example_command_lines Examples
///
/// All of the examples below assume:
/// - The ExaSp2 bin directory is located in ../bin
///
/// Running in the examples directory will satisfy these requirements.
///
/// ------------------------------
///
/// The canonical base simulation, is 
///
///     $ mpirun -np 1 ../bin/ExaSp2Parallel 
/// 
/// Or, if the code was built without MPI:
///
///     $ ../bin/ExaSp2
///
/// ------------------------------
///
/// \subsubsection cmd_examples_hmatrix Specify H matrix
///
/// To run with a given H matrix(in Matrix Market format), specify -f:
///
///     $ ../bin/ExaSp2-parallel -f my_hmatrix.mtx 
///
/// ------------------------------
///
/// To run with a chosen size a matrix generated:
///
///     $ ../bin/ExaSp2-Parallel -N 10000 -M 124
///
/// ------------------------------
///
/// To run with a chosen SP2 loop tolerance:
///
///     $ ../bin/ExaSp2-Parallel -N 20000 -M 40 -i 1.0E-08 
/// 
///
/// ------------------------------
///

/// \details Initialize a Command structure with default values, then
/// parse any command line arguments that were supplied to overwrite
/// defaults.
///
/// \param [in] argc the number of command line arguments
/// \param [in] argv the command line arguments array
Command parseCommandLine(int argc, 
                         char** argv)
{
   Command cmd;

   memset(cmd.hmatName, 0, 1024);
   cmd.N = 1600;
   cmd.M = 1600;
   cmd.mtype = 2;
   cmd.dout = 0;
   cmd.gen = 0;
   cmd.minsp2iter = 25;
   cmd.maxsp2iter = 100;
   cmd.nsteps = 18;
   cmd.osteps = 0;
   cmd.debug = 0;
   cmd.nocc = 0.0;
   cmd.eps = 1.0E-05;
   cmd.idemTol = 1.0E-14;
   cmd.bndfil = 0.5;
   cmd.beta = 0.0;
   cmd.tscale = 1.0;
   cmd.occLimit = 1.0E-09;
   cmd.traceLimit = 1.0E-12;

   int help=0;
   // add arguments for processing.  Please update the html documentation too!
   addArg("help",       'h', 0, 'i',  &(help),             0,             "print this message");
   addArg("hmatName",   'f', 1, 's',  cmd.hmatName,   sizeof(cmd.hmatName), "H matrix file name");
   addArg("N",          'n', 1, 'i',  &(cmd.N),            0,             "rows in matrix");
   addArg("M",          'm', 1, 'i',  &(cmd.M),            0,             "non-zeroes per row");
   addArg("mtype",      'y', 1, 'i',  &(cmd.mtype),        0,             "matrix type (1-dense,2-ellpack)");
   addArg("minIter",    'w', 1, 'i',  &(cmd.minsp2iter),   0,             "min sp2 iters");
   addArg("maxIter",    'x', 1, 'i',  &(cmd.maxsp2iter),   0,             "max sp2 iters");
   addArg("nsteps",     's', 1, 'i',  &(cmd.nsteps),       0,             "num sp2 iters");
   addArg("occSteps",   'c', 1, 'i',  &(cmd.osteps),       0,             "num occ iters");
   addArg("gen",        'g', 1, 'i',  &(cmd.gen),          0,             "generate H matrix");
   addArg("dout",       'w', 1, 'i',  &(cmd.dout),         0,             "write out density matrix");
   addArg("debug",      'd', 1, 'i',  &(cmd.debug),        0,             "write out debug messages");
   addArg("nocc",       'o', 1, 'd',  &(cmd.nocc),         0,             "number of occupied states");
   addArg("bndfil",     'b', 1, 'd',  &(cmd.bndfil),       0,             "bndfil");
   addArg("beta",       'k', 1, 'd',  &(cmd.beta),         0,             "beta=1/KBT");
   addArg("eps",        'e', 1, 'd',  &(cmd.eps),          0,             "threshold for sparse math");
   addArg("idemtol",    'i', 1, 'd',  &(cmd.idemTol),      0,             "threshold for SP2 loop");
   addArg("tscale",     't', 1, 'd',  &(cmd.tscale),       0,             "scaling factor");
   addArg("traceLimit", 'a', 1, 'd',  &(cmd.traceLimit),   0,             "trace limit");
   addArg("occLimit",   'r', 1, 'd',  &(cmd.occLimit),     0,             "occ err limit");

   processArgs(argc,argv);

   // If user didn't set hmatName, set to generate H matrix.
   if (strlen(cmd.hmatName) == 0) 
   {
     cmd.gen = 1;
   }
      
   if (help)
   {
      printArgs();
      freeArgs();
#ifdef DO_MPI
      destroyParallel();
#endif
      exit(0);
   }
   freeArgs();

   return cmd;
}
