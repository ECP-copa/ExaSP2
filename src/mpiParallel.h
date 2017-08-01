/// \file
/// Wrappers for MPI functions.

#ifndef _MPI_PARALLEL_H_
#define _MPI_PARALLEL_H_

#include "mytype.h"

/// Structure for use with MPI_MINLOC and MPI_MAXLOC operations.
typedef struct RankReduceDataSt
{
   double val;
   int rank;
} RankReduceData;

/// Return total number of processors.
int getNRanks(void);

/// Return local rank.
int getMyRank(void);

/// Return non-zero if printing occurs from this rank.
int printRank(void);

/// Print a timestamp and message when all tasks arrive.
void timestampBarrier(const char* msg);

/// Wrapper for MPI_Init.
void initParallel(int *argc, 
                  char ***argv);

/// Wrapper for MPI_Finalize.
void destroyParallel(void);

/// Wrapper for MPI_Barrier(MPI_COMM_WORLD).
void barrierParallel(void);

/// Wrapper for MPI_Sendrecv.
int sendReceiveParallel(const void* sendBuf, 
                        const int sendLen, 
                        const int dest,
                        void* recvBuf, 
                        const int recvLen, 
                        const int source);

/// Wrapper for MPI_Send.
int sendParallel(const void* sendBuf, 
                 const int sendLen, 
                 const int dest);

/// Wrapper for MPI_Isend, non-blocking send.
int isendParallel(const void* sendBuf, 
                  const int sendLen, 
                  const int dest);

/// Wrapper for MPI_Recv from any processor.
int recvAnyParallel(void* recvBuf, 
                    const int recvLen);

/// Wrapper for MPI_Irecv from any processor, non-blocking receive.
int irecvAnyParallel(void* recvBuf, 
                     const int recvLen);

/// Wrapper for MPI_Recv.
int recvParallel(void* recvBuf, 
                 const int recvLen,
                 const int source);

/// Wrapper for MPI_Recv, non-blocking receive.
int irecvParallel(void* recvBuf, 
                  const int recvLen, 
                  const int source);

/// Wrapper for MPI_Wait on non-blocking receive.
int waitIrecv(int rind);

/// Wrapper for MPI_Wait on non-blocking send.
int waitIsend(int rind);

/// Wrapper for MPI_Test on non-blocking receive.
int testIrecv(int rind);

/// Wrapper for MPI_Test on non-blocking send.
int testIsend(int rind);

/// Wrapper for MPI_Allreduce integer sum.
void addIntReduce2(int* value0, 
                   int* value1);
void addIntParallel(const int* sendBuf, 
                    int* recvBuf, 
                    const int count);

/// Wrapper for MPI_Allreduce real sum.
void addRealReduce2(real_t* value0, 
                    real_t* value1);
void addRealParallel(const real_t* sendBuf, 
                     real_t* recvBuf, 
                     const int count);

/// Wrapper for MPI_Allreduce double sum.
void addDoubleParallel(const double* sendBuf, 
                       double* recvBuf, 
                       const int count);

/// Wrapper for MPI_Allreduce real min.
void minRealReduce(real_t* value);
void minRealParallel(const real_t* sendBuf, 
                     real_t* recvBuf, 
                     const int count);

/// Wrapper for MPI_Allreduce integer max.
void maxIntReduce2(int* value0, 
                   int* value1);
void maxIntParallel(const int* sendBuf, 
                    int* recvBuf, 
                    const int count);

/// Wrapper for MPI_Allreduce real max.
void maxRealReduce(real_t* value);
void maxRealParallel(const real_t* sendBuf, 
                     real_t* recvBuf, 
                     const int count);

/// Wrapper for MPI_Allreduce double min with rank.
void minRankDoubleParallel(const RankReduceData* sendBuf, 
                           RankReduceData* recvBuf, 
                           const int count);

/// Wrapper for MPI_Allreduce double max with rank.
void maxRankDoubleParallel(const RankReduceData* sendBuf, 
                           RankReduceData* recvBuf, 
                           const int count);

/// Wrapper for MPI_Bcast
void bcastParallel(const void* buf, 
                   const int len, 
                   const int root);

///  Return non-zero if code was built with MPI active.
int builtWithMpi(void);

#endif

