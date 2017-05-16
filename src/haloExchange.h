/// \file
/// Communicate halo data such as "ghost" rows with neighboring tasks.

#ifndef __HALO_EXCHANGE_
#define __HALO_EXCHANGE_

#include <assert.h>

#include "matrixMath.h"
#include "decomposition.h"
#include "mytype.h"

typedef struct DataExchangeSt
{
  int maxExchange;     //!< max number of processors for halo
  int bufferSize;      //!< max size in bytes for send/receive buffer
  int exchangeCount;   //!< number of processors to send/reeive

  int* exchangeProc;   //!< array of halo procs to send/receive
  int* rlist;          //!< array of request indeces for non-blocking receives

  char* sendBuf;       //!< send buffer
  char** recvBuf;      //!< recv buffer
} 
DataExchange;

/// Create a DataExchange.
DataExchange* initDataExchange(struct DomainSt* domain);

/// DataExchange destructor.
void destroyDataExchange(struct DataExchangeSt* dataExchange);

/// Update data to be exchanged.
void updateData(struct DataExchangeSt* DataExchange, struct DataMatrixSt* xmatrix, struct DomainSt* domain);

/// Exchange setup - post reads
void exchangeSetup(struct DataExchangeSt* dataExchange, struct DataMatrixSt* xmatrix, struct DomainSt* domain);

/// Execute a data exchange.
void exchangeData(struct DataExchangeSt* dataExchange, struct DataMatrixSt* xmatrix, struct DomainSt* domain);

/// Add exchange proc if not previously added
void addExchangeProc(struct DataExchangeSt* dataExchange, int rproc);
int isExchangeProc(struct DataExchangeSt* dataExchange, int rproc);

/// Buffer functions
int loadBuffer(char* buf, struct DataMatrixSt* xmatrix, struct DomainSt* domain);
void unloadBuffer(char* buf, int bufSize, struct DataMatrixSt* xmatrix, struct DomainSt* domain);

/// Gather functions
void gatherData(struct DataExchangeSt* dataExchange, struct DataMatrixSt* xmatrix, struct DomainSt* domain);
void allGatherData(struct DataExchangeSt* dataExchange, struct DataMatrixSt* xmatrix, struct DomainSt* domain);

#endif
