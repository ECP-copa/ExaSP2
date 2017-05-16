/// \file
/// Communicate halo data such as "ghost" rows with neighboring tasks.

#ifdef DATAEX_HALO

#include "haloExchange.h"

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "matrixMath.h"
#include "decomposition.h"
#include "parallel.h"
#include "performance.h"
#include "constants.h"
#include "mytype.h"

/// A structure to package data for a single row to pack into a
/// send/recv buffer.
typedef struct NonZeroMsgSt
{
   int irow;
   int icol;
   real_t val;
}
NonZeroMsg;

/// \details
DataExchange* initDataExchange(struct DomainSt* domain)
{
   DataExchange* hh = (DataExchange*)malloc(sizeof(DataExchange));
   
   hh->maxExchange = domain->totalProcs;

   hh->exchangeCount = 0;
   hh->exchangeProc = (int*)malloc(hh->maxExchange*sizeof(int));

   hh->bufferSize = domain->localRowExtent * domain->totalCols * (
                    2 * sizeof(int) + sizeof(real_t)); // row, col, value

   if (printRank() && debug == 1)
     printf("bufferSize = %d\n", hh->bufferSize);

   hh->sendBuf = (char*)malloc(hh->bufferSize*sizeof(char));
   hh->recvBuf = (char**)malloc(getNRanks()*sizeof(char*));
   for (int i = 0; i < getNRanks(); i++)
   {
     hh->recvBuf[i] = (char*)malloc(hh->bufferSize*sizeof(char));
   }

   return hh;
}

/// \details
void destroyDataExchange(struct DataExchangeSt* dataExchange)
{
  free(dataExchange->exchangeProc);
  free(dataExchange->sendBuf);

  for (int i = 0; i < getNRanks(); i++)
  {
    free(dataExchange->recvBuf[i]);
  }
  free(dataExchange->recvBuf);

  free(dataExchange);
}

/// Setup for data exchange - post non-blocking reads
void exchangeSetup(struct DataExchangeSt* dataExchange, struct DataMatrixSt* spmatrix, struct DomainSt* domain)
{
  // Update halo processors for matrix
  updateData(dataExchange, spmatrix, domain);

  // Post receives from halo processors and
  // Send local row to halo processors
  if (dataExchange->exchangeCount > 0)
  {
    // Post non-blocking receives
    dataExchange->rlist = (int*)malloc(dataExchange->exchangeCount*sizeof(int));
    for (int i = 0; i < dataExchange->exchangeCount; i++)
    {
      dataExchange->rlist[i] = irecvAnyParallel(dataExchange->recvBuf[i], dataExchange->bufferSize);
    }
  }
}    

/// This is the function that does the heavy lifting for the
/// communication of halo data.
void exchangeData(struct DataExchangeSt* dataExchange, struct DataMatrixSt* spmatrix, struct DomainSt* domain)
{
  if (dataExchange->exchangeCount > 0)
  {
    // Send local rows to each halo processor
    int nSendLen = loadBuffer(dataExchange->sendBuf, spmatrix, domain);
    for (int i = 0; i < dataExchange->exchangeCount; i++)
    {
      int nSend = sendParallel(dataExchange->sendBuf, nSendLen, dataExchange->exchangeProc[i]);
      collectCounter(sendCounter, nSendLen);
    }

    // Receive remote rows from each halo processor
    for (int i = 0; i < dataExchange->exchangeCount; i++)
    {
      int nRecv = waitIrecv(dataExchange->rlist[i]);
      unloadBuffer(dataExchange->recvBuf[i], nRecv, spmatrix, domain); 
      collectCounter(recvCounter, nRecv);
    }
  }

}

/// \details
/// Determine processors in halo
void updateData(struct DataExchangeSt* dataExchange, struct DataMatrixSt* spmatrix, struct DomainSt* domain)
{
  dataExchange->exchangeCount = 0;

  for (int i = domain->localRowMin; i < domain->localRowMax; i++)
  {
    for (int j = 0; j < spmatrix->iia[i]; j++)
    {
      int rnum = spmatrix->jja[i][j];
      if (rnum < domain->localRowMin ||
          rnum >= domain->localRowMax)
      {
        int rowProc = processorNum(domain, rnum);
        addExchangeProc(dataExchange, rowProc);
      }
    }
  }
}

/// \details
/// Gather sparse matrix data to processor 0
void gatherData(struct DataExchangeSt* dataExchange, struct DataMatrixSt* dataMatrix, struct DomainSt* domain)
{
  int myRank = getMyRank();

  // If rank 0, read all blocks that have not been received as halos
  if (myRank == 0)
  {
    free(dataExchange->rlist);
    dataExchange->rlist = (int*)malloc((getNRanks() - dataExchange->exchangeCount)*sizeof(int));
    int ir = 0;
    for (int i = 1; i < getNRanks(); i++)
    {
      if (!isExchangeProc(dataExchange, i))
      {
        dataExchange->rlist[ir] = irecvAnyParallel(dataExchange->recvBuf[ir], dataExchange->bufferSize);
        ir++;
      }
    }

    for (int i = 0; i < ir; i++)
    {
      int nRecv = waitIrecv(dataExchange->rlist[i]);
      unloadBuffer(dataExchange->recvBuf[i], nRecv, dataMatrix, domain);
      collectCounter(recvCounter, nRecv);
    }  
  }

  // Else send block if wasn't sent as halo
  else 
  {
    if (!isExchangeProc(dataExchange, 0))
    {
      int nSendLen = loadBuffer(dataExchange->sendBuf, dataMatrix, domain);
      int nSend = sendParallel(dataExchange->sendBuf, nSendLen, 0);
      collectCounter(sendCounter, nSendLen);
    }
  }
}

/// \details
/// Gather sparse matrix data to processor 0
void allGatherData(struct DataExchangeSt* dataExchange, struct DataMatrixSt* dataMatrix, struct DomainSt* domain)
{
  int myRank = getMyRank();

  // Post reads for all blocks that have not been received as halos
  free(dataExchange->rlist);
  dataExchange->rlist = (int*)malloc((getNRanks() - dataExchange->exchangeCount)*sizeof(int));
  int ir = 0;
  for (int i = 0; i < getNRanks(); i++)
  {
    if (i != myRank && !isExchangeProc(dataExchange, i))
    {
      dataExchange->rlist[ir] = irecvAnyParallel(dataExchange->recvBuf[ir], dataExchange->bufferSize);
      ir++;
    }
  }

  // Send to other ranks that aren't halos
  int nSendLen = loadBuffer(dataExchange->sendBuf, dataMatrix, domain);
  for (int i = 0; i < getNRanks(); i++)
  {
    if ( i != myRank && !isExchangeProc(dataExchange, i))
    {
      int nSend = sendParallel(dataExchange->sendBuf, nSendLen, i);
      collectCounter(sendCounter, nSendLen);
    }
  }

  for (int i = 0; i < ir; i++)
  {
    int nRecv = waitIrecv(dataExchange->rlist[i]);
    unloadBuffer(dataExchange->recvBuf[i], nRecv, dataMatrix, domain);
    collectCounter(recvCounter, nRecv);
  }
}

/// \details
/// Check if proc is a halo proc
int isExchangeProc(struct DataExchangeSt* dataExchange, int rproc)
{
  for (int i = 0; i < dataExchange->exchangeCount; i++)
  {
    if (dataExchange->exchangeProc[i] == rproc) return 1;
  }
  return 0;
}

/// \details
/// Add halo processor number to list if not currently present
void addExchangeProc(struct DataExchangeSt* dataExchange, int rproc)
{
  if (dataExchange->exchangeCount == 0)
  {
    dataExchange->exchangeProc[0] = rproc;
    dataExchange->exchangeCount = 1;

  }
  else
  {
    if (isExchangeProc(dataExchange, rproc) == 1) return;

    dataExchange->exchangeProc[dataExchange->exchangeCount] = rproc;
    dataExchange->exchangeCount++;
  }
}


/// \details
/// The loadBuffer function for a halo exchange of row data.
int loadBuffer(char* buf, struct DataMatrixSt* xmatrix, struct DomainSt* domain)
{
  NonZeroMsg* rbuf = (NonZeroMsg*) buf;

  int nBuf = 0;
  for (int i = domain->localRowMin; i < domain->localRowMax; i++)
  {
    for (int j = 0; j < xmatrix->iia[i]; j++)
    {
      rbuf[nBuf].irow = i;
      rbuf[nBuf].icol = xmatrix->jja[i][j];
      rbuf[nBuf].val = xmatrix->val[i][j];

      nBuf++;
    }
  }

  return nBuf*sizeof(NonZeroMsg);  
}

/// \details
/// The unloadBuffer function for a halo exchange of row data.
void unloadBuffer(char* buf, int bufSize, struct DataMatrixSt* xmatrix, struct DomainSt* domain)
{
  NonZeroMsg* rbuf = (NonZeroMsg*) buf;

  int nBuf = bufSize / sizeof(NonZeroMsg);
  assert(bufSize % sizeof(NonZeroMsg) == 0);

  // Assume all non-zeros for a row areastored contiguously
  int rcurrent = -1;
  int j = 0;
  for (int i = 0; i < nBuf; i++)
  {
    int irow = rbuf[i].irow;

    if (irow != rcurrent)
    {
      xmatrix->iia[irow] = 0;
      rcurrent = irow;
      j = 0;
    }
    xmatrix->iia[irow]++;
    xmatrix->jja[irow][j] = rbuf[i].icol;
    xmatrix->val[irow][j] = rbuf[i].val;
    j++;

  }

}

#endif
