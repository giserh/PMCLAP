/*
 * Probabilistic Maximal Covering Location Allocation Problem (PMCLAP)
 * Queuing Maximal Covering Location-Allocation Model (QM-CLAM)
 *  Probabilistic, Maximal Covering Location-Allocation Models for Congested Systems
 *   Vladimir Marianov & Daniel Serra
 *    Journal of Regional Science,
 *    Vol. 38, No. 3, 1998, pp. 401-424
 */

#ifndef _PMCLAP_H
#define _PMCLAP_H

#include "point.h"
#include <ilcplex/ilocplex.h>

ILOSTLBEGIN

typedef IloArray<IloBoolVarArray> BoolVarMatrix;
typedef IloArray<BoolVarMatrix> BoolVarArrayMatrix;
typedef IloArray<IloNumArray> NumMatrix;
typedef IloArray<IloIntArray> IntMatrix;

#define QUEUE_SIZE 0
#define WAITING_TIME 1

void QM_CLAM
(instance*, /* Instance */
 IloInt,    /* Number of facilities */
 IloBool,   /* Constraint type queue_size 0 or waitin_time 1 */
 IloNum,    /* Congestion parameter */
 IloNum,    /* Minimum probability */
 IloNum,    /* Service parameter*/
 IloNum     /* Portion of demand to rate parameter */
);

IloNum RHS_Queue_Size(IloNum,IloNum,IloNum);
IloNum RHS_Waiting_time(IloNum,IloNum,IloNum);
void usage();

#endif

/* eof */
