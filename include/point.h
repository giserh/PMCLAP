
#ifndef _POINT_H_
#define _POINT_H_

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

struct point {
  float x,y;
  /* (x,y) coordinates of the point */
  float demand;
  /* demand of the point */
  bool *N_i;
  /* array that indicates whether the point j is in range S */
};

struct instance {
  int n;
  /* numbers of nodes */
  float S;
  /* maximum distance */
  struct point *points;
  /* set of nodes */
};

typedef struct point point;
typedef struct instance instance;

float dist(point*,point*);
instance* read_points(const char*);

#endif

/* eof */
