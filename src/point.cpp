
#include "point.h"


float dist(point *i,point *j){
  return sqrt((i->x - j->x)*(i->x - j->x) + (i->y - j->y) * (i->y - j->y));
}

instance* read_points(const char *filename){
  int i;
  instance *I;
  point *puntos;
  ifstream file(filename);
  
  if (!file) {
    cerr << "ERROR: colud not open file '" << filename
	 << "' for reading" << endl;
    return NULL;
  }

  I = new instance;
  file >> I->n >> I->S;
  puntos = new point[I->n];
  for (i = 0;i < I->n;i++)
    file >> puntos[i].x >> puntos[i].y >> puntos[i].demand;
  I->points = puntos;
  file.close();
  
  return I;
}
