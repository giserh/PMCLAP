
#include "PMCLAP.h"

int main(int argc, char **argv){
  IloInt i,j;
  // Datos del problema
  instance *I;

  IloNum mu = 72;
  IloNum f = 0.015;

  IloInt p = 2;
  IloBool constraint_type = WAITING_TIME;
  IloNum congestion_parameter = 0;
  IloNum alpha = 85;
  const char* filename;
  if(argc > 5){
    filename = argv[1];
    p = atoi(argv[2]);
    constraint_type = atoi(argv[3]);
    congestion_parameter = atoi(argv[4]);
    alpha = atoi(argv[5]);
  }
  else {
    usage();
    return 0;
  }

  // Comienza lectura de archivo
  I = read_points(filename);
  if (I == NULL) {
    return 1;
  }

  IloInt n = I->n;
  IloNum S = I->S;
  point *puntos = I->points;
  // Determinar los ejes y los conjuntos N_i
  for (i = 0;i < n;i++) {
    puntos[i].N_i = new bool[n];
    for (j = 0;j < n;j++) puntos[i].N_i[j] = false;
  }

  bool **A;
  A = new bool*[n];
  for (i = 0;i < n;i++) {
    A[i] = puntos[i].N_i;
  }

  for (i = 0;i < n;i++) {
    A[i][i] = true;
    for (j = i + 1;j < n;j++) {
      if (dist(&(puntos[i]),&(puntos[j])) <= S) {
	A[i][j] = true;
	A[j][i] = true;
      }
    }
  }
  delete[] A;

  switch (n) {
  case 30:
    mu = 72;
    f = (constraint_type == WAITING_TIME ? 0.006 : 0.015);
    break;
  case 324:
  case 818:
  default:
    mu = 96;
    f = 0.01;
    break;
  }

  QM_CLAM(I,p,constraint_type,congestion_parameter,alpha,mu,f);
  
  cout << "QM-CLAM con p = " << p;
  if(constraint_type == WAITING_TIME)
    cout << " T = " << congestion_parameter;
  else
    cout << " b = " << congestion_parameter;
  cout << " alpha = " << alpha << endl;
  cout << "mu = " << mu << " f = " << f << " RHS = " 
       << (constraint_type == WAITING_TIME ? 
	   RHS_Waiting_time(mu,congestion_parameter,alpha) :
	   RHS_Queue_Size(mu,congestion_parameter,alpha)) << endl;

  for (i = 0;i < n;i++) delete[] puntos[i].N_i;
  delete[] puntos;
  delete I;

  return 0;
}	  

void QM_CLAM
(instance* I, // Set of points
 IloInt p, // facilities
 IloBool constraint_type, // queue size or waiting time
 IloNum congestion_parameter, // congestion parameter
 IloNum alpha, // minimum probability
 IloNum mu, // rate parameter
 IloNum f) // 
{
  IloEnv env;
  try {
    int i,j;
    IloInt n = I->n;
    point *puntos = I->points;
    // Comienza definicion del Modelo PMCLAP (o QM-CLAM)
#ifdef DEBUG
    ofstream log ("PMCLAP.log",ofstream::app);
    log  << "Comienza definicion de modelo" << endl;
#endif
    IloModel PMCLAP(env);
    
    // Variables
#ifdef DEBUG
    log << "Comienza creacion de variables" << endl;
#endif

    IloBoolVarArray y(env);
    BoolVarMatrix x(env,n);
    
    char VarName[16];
    for(i = 0;i < n;i++){
      sprintf(VarName,"y%d",i+1);
      y.add(IloBoolVar(env,VarName));
    }

#ifdef DEBUG
    log << "Termina creacion de variables $y_j$" << endl;
#endif
    // Nombre las variables para facil identificacion
    for(i = 0;i < n;i++){
      IloBoolVarArray x_i(env);
      for(j = 0;j < n;j++){
	sprintf(VarName,"x_%d_%d",i+1,j+1);
	x_i.add(IloBoolVar(env,VarName));
      }
      x[i] = x_i;
    }

    // ++ Restricciones ++
#ifdef DEBUG
    log << "Comienza definicion de restricciones" << endl;
#endif

    // Solo una instalacion puede servir al cliente en el nodo i
    bool *N_i;
    for(i = 0;i < n;i++){
      IloExpr Cover(env);
      N_i = puntos[i].N_i;
      for (j = 0;j < n;j++) {
	if (N_i[j]) Cover += x[i][j];
      }
      PMCLAP.add(Cover <= 1);
      Cover.end();
    }
#ifdef DEBUG
    log << "Termina restricciones $\\sum_j in N_i x_ij <= 1$" << endl;
#endif

    // Instalaciones a abrir
    PMCLAP.add(IloSum(y) == p);
#ifdef DEBUG
    log << "Termina restricciones $\\sum_j in  N y_i = p$" << endl;
#endif

    // Relacionar variables de localizacion y asignacion
    for(i = 0;i < n;i++){
      N_i = puntos[i].N_i;
      for (j = 0;j < n;j++) {
	if (N_i[j]) PMCLAP.add(x[i][j] <= y[j]);
      }
    }
#ifdef DEBUG
    log << "Termina restricciones $x_ij <= y_j$" << endl;    
#endif

    // Restricciones de espera
    IloNum fi;
    IloNum RHS;
    switch (constraint_type)
      {
      case QUEUE_SIZE:
#ifdef DEBUG
	log << "Comienza creacion de restricciones QUEUE SIZE" << endl;
#endif
	RHS = RHS_Queue_Size(mu,congestion_parameter,alpha);
	break;
      case WAITING_TIME:
#ifdef DEBUG
	log << "Comienza creacion de restricciones WAITING TIME" << endl;
#endif
	RHS = RHS_Waiting_time(mu,congestion_parameter,alpha);
	break;
      default:
	cerr << "Congestion type parameter unknown" << endl;
	break;
      }

    bool *N_j;
    for(j = 0;j < n;j++){
      N_j = puntos[j].N_i;
      IloExpr congestion_constraint(env);
      for (i = 0;i < n;i++) {
	if (N_j[i]) {
	  fi = f * puntos[i].demand;
	  congestion_constraint += fi * x[i][j];
	}
      }
      PMCLAP.add(congestion_constraint <= RHS);
      congestion_constraint.end();
    }

    // Variables iguales a 0
#ifdef DEBUG
    log << "Comienza creacion de restricciones iguales a 0" << endl;
#endif
    for (i = 0;i < n;i++) {
      N_i = puntos[i].N_i;
      for (j = i + 1;j < n;j++) {
	if(!N_i[j]) {
	  PMCLAP.add(x[i][j] == 0);
	  PMCLAP.add(x[j][i] == 0);
	}
      }
    }
    
    // Funcion Objetivo
#ifdef DEBUG
    log << "Crea funcion Objetivo" << endl;
#endif
    IloExpr Obj(env);
    for(i = 0;i < n;i++){
      N_i = puntos[i].N_i;
      for(j = 0;j < n;j++){
	if (N_i[j]) Obj += puntos[i].demand * x[i][j];
      }
    }
    PMCLAP.add(IloMaximize(env,Obj));
    Obj.end();

    // Resuelver modelo
    IloCplex cplex(PMCLAP);
    cplex.exportModel("QM-CLAM.lp");
    
    // Ajusta prioridades de variables
    for (j = 0;j < n;j++) {
      cplex.setPriority(y[j],1);
    }

    if(cplex.solve()){
      cout << "Solution status: " << cplex.getStatus() << endl;
      cout << "Maximum profit = " << cplex.getObjValue() << endl;
      for(j = 0;j < n;j++){
	if(cplex.getValue(y[j]) > 0) cout << j+1 << " ";
	//cout << cplex.getValue(y[j]);
      }
      cout << endl;
      /*for(i = 0;i < n;i++){
	cout << "-";
      }
      cout << endl;
      for(i = 0;i < n;i++){
	for(j = 0;j < n;j++){
	  cout << cplex.getValue(x[i][j]);
	}
	cout << endl;
      }/**/
      gnuplot(I,&cplex,&y,&x,p,constraint_type,congestion_parameter,alpha);
    }
    else {
      cout << "No solution found" << endl;
    }
    
  }
  catch (IloException& e) {
    cerr << "Concert exception caught: " << e << endl;
  }
  catch (...) {
    cerr << "Unknown exception caught" << endl;
  }
  env.end();
  
}

IloNum RHS_Queue_Size(IloNum mu,IloNum b,IloNum alpha){
  return mu * pow(1-alpha/100,1./(b+2));
}

IloNum RHS_Waiting_time(IloNum mu ,IloNum Tao,IloNum alpha){
  return mu + 60*24*log(1-alpha*0.01)/Tao;
}

void usage() {
  cout << "Usage:" << endl;
  cout << "\t./PMCLAP <filename> <p> <constraint_type>  <congestion_parameter> <alpha>" << endl;
  cout << "where:" << endl;
  cout << "\t<filename> is a file with the format:" << endl
       << "\t\tn S" << endl
       << "\t\tx1 y1 d1" << endl 
       << "\t\t..." << endl
       << "\t\txn yn dn" << endl;
  cout << "\t<p> is the number of centers to open" << endl;
  cout << "\t<constraint_type> is" << endl
       << "\t\t0: for queue size constraint" << endl
       << "\t\t1: for waiting time constraint" << endl;
  cout << "\t<congestion_parameter> is" << endl
       << "\t\tb the maximium queue lenght if queue size constraint was selected" << endl
       << "\t\tTao the maximum waiting time if waiting time constraint was selected" << endl;
  cout << "\t<alpha> is the minimum probability of at most" << endl
       << "\t\ta queue with b clients or" << endl
       << "\t\ta waiting time of Tao minutes" << endl;
}

void gnuplot(instance* I,IloCplex *cplex,IloBoolVarArray *y,BoolVarMatrix *x,IloInt p,IloBool constraint_type,IloNum congestion_parameter,IloNum alpha){
  IloInt i,j;
  IloInt n = I->n;
  point *puntos = I->points;

  char outfilename[32],centersfilename[32];
  FILE *gnuPipe = popen("gnuplot","w");

  sprintf(outfilename,"../../gnuplot/Ejes_%d.dat",n);
  ofstream outfile(outfilename);
  sprintf(centersfilename,"../../gnuplot/Centros_%d.dat",n);
  ofstream centros(centersfilename);
  
  for(i = 0;i < n;i++){
    if (cplex->getValue((*y)[i]))
      centros << puntos[i].x << " " << puntos[i].y << endl;
  }
  centros.close();

  for(i = 0;i < n;i++){
    for(j = 0;j < n;j++){
      if (i != j && puntos[i].N_i[j] && cplex->getValue((*x)[i][j]))
	outfile << puntos[i].x << " " << puntos[i].y << " " 
		<< puntos[j].x - puntos[i].x << " " << puntos[j].y - puntos[i].y << endl;
    }
  }
  outfile.close();

  fprintf(gnuPipe,"set term svg\n");
  fprintf(gnuPipe,"set output '../../gnuplot/PMCLAM_%d_%d_%d_%.0f_%.0f.svg'\n",n,p,constraint_type,congestion_parameter,alpha);
  fprintf(gnuPipe,"unset key\n");
  fprintf(gnuPipe,"unset border\n");
  fprintf(gnuPipe,"unset yzeroaxis\n");
  fprintf(gnuPipe,"unset xtics\n");
  fprintf(gnuPipe,"unset ytics\n");
  fprintf(gnuPipe,"unset ztics\n");
  fprintf(gnuPipe,"set title \"Servicio de %.0f\n",cplex->getObjValue());
  //fprintf(gnuPipe,"set style arrow 1 nohead lw 2\n");
  //fprintf(gnuPipe,"set arrow arrowstyle 1\n");
  fprintf(gnuPipe,"plot ");
  fprintf(gnuPipe,"'../PMCLAP/Instancias/Q_MCLP_%d.txt' every ::1 using 1:2 with points lc rgb \"black\"",n);
  fprintf(gnuPipe,", '../../gnuplot/Ejes_%d.dat' using 1:2:3:4 with vectors nohead linecolor rgb \"dark-blue\"",n);
  fprintf(gnuPipe,", '../../gnuplot/Centros_%d.dat' using 1:2 with points lc rgb \"red\"",n);

  fprintf(gnuPipe,"\n");
  pclose(gnuPipe);
  
}
