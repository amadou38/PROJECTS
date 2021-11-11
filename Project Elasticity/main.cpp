#include "headers.hpp"

// Driver Code

int myRank;
int nbTasks;

int main(int argc, char* argv[])
{
  // 1. Initialize MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
  MPI_Comm_size(MPI_COMM_WORLD, &nbTasks);
  
  double t1 = MPI_Wtime();

  // 2. Read the mesh, and build lists of nodes for MPI exchanges, local numbering
  Mesh m;
  readMsh(m, "benchmark/DOM/dom1x545x1.msh");

  //cout.precision(2);
  cout << "\n" << endl;
      
  // Test constantes
  //E1 = 126; E2 = 11; G = 6.6; nu1 = 0.28; nu2 = (E2/E1)*nu1; c = 1-nu1*nu2;  % constants
  //E1 = 11.6; E2 = 0.9; G = 0.75; nu1 = 0.37; nu2 = (E2/E1)*nu1; c = 1-nu1*nu2; rho = 390;  // constants Epicea
  double E1 = 13.7, E2 = 2.24, G = 1.61, nu1 = 0.3, nu2 = (E2/E1)*nu1, c = 1-nu1*nu2, rho = 750;  // constants Hetre
  //E1 = 11.5; E2 = 0.47; G = 0.5; nu1 = 0.005; nu2 = (E2/E1)*nu1; c = 1-nu1*nu2; rho = 392;  // constants Epicea
  //E1 = 8.86; E2 = 0.54; G = 1.6; nu1 = 0.005; nu2 = (E2/E1)*nu1; c = 1-nu1*nu2; rho = 691;   // constants Sapin

  buildListsNodesMPI(m);

  // 3. Build problem (fields and system)
  double t2 = MPI_Wtime();
  Problem p;
  Matrix C(3,3);

  for (int i = 0; i < C.rows(); ++i)
    for (int j = 0; j < C.cols(); ++j)
        C(i,j) = 0;
  C(0,0) = E1/c; C(0,1) = E1*nu2/c; 
  C(1,0) = E2*nu1/c; C(1,1) = E2/c; C(2,2) = G;
  //cout << "Matrix C (stress): \n" << C << endl;
  
  Vector Uex = SolEx(m);

  elasto(p, m, C);
  buildDirichletBC(p, m, Uex);

  // 4. Solve problem  
  double t3 = MPI_Wtime();

  //Compute weight tab for scalar product
  Vector w = initWeightTab(p, m);
 
  double tol = 1e-6;
  int maxit = 1000;
  GC(p.A, p.F, p.U, m, w, tol, maxit);
  // jacobi(p.A, p.F, p.U, m, w, tol, maxit);
  double t4 = MPI_Wtime();
  Vector U1(m.nbOfNodes), U2(m.nbOfNodes), U1ex(m.nbOfNodes), U2ex(m.nbOfNodes);
  for (int i = 0; i < m.nbOfNodes; ++i)
  {
    U1(i) = p.U(i);
    U2(i) = p.U(m.nbOfNodes+i);
    U1ex(i) = Uex(i);
    U2ex(i) = Uex(m.nbOfNodes+i);
  }
  
  Vector uErr = p.U - Uex;
  Vector MuErr = uErr, MUex = Uex;
  Vector KuErr = p.A*uErr, KUex = p.A*Uex;
  exchangeAddInterfMPI(MuErr, m);
  exchangeAddInterfMPI(KuErr, m);
  exchangeAddInterfMPI(MUex, m);
  exchangeAddInterfMPI(KUex, m);
  double Linf = max(Abs(uErr))/max(Abs(Uex));
  double L2relative = norm(uErr)/norm(Uex);
  double t5 = MPI_Wtime();

  if (myRank == 0)
  {
    printf("\nGlobal time              -> %f s\n", t5-t1);
    printf("\nProblem building time    -> %f s\n", t3-t2);
    printf("\nProblem solving time     -> %f s\n", t4-t3);
    cout << "\nL2relative norm: " << L2relative << endl;
    cout << "\nLinf norm: " << Linf << endl;
  }
  
  // saveToMsh(U1, m, "solNum1", "benchmark/solNum1.msh");
  // saveToMsh(U1ex, m, "solRef1", "benchmark/solExa1.msh");
  // saveToMsh(U2, m, "solNum2", "benchmark/solNum2.msh");
  // saveToMsh(U2ex, m, "solRef2", "benchmark/solExa2.msh");
  // saveToMsh(uErr, m, "solErr", "benchmark/solErr.msh");

  saveSol(U1, "benchmark/solNum1.msh");
  saveSol(U1ex, "benchmark/solExa1.msh");
  saveSol(U2, "benchmark/solNum2.msh");
  saveSol(U2ex, "benchmark/solExa2.msh");
  saveSol(uErr, "benchmark/solErr.msh");
  
  // 6. Finilize MPI
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();

  return 0;
}