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
  readMsh(m, "benchmark/mesh1.msh");
  //buildListsNodesMPI(m);
  //cout.precision(2);
  cout << "\n" << endl;

  buildListsNodesMPI(m);
  
  // 3. Build problem (fields and system)
  double t2 = MPI_Wtime();
  Problem p;

  //Compute weight tab for scalar product
  Acoustics(p, m);
  Vector Uex = SolEx(m);
  buildDirichletBC(p, m, Uex);
  
  // 4. Solve problem
  double t3 = MPI_Wtime();
  
  //Compute weight tab for scalar product
  Vector w = initWeightTab(p, m);

  double tol = 1e-6;
  int maxit = 1000;
  GC(p.A, p.F, p.U, m, w, tol, maxit);
  //jacobi(p.A, p.F, p.U, m, w, tol, maxit);
  double t4 = MPI_Wtime();

  Vector uErr = p.U - Uex;
  Vector MuErr = p.M0*uErr;
  Vector KuErr = p.M1*uErr;
  exchangeAddInterfMPI(MuErr, m);
  exchangeAddInterfMPI(KuErr, m);
  double L2Norm = sqrt(scalarProductParallel(uErr, MuErr, w));
  double H1Norm = sqrt(scalarProductParallel(uErr, KuErr, w));
  double t5 = MPI_Wtime();

  if (myRank == 0)
  {
    printf("\nGlobal time              -> %f s\n", t5-t1);
    printf("\nProblem building time    -> %f s\n", t3-t2);
    printf("\nProblem solving time     -> %f s\n", t4-t3);
    cout << "\nL2 norm: " << L2Norm << endl;
    cout << "\nH1 norm: " << H1Norm << endl;
  }
  
  saveToMsh(p.U, m, "solNum", "benchmark/solNum.msh");
  saveToMsh(Uex, m, "solRef", "benchmark/solExa.msh");
  saveToMsh(uErr, m, "solErr", "benchmark/solErr.msh");

  // 6. Finilize MPI
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();

  return 0;
}