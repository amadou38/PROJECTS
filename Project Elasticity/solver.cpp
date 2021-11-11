#include "headers.hpp"

extern int myRank;
extern int nbTasks;

//================================================================================
// Solution of the system Au=b with Jacobi
//================================================================================

void jacobi(SpMatrix& A, Vector& b, Vector& u, Mesh& m, Vector& w, double tol, int maxit)
{
  if(myRank == 0)
    printf("== jacobi\n");
  
  // Compute the solver matrices
  int size = A.rows();
  Vector Mdiag(size);
  SpMatrix N(size, size);
  for(int k=0; k<A.outerSize(); ++k){
    for(SpMatrix::InnerIterator it(A,k); it; ++it){
      if(it.row() == it.col())
        Mdiag(it.row()) = it.value();
      else
        N.coeffRef(it.row(), it.col()) = -it.value();
    }
  }
  exchangeAddInterfMPI(Mdiag, m);
  
  // Jacobi solver
  double bNorm = sqrt(scalarProductParallel(b, b, w));
  double residuNorm = tol*2*bNorm;
  int it = 0;
  Vector Nu, Au, r;
  //while (it < maxit){
  while (residuNorm > tol * bNorm && it < maxit){
    
    // Compute N*u
    Nu = N*u;
    exchangeAddInterfMPI(Nu, m);
    // Nu = RowMatrixVectorMultiply(N, u);
    // Update field
    u = Mdiag.cwiseInverse().cwiseProduct(Nu+b);
    
    // Update residual and iterator
    if((it % 100) == 0){
      Au = A*u;
      exchangeAddInterfMPI(Au, m);
      // Au = RowMatrixVectorMultiply(A, u);
      r = b - Au;
      residuNorm = sqrt(scalarProductParallel(r, r, w));
      if(myRank == 0)
        printf("\r%i %e\n", it, residuNorm);
    }
    it++;
  }
  
  if(myRank == 0){
    printf("\r   -> final iteration: %i (prescribed max: %i)\n", it, maxit);
    printf("   -> final residual: %e (prescribed tol: %e)\n", residuNorm, tol);
  }
}

double scalarProductParallel(Vector& a, Vector& b, Vector& weight){
    double scalarProductLocal = 0;
    double scalarProduct = 0;
    
    for(int i=0; i<a.rows(); i++){
       scalarProductLocal += a(i)*b(i) / (weight(i)+1);
    }

    MPI_Allreduce(&scalarProductLocal, &scalarProduct, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return scalarProduct;
}

void GC(SpMatrix& A, Vector& b, Vector& u, Mesh& m, Vector& w, double tol, int maxit)
{
  if(myRank == 0)
    printf("== Gradient conjugue\n");

  double alpha, beta, residuNorm2_k_1;
  Vector Ap;

  // Init residual
  Vector Au = A*u;
  exchangeAddInterfMPI(Au, m);
  Vector r = b - Au;
  double residuNorm2_k = scalarProductParallel(r, r, w);
  // Init p
  Vector p = r;

  double bNorm = sqrt(scalarProductParallel(b, b, w));
  
  // GC
  int it = 0;
  //while (it < maxit){
  while (sqrt(residuNorm2_k) > tol * bNorm && it < maxit){
    // Compute A*p
    Ap = A*p;
    exchangeAddInterfMPI(Ap, m);
    
    // Compute alpha_k
    alpha = residuNorm2_k/scalarProductParallel(Ap, p, w);

    // Update x_k to x_k+1
    u = u + alpha * p;

    // Update r_k to r_k+1
    r = r - alpha * Ap;

    // Compute beta_k
    residuNorm2_k_1 = scalarProductParallel(r, r, w);
    beta = residuNorm2_k_1/residuNorm2_k;

    // Update p_k to p_k+1
    p = r + beta * p;

    residuNorm2_k = residuNorm2_k_1;

    // Update residual and iterator
    if((it % 1000) == 0){
      if(myRank == 0)
        printf("%i %e\n", it, sqrt(residuNorm2_k));
    }
    it++;
  }
  if(myRank == 0){
    printf("\r   -> final iteration: %i (prescribed max: %i)\n", it, maxit);
    printf("   -> final residual: %e (prescribed tol: %e)\n", residuNorm2_k, tol);
  }
}

Vector initWeightTab(Problem& p, Mesh& m){
  int size = p.A.rows();
  int numToExch, nTask, nExch;
  Vector weight(size);
  for(int i=0; i<size; i++){
    weight(i) = 0;
  }
  for(nTask=0; nTask<nbTasks; nTask++){
    numToExch = m.numNodesToExch(nTask);
    if(numToExch > 0){
      for(nExch=0; nExch<numToExch; nExch++){
        weight(m.nodesToExch(nTask,nExch)-1) += 1;
        weight(m.nbOfNodes+m.nodesToExch(nTask,nExch)-1) += 1;
      }
    }
  }

  return weight;
}