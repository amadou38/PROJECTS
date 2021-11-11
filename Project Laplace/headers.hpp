#ifndef HEADERS_HPP
#define HEADERS_HPP

/*#include <Eigen/Eigen>
#include <Spectra/SymEigsShiftSolver.h>
#include <iostream>
#include <fstream>
#include <mpi.h>
#include <stdio.h>
#include <math.h>
#include <Eigen/SVD>

using namespace Spectra;*/

#include <iostream>
#include <fstream>
#include <Eigen/Eigen>
#include <mpi.h>
#include <math.h>
#include <stdio.h>

using namespace std;

//================================================================================
// SPECIAL TYPES
//================================================================================

// Types for dense storage
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> Vector;
typedef Eigen::Matrix<double, 1, Eigen::Dynamic> Vector2;
typedef Eigen::Matrix<int,    Eigen::Dynamic, 1> IntVector;
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Matrix;
typedef Eigen::Matrix<int,    Eigen::Dynamic, Eigen::Dynamic> IntMatrix;

// Type for sparse storage
typedef Eigen::SparseMatrix<double> SpMatrix;
typedef Eigen::SparseMatrix<int> SpIntMatrix;

// Stucture for mesh
struct Mesh
{

  // Infos about nodes
  int nbOfNodes;  // number of nodes
  Matrix coords;  // coordinates of each node  (Size: nbOfNodes x 3)
  
  // Infos about elements
  int nbOfElems;              // number of elements
  int nbOfEdges;              // number of edges
  int nbOfEneighs;            // number of element neighboords
  int nbOfConns;              // number of connections (normal edges)
  IntVector nbOfEachElts;     // number of each types of elements
  IntMatrix elemNum;          //
  SpIntMatrix elemNodes;      // associated nodes (max 12. for triangle)
  IntMatrix edgNum;
  SpIntMatrix edgNodes;      // associated nodes (max 3. for triangle)
  IntMatrix eneighNum;
  SpIntMatrix eneighNodes;      // associated nodes (max 3. for triangle)
  IntMatrix connNum;
  SpIntMatrix connNodes;      // associated nodes (max 3. for triangle)
  // We are here
  int nbOfBndNodes;           // number of nodes on the boundary
  int nbOfIntNodes;           // number of interior nodes
  int nbOfFreeedg;            // number of edges from the boundary

  IntVector BndNodes;         // nodes from the boundary      (Size: nbOfBndNodes)
  IntVector IntNodes;         // Interior nodes               (Size: nbOfIntNodes)
  IntVector Freeedg;          // edges from the boundary      (Size: nbOfFreeedg)

  // Infos for parallel computations
  IntVector numNodesToExch;   // number of nodal values to exchanges between the current proc and each other proc  (Size: nbTasks)
  IntMatrix nodesToExch;      // list of nodal values to exchanges between the current proc and each other proc    (Size: nbTasks x max(numNodesToExch) )

  IntVector Elem2Nodes(int elt){
    IntVector Elem(elemNum(elt,1));
    for (int i = 0; i < elemNum(elt,1); ++i)
      Elem(i) = elemNodes.coeffRef(elt,i);
    return Elem;
  }
  IntVector E2edg(int elt){
    IntVector Elem(edgNum(elt,1));
    for (int i = 0; i < edgNum(elt,1); ++i)
      Elem(i) = edgNodes.coeffRef(elt,i);
    return Elem;
  }
  IntVector E2neigh(int elt){
    IntVector Elem(eneighNum(elt,1));
    for (int i = 0; i < eneighNum(elt,1); ++i)
      Elem(i) = eneighNodes.coeffRef(elt,i);
    return Elem;
  }
  IntVector E2conn(int elt){
    IntVector Elem(connNum(elt,1));
    for (int i = 0; i < connNum(elt,1); ++i)
      Elem(i) = connNodes.coeffRef(elt,i);
    return Elem;
  }
};

// Stucture for problem
struct Problem
{
  SpMatrix M0;   // stiffness matrix
  SpMatrix M1;   // mass matrix
  SpMatrix ML2;
  SpMatrix A;   // system matrix
  Vector U;
  Vector F;     // RHS
};

//================================================================================
// FUNCTIONS
//================================================================================

//==== Functions 

// Function to find the double integral value
double F(double(*f)(int ind, double x, double y, Vector Xe, double he), double x, double y, int indX, int indY, Vector Xe, double he);
// double givenFunction(double x, double y);

double basis(int ind, double x, double y, Vector Xe, double he);
void quadrature(int i, Matrix& Pts, Vector& Ws);
void Quad2D(int prec, Matrix T, Vector& W, Matrix& X);
double det(Matrix A, int N);
Matrix transpose(Matrix B);
Matrix inv1(Matrix A);
double integral_fct(double(*f)(int ind, double x, double y, Vector Xe, double he),  int indX, int indY, Matrix Q, Vector Xe, double he, int prec);
double sum(Vector X);
double trace(Matrix A);
double area(Matrix Elt);
Vector switch_prod(Matrix X);
Vector Abs(Vector X);
Vector centroid(Matrix Elt);
Vector Prod_vec(Vector X, Vector Y);
Vector Prod_vecMat(Vector X, Matrix M);
double max(Vector X);
double norm(Vector X);
double Prodsc(Vector X, Vector Y);
double diameter(Matrix Elt);
Matrix eNod2Coords(Mesh& m, int elt);
Matrix Consistence(int ne, int np, IntMatrix pol, Matrix Verts, Vector Xe, double he);
Matrix Stability(int ne, int np, IntMatrix pol, Matrix Verts, Vector Xe, double he);
Matrix RHS_P(int ne, int np, IntMatrix pol, Matrix Verts, Vector Xe, double he);
Matrix LHS(int np, Matrix Verts, Vector Xe, double he);
Matrix dof(int ne, int np, IntMatrix pol, Matrix Verts, Vector Xe, double he);
Matrix Projection(int ne, int np, IntMatrix pol, Matrix Verts, Vector Xe, double he);
Vector Source(int ne, Vector Xe, double Area);
Vector SolEx(Mesh m);
double getDeterminant(Matrix vect);
Matrix getTranspose(Matrix matrix1);
Matrix getCofactor(Matrix vect);
Matrix inv(Matrix vect);

//==== Functions in 'mesh.cpp'

// Read the mesh from a gmsh-file (.msh) and store in a mesh-structure 'm'
void readMsh(Mesh& m, string fileName);
void readInfo(int& nbOf, IntMatrix& Num, SpIntMatrix& Nodes, string fileName, string startName);
void resizeElts(Mesh& m);
// Write a solution 'u' in a gmsh-file (.msh)
void saveToMsh(Vector& u, Mesh& m, string name, string fileName);
void saveSol(Vector& u, string fileName);
//==== Functions in 'parallel.cpp'

// Build the local numbering and list of nodes for MPI communications
void buildListsNodesMPI(Mesh& m);

// MPI-parallel exchange/add the interface terms
void exchangeAddInterfMPI(Vector& vec, Mesh& m);
void exchangeAddInterfMPI2(Vector& vec, Mesh& m);

Vector RowMatrixVectorMultiply(Matrix& A, Vector& b);
//==== Functions in 'problem.cpp'

// Compute the matrices of the linear system
void Laplace(Problem& p, Mesh& m);

// Special treatment of the matrices for a Dirichlet boundary condition
void buildDirichletBC(Problem& p, Mesh& m, Vector& uExa);

//==== Functions in 'solver.cpp'

// Solution of the system Au=b with Jacobi
void jacobi(SpMatrix& A, Vector& b, Vector& u, Mesh& m, Vector& w, double tol, int maxit);
void GC(SpMatrix& A, Vector& b, Vector& u, Mesh& m, Vector& w, double tol, int maxit);
double scalarProductParallel(Vector& a, Vector& b, Vector& w);
Vector initWeightTab(Problem& p, Mesh& m);

#endif /* HEADERS_HPP */