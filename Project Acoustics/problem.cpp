#include "headers.hpp"

extern int myRank;
extern int nbTasks;

//================================================================================
// Compute the matrices of the linear system
//================================================================================

void Acoustics(Problem& p, Mesh& m)
{
  if(myRank == 0)
    printf("== build linear system\n");
  
  int nd = m.nbOfTotEdges;
  int np = 2;
  cout.precision(4);
  p.M0.resize(nd, nd);
  p.M1.resize(nd, nd);
  p.A.resize(nd, nd);
  p.ML2.resize(nd,nd);
  p.U.resize(nd);
  p.F.resize(nd);
  for (int i = 0; i < nd; ++i){
    p.F(i) = 0;
    p.U(i) = 0;
    for (int j = 0; j < nd; ++j)
    {
      p.M0.coeffRef(i,j) = 0;
      p.M1.coeffRef(i,j) = 0;
      p.A.coeffRef(i,j) = 0;
    }
  }
  for (int iElt = 0; iElt < m.nbOfElems; ++iElt)
  {
    IntVector vtx_id = m.E2conn(iElt);
    Matrix Verts = eNod2Coords(m, iElt);
    int ne = vtx_id.rows();
    Vector Xe = centroid(Verts);
    double he = diameter(Verts);
    double Area = area(Verts);
    double sigma_E = 1;
    IntMatrix pol(2,2);
    pol(0,0) = 1; pol(0,1) = 0; 
    pol(1,0) = 0; pol(1,1) = 1;
    Matrix M0_E = Consistence(ne, np, Verts, Xe, Area, he);
    Matrix M1_E = Stability(ne, np, pol, Verts, Xe, Area, he, sigma_E);
    Matrix Mass_E = M0_E + M1_E;
    Matrix A_E = Divmatrix(ne, Area);
    //Vector F_E = Source(ne, Xe, Area);
    for (int i = 0; i < ne; ++i)
    {
      int I = abs(vtx_id(i))-1;
      for (int j = 0; j < ne; ++j)
      {
        int J = abs(vtx_id(j))-1;
        double sign = vtx_id(i)*vtx_id(j)/((I+1)*(J+1));
        p.M0.coeffRef(I,J) = p.M0.coeffRef(I,J) + sign*Mass_E(i,j);
        p.M1.coeffRef(I,J) = p.M1.coeffRef(I,J) + sign*(A_E(i,j) + Mass_E(i,j));
      }
    //p.F(I) = p.F(I) + F_E(i);
    }
  
  }
  double alpha = 1;
  p.A = p.M1 + alpha*p.M0;
  //p.F = p.M0*p.F;
  //exchangeAddInterfMPI(p.F, m);
}

Vector Source(int ne, Vector Xe, double Area)
{
  Vector F_E(ne);

  double u = (1+2*M_PI*M_PI)*sin(M_PI*Xe(0))*sin(M_PI*Xe(1));

  for (int i = 0; i < ne; ++i)
    F_E(i) = u*Area/ne;

  return F_E;
}

Vector SolEx(Mesh m)
{
  Vector Uex(m.nbOfNodes);
  for (int i = 0; i < m.nbOfNodes; ++i)
    Uex(i) = sin(M_PI*m.coords(i,0))*sin(M_PI*m.coords(i,1));
  
  return Uex;
}
//================================================================================
// Special treatment of the matrices for a Dirichlet boundary condition
//================================================================================

void buildDirichletBC(Problem& p, Mesh& m, Vector& uExa)
{
  if(myRank == 0)
    printf("== build Dirichlet BC\n");
  
  // Build 'relevement' function and mask (0 = interior node, 1 = boundary node)
  Vector g(m.nbOfNodes);
  IntVector IndBnd(m.nbOfBndNodes);
  for (int i = 0; i < m.nbOfBndNodes; ++i)
    IndBnd(i) = m.BndNodes(i)-1;

  for(int nLoc=0; nLoc<m.nbOfNodes; ++nLoc)
    g(nLoc) = 0.;
  for(int iBnd=0; iBnd<m.nbOfBndNodes; ++iBnd)
    g(IndBnd(iBnd)) = uExa(IndBnd(iBnd));
  
  // Modification of the system with 'relevement'
  Vector Ag = p.A*g;
  p.F = p.F - Ag;
  
  // Pseudo-reduction
  for(int iBnd=0; iBnd<m.nbOfBndNodes; ++iBnd){
    for(int i=0; i<m.nbOfNodes; ++i){
      if(IndBnd(iBnd) == i){
        for(int j=0; j<m.nbOfNodes; ++j){
          if(i == j){
            p.A.coeffRef(i, j) = 1.;
            p.F(i) = g(i);
          }
          else
            p.A.coeffRef(i, j) = 0.;
        }
      }
    }
  }
}
