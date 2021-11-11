#include "headers.hpp"

extern int myRank;
extern int nbTasks;

//================================================================================
// Compute the matrices of the linear system
//================================================================================

void elasto(Problem& p, Mesh& m, Matrix C)
{
  if(myRank == 0)
    printf("== build linear system\n");
  
  int nd = 2*m.nbOfNodes;
  int np = 3;

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
    IntVector vtx_id = m.Elem2Nodes(iElt);
    Matrix Verts = eNod2Coords(m, iElt);
    int ne = Verts.rows();
    Vector Xe = centroid(Verts);
    double he = diameter(Verts);
    double Area = area(Verts);
    Matrix M0_E = Consistence(C, ne, np, Verts, Xe, he);
    Matrix M1_E = 0.5*trace(M0_E)*Stability(C, ne, np, Verts, Xe, he);
    Vector F_E = Source(C, ne, Xe, Area);
  
    for (int i = 0; i < ne; ++i)
    {
      int I = vtx_id(i)-1, Ip = m.nbOfNodes + I; 
      for (int j = 0; j < ne; ++j)
      {
        int J = vtx_id(j) - 1, Jp = m.nbOfNodes + J;
        p.M0.coeffRef(I,J) = p.M0.coeffRef(I,J) + M0_E(2*i,2*j);
        p.M0.coeffRef(Ip,J) = p.M0.coeffRef(Ip,J) + M0_E(2*i+1,2*j);
        p.M0.coeffRef(I,Jp) = p.M0.coeffRef(I,Jp) + M0_E(2*i,2*j+1);
        p.M0.coeffRef(Ip,Jp) = p.M0.coeffRef(Ip,Jp) + M0_E(2*i+1,2*j+1);
        p.M1.coeffRef(I,J) = p.M1.coeffRef(I,J) + M1_E(2*i,2*j);
        p.M1.coeffRef(Ip,J) = p.M1.coeffRef(Ip,J) + M1_E(2*i+1,2*j);
        p.M1.coeffRef(I,Jp) = p.M1.coeffRef(I,Jp) + M1_E(2*i,2*j+1);
        p.M1.coeffRef(Ip,Jp) = p.M1.coeffRef(Ip,Jp) + M1_E(2*i+1,2*j+1);
      /*p.ML2.coeffRef(I,J) = p.ML2.coeffRef(I,J) + ML2_E(2*i,2*j);
        p.ML2.coeffRef(Ip,J) = p.ML2.coeffRef(Ip,J) + ML2_E(2*i+1,2*j);
        p.ML2.coeffRef(I,Jp) = p.ML2.coeffRef(I,Jp) + ML2_E(2*i,2*j+1);
        p.ML2.coeffRef(Ip,Jp) = p.ML2.coeffRef(Ip,Jp) + ML2_E(2*i+1,2*j+1);
      */
      }
    p.F(I) = p.F(I) + F_E(2*i);
    p.F(Ip) = p.F(Ip) + F_E(2*i+1);
    }
  }
  p.A = p.M0 + p.M1;
  exchangeAddInterfMPI(p.F, m);
}

Vector Source(Matrix C, int ne, Vector Xe, double Area)
{
  Vector F_E(2*ne);

  double u1 = M_PI*M_PI*sin(M_PI*Xe(0))*sin(M_PI*Xe(1));
  double u = M_PI*M_PI*cos(M_PI*Xe(0))*cos(M_PI*Xe(1));
  double f1 = u1*(C(0,0) + C(2,2)) - u*(C(0,1) + C(2,2));
  double f2 = u1*(C(1,1) + C(2,2)) - u*(C(0,1) + C(2,2));

  for (int i = 0; i < ne; ++i)
  {
    F_E(2*i) = f1*Area/ne;
    F_E(2*i+1) = f2*Area/ne;
  }

  return F_E;
}

Vector SolEx(Mesh m)
{
  Vector Uex(2*m.nbOfNodes);
  for (int i = 0; i < m.nbOfNodes; ++i)
  {
    Uex(i) = sin(M_PI*m.coords(i,0))*sin(M_PI*m.coords(i,1));
    Uex(m.nbOfNodes + i) = sin(M_PI*m.coords(i,0))*sin(M_PI*m.coords(i,1));
  }
  
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
  Vector g(2*m.nbOfNodes);
  IntVector IndBnd(2*m.nbOfBndNodes);
  for (int i = 0; i < m.nbOfBndNodes; ++i)
  {
    IndBnd(i) = m.BndNodes(i)-1;
    IndBnd(m.nbOfBndNodes+i) = m.nbOfNodes+IndBnd(i);
  }
  for(int nLoc=0; nLoc<2*m.nbOfNodes; ++nLoc)
    g(nLoc) = 0.;
  for(int iBnd=0; iBnd<2*m.nbOfBndNodes; ++iBnd)
    g(IndBnd(iBnd)) = uExa(IndBnd(iBnd));
  
  // Modification of the system with 'relevement'
  Vector Ag = p.A*g;
  p.F = p.F - Ag;
  
  // Pseudo-reduction
  for(int iBnd=0; iBnd<2*m.nbOfBndNodes; ++iBnd){
    for(int i=0; i<2*m.nbOfNodes; ++i){
      if(IndBnd(iBnd) == i){
        for(int j=0; j<2*m.nbOfNodes; ++j){
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
