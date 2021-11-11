#include "headers.hpp"

extern int myRank;
extern int nbTasks;

Matrix Consistence(int ne, int np, Matrix Verts, Vector Xe, double Area, double he)
{
	Matrix Proj = Projection(ne, np, Verts, Xe, Area, he);
	Matrix Pab = LHS(np, Area, he);

	return transpose(Proj)*Pab*Proj;
}

Matrix Stability(int ne, int np, IntMatrix pol, Matrix Verts, Vector Xe, double Area, double he, double sigma_E)
{
	Matrix Proj = Projection(ne, np, Verts, Xe, Area, he);		
	Matrix D = dof(ne, np, pol, Verts, Xe, he);
	Matrix eye(ne,ne);
	
	for (int i = 0; i < eye.rows(); ++i)
		for (int j = 0; j < eye.cols(); ++j)
				eye(i,j) = (i==j) ? 1:0 ;

	return sigma_E*transpose(eye - D*Proj)*(eye - D*Proj);
}

Matrix RHS_P(int ne, int np, Matrix Verts, Vector Xe, double he)
{
	Matrix B(np, ne);

	for (int k = 0; k < ne; ++k)
	{
		Vector Vert(2);
		Vert(0) = Verts(k,0); Vert(1) = Verts(k,1);
		int knext = (k+1)%ne;
    Vector Next(2);
    Matrix Middle(1,2);
		Next(0) = Verts(knext,0); Next(1) = Verts(knext,1);
		Middle(0,0) = 0.5*(Vert(0) + Next(0)); Middle(0,1) = 0.5*(Vert(1) + Next(1));
		for (int j = 0; j < np; ++j)
			B(j,k) = basis(j, Middle(0,0), Middle(0,1), Xe, he) - integral_fct(basis, j, 2, Middle, Xe, he, 3);
	}

	return B;
}

Matrix Divmatrix(int ne, double Area)
{
	Matrix A_E(ne,ne);
	
	for (int i = 0; i < ne; ++i)
		for (int j = 0; j < ne; ++j)
			A_E(i,j) = 1./Area;

	return A_E;
}

Matrix LHS(int np, double Area, double he)
{
	Matrix Pab(np,np);

	for (int i = 0; i < np; ++i)
		for (int j = 0; j < np; ++j)
			Pab(i,j) = (i==j) ? Area/pow(he,2):0;

	return Pab;
}

Matrix dof(int ne, int np, IntMatrix pol, Matrix Verts, Vector Xe, double he)
{
	Matrix D(ne, np);
	
	for (int k = 0; k < ne; ++k)
	{
		Vector Vert(2);
		Vert(0) = Verts(k,0); Vert(1) = Verts(k,1);
		int knext = (k+1)%ne;
    Vector Next(2), Vec_normal(2);
		Next(0) = Verts(knext,0); Next(1) = Verts(knext,1);
		double avg = norm(Vert - Next);
		Vec_normal(0) = Next(1) - Vert(1); Vec_normal(1) = Vert(0) - Next(0);
		Vec_normal = Vec_normal/avg;
		for (int j = 0; j < np; ++j)
		{
			Vector polv(2);
			polv(0) = pol(j,0); polv(1) = pol(j,1);
			D(k,j) = avg*Prodsc(polv, Vec_normal)/he;
		}
	}
	
	return D;
}

Matrix Projection(int ne, int np, Matrix Verts, Vector Xe, double Area, double he)
{
	Matrix B = RHS_P(ne, np, Verts, Xe, he);
	Matrix G = LHS(np, Area, he);

	return inv(G)*B;
}

Vector Abs(Vector X)
{
	Vector Abs_X(X.rows());
	for (int i = 0; i < X.rows(); i++)
		Abs_X(i) = abs(X(i));

	return Abs_X;
}

double det(Matrix A, int N)
{
  double d = 0;

  if (N <= 2)
      d = A(0,0)*A(1,1) - A(0,1)*A(1,0);
  else
    for(int i = 0; i < N; i++)
      d += (A(0,i) * (A(1,(i+1)%N) * A(2,(i+2)%N) - A(1,(i+2)%N) * A(2,(i+1)%N)));

    return d;
}
