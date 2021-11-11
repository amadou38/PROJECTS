#include "headers.hpp"

extern int myRank;
extern int nbTasks;

Matrix Consistence(int ne, int np, IntMatrix pol, Matrix Verts, Vector Xe, double he)
{
	Matrix Proj = Projection(ne, np, pol, Verts, Xe, he);	
	Matrix H = LHS(np, Verts, Xe, he);

	return transpose(Proj)*H*Proj;
}

Matrix Stability(int ne, int np, IntMatrix pol, Matrix Verts, Vector Xe, double he)
{
	Matrix Proj = Projection(ne, np, pol, Verts, Xe, he);		
	Matrix D = dof(ne, np, pol, Verts, Xe, he);
	Matrix H = LHS(np, Verts, Xe, he);
	Matrix eye(ne,ne);
	
	for (int i = 0; i < eye.rows(); ++i)
		for (int j = 0; j < eye.cols(); ++j)
				eye(i,j) = (i==j) ? 1:0 ;

	return H(0,0)*transpose(eye - D*Proj)*(eye - D*Proj);
}

Matrix RHS_P(int ne, int np, IntMatrix pol, Matrix Verts, Vector Xe, double he)
{
	Matrix B(np, ne);
	for (int i = 0; i < B.cols(); ++i)
		B(0,i) = 1./ne;
	for (int k = 0; k < ne; ++k)
	{
		Vector Vert(2);
		Vert(0) = Verts(k,0); Vert(1) = Verts(k,1);
		int kprev = (k-1)%ne, knext = (k+1)%ne;
        if (k == 0)
          kprev = ne-1;
    Vector Prev(2), Next(2), Vec_normal(2);
		Prev(0) = Verts(kprev,0); Prev(1) = Verts(kprev,1);
		Next(0) = Verts(knext,0); Next(1) = Verts(knext,1);
		double avg = norm(Prev - Next);
		Vec_normal(0) = Next(1) - Prev(1); Vec_normal(1) = Prev(0) - Next(0);
		for (int j = 1; j < np; ++j)
		{
			Vector monomial_grad(2);
			monomial_grad(0) = pol(j,0)/he; monomial_grad(1) = pol(j,1)/he;
			B(j,k) = 0.5*Prodsc(monomial_grad, Vec_normal);
		}
	}

	return B;
}

Matrix LHS(int np, Matrix Verts, Vector Xe, double he)
{
	Matrix H(np,np);
	int prec = 1, prec1 = 0;
	for (int i = 0; i < np; ++i)
		for (int j = 0; j < np; ++j){
			if ((i == 0) && (j==0))
        prec1 = 1;
      else
        prec1 = prec;
			H(i,j) = integral_fct(basis, i, j, Verts, Xe, he, prec1);
		}

	return H;
}

Matrix dof(int ne, int np, IntMatrix pol, Matrix Verts, Vector Xe, double he)
{
	Matrix D(ne, np);
	for (int i = 0; i < D.rows(); ++i)
		D(i,0) = 1;
	for (int k = 0; k < ne; ++k)
	{
		Vector Vert(2);
		Vert(0) = Verts(k,0); Vert(1) = Verts(k,1);
		for (int j = 1; j < np; ++j)
		{
			Vector polv(2);
			polv(0) = pol(j,0); polv(1) = pol(j,1);
			D(k,j) = Prodsc(Vert - Xe, polv)/he;
		}
	}
	
	return D;
}

Matrix Projection(int ne, int np, IntMatrix pol, Matrix Verts, Vector Xe, double he)
{
	Matrix D = dof(ne, np, pol, Verts, Xe, he);
	Matrix B = RHS_P(ne, np, pol, Verts, Xe, he);
	Matrix H = LHS(np, Verts, Xe, he);
	Matrix Proj = inv(B*D)*B;

	Matrix C = H*Proj;
	Matrix Proj0 = inv(H)*C;

	return Proj0;
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
