#include "headers.hpp"

extern int myRank;
extern int nbTasks;

Matrix Consistence(Matrix C, int ne, int np, Matrix Verts, Vector Xe, double he)
{
	Matrix Proj = Projection(C, ne, np, Verts, Xe, he);	

	Matrix G = LHS(C, ne, np, Verts, Xe, he);
	for (int i = 0; i < 3; ++i)
		 	for (int j = 0; j < G.cols(); ++j)
		 		G(i,j) = 0;

	return transpose(Proj)*G*Proj;
}

Matrix Stability(Matrix C, int ne, int np, Matrix Verts, Vector Xe, double he)
{
	Matrix Proj = Projection(C, ne, np, Verts, Xe, he);	
	Matrix D = dof(ne, np, basis, Verts, Xe, he);
	Matrix eye(2*ne,2*ne);

	for (int i = 0; i < eye.rows(); ++i)
		for (int j = 0; j < eye.cols(); ++j)
				eye(i,j) = (i==j) ? 1:0 ;

	return transpose(eye - D*Proj)*(eye - D*Proj);
}

Matrix RHS_P(Matrix C, int ne, int np, Matrix Verts, Vector Xe, double he)
{
	Matrix B(2*np,2*ne);
	for (int i = 0; i < B.rows(); ++i)
		for (int j = 0; j < B.cols(); ++j)
			B(i,j) = 0;
	Matrix Gs(3,2*np);
	for (int i = 0; i < Gs.rows(); ++i)
		for (int j = 0; j < Gs.cols(); ++j)
			Gs(i,j) = 0;
	Gs(2,3) = 2/he; Gs(0,4) = 1/he; Gs(1,5) = 1/he;
	Gs = C*Gs;
	for (int j = 0; j < 3; ++j)
	{
		IntVector ind1(2), ind2(2);
		ind1(0) = 0; ind1(1) = j;
		ind2(0) = 1; ind2(1) = j; 
		for (int k = 0; k < ne; ++k)
		{
			B(j, 2*k) = basis(ind1, Verts(k,0), Verts(k,1), Xe, he)/ne;
			B(j, 2*k+1) = basis(ind2, Verts(k,0), Verts(k,1), Xe, he)/ne;	
		}
	}
	Vector Prev(2), Next(2);
	Vector Vec_normal(2);
	Matrix Mat_normal(3,2);
	for (int i = 0; i < Mat_normal.rows(); ++i)
		for (int j = 0; j < Mat_normal.cols(); ++j)
			Mat_normal(i,j) = 0;
	Vector V_const1(2), V_const2(2), V_Gs(3);
	V_const1(0) = 0.5; V_const2(1) = 0.5; V_const1(1) = 0; V_const2(0) = 0;
	for (int k = 0; k < ne; ++k)
	{
		int kprev = (k-1)%ne, knext = (k+1)%ne;
        if (k == 0)
          kprev = ne-1;
		Prev(0) = Verts(kprev,0); Prev(1) = Verts(kprev,1);
		Next(0) = Verts(knext,0); Next(1) = Verts(knext,1);
		double avg = norm(Prev - Next);
		Vec_normal(0) = Next(1) - Prev(1); Vec_normal(1) = Prev(0) - Next(0);
		Vec_normal = Vec_normal/avg;
		Mat_normal(0,0) = Vec_normal(0); Mat_normal(2,0) = Vec_normal(1);
		Mat_normal(1,1) = Vec_normal(1); Mat_normal(2,1) = Vec_normal(0);
		for (int j = 3; j < 2*np; ++j)
		{
			for (int i = 0; i < 3; ++i)
				V_Gs(i) = Gs(i,j);
			Matrix V = Prod_vecMat(V_Gs, Mat_normal);
			B(j, 2*k) = B(j, 2*k) + avg*Prodsc(V_const1,V);
			B(j, 2*k+1) = B(j, 2*k+1) + avg*Prodsc(V_const2,V);
		}
	}

	return B;
}

Matrix LHS(Matrix C, int ne, int np, Matrix Verts, Vector Xe, double he)
{
	Matrix D = dof(ne, np, basis, Verts, Xe, he);
	Matrix B = RHS_P(C, ne, np, Verts, Xe, he);
	Matrix G = B*D;
	Matrix V(2*ne,2*np);
	for (int k = 0; k < ne; ++k)
	{
		IntVector Ind(2);
		for (int i = 0; i < 2; ++i){
			Ind(0) = i; 
			for (int j = 0; j < 2*np; ++j){
				Ind(1) = j;
				V(2*k+i,j) = basis(Ind, Verts(k,0), Verts(k,1), Xe, he);
			}
		}
	}
	Matrix Bs1(V.rows(),V.cols());
	Vector Vs1(V.cols()), Sum1(V.rows());
	for (int i = 0; i < V.cols(); ++i)
		Vs1(i) = 0;
	for (int i = 0; i < V.rows(); ++i)
		Sum1(i) = 0;

	for (int i = 0; i < 3; ++i)
	{
		for (int j = 0; j < V.rows(); ++j)
			for (int k = 0; k < V.cols(); ++k)
				Bs1(j,k) = V(j,i)*V(j,k);
		for (int j = 0; j < V.cols(); ++j)
		{
			for (int k = 0; k < V.rows(); ++k)
				Sum1(k) = Bs1(k,j);
			Vs1(j) = sum(Sum1)/ne;
			G(i,j) = Vs1(j);
		}
	}
	for (int i = 0; i < G.rows(); ++i)
		for (int j = 0; j < G.cols(); ++j)
			if (G(i,j) < 1e-12 && G(i,j) > -1e-12)
				G(i,j) = 0;

	return G;
}

Matrix dof(int ne, int np, double (*f)(IntVector Ind, double x, double y, Vector Xe, double he), Matrix Verts, Vector Xe, double he)
{
	Matrix D(2*ne, 2*np);
	IntVector Ind(2);
	for (int k = 0; k < ne; k++)
	{
		for (int i = 0; i < 2; i++)
		{
			Ind(0) = i;
			for (int j = 0; j < 2*np; j++)
			{
				Ind(1) = j;
				D(2*k+i, j) = f(Ind, Verts(k,0), Verts(k,1), Xe, he);
			}
		}
	}
	return D;
}

Matrix Projection(Matrix C, int ne, int np, Matrix Verts, Vector Xe, double he)
{
	Matrix B = RHS_P(C, ne, np, Verts, Xe, he);
	Matrix G = LHS(C, ne, np, Verts, Xe, he);
	Matrix Proj = inv(G)*B;

	return Proj;
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
