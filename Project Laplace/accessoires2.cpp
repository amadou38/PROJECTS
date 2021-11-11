#include "headers.hpp"

extern int myRank;
extern int nbTasks;

Matrix eNod2Coords(Mesh& m, int elt)
{
	Matrix Coords(m.elemNum(elt,1),2);
	IntVector Elem = m.Elem2Nodes(elt);
	for (int i = 0; i < Elem.rows() ; ++i)
	{
		Coords(i,0) = m.coords(Elem(i)-1,0);
		Coords(i,1) = m.coords(Elem(i)-1,1);
	}

	return Coords;
}

void quadrature(int i, Matrix& Pts, Vector& Ws)
{
	if (i == 1)
	{
		Pts(i-1,0) = 1./3; Pts(i-1,1) = 1./3;
		Ws(i-1) = 0.5;
	}
	if (i == 3)
	{
		double a = 2./3, b = 1./6;
		Pts(0,0) = b; Pts(0,1) = b; 
		Pts(1,0) = a; Pts(1,1) = b; 
		Pts(2,0) = b; Pts(2,1) = a;
		for (int k = 0; k < i; k++)
		 	Ws(k) = 1./6;
	}
}

void Quad2D(int prec, Matrix T, Vector& W, Matrix& X)
{
	cout.precision(2);
	Matrix Pts(prec,2); Vector Ws(prec); 
	quadrature(prec, Pts, Ws);
	Vector A(2);
	Matrix ABC(2,2);

	for (int i = 0; i < T.rows(); i++)
	{
		A(i) = T(i,0);
		ABC(i,0) = T(i,1) - T(i,0);
		ABC(i,1) = T(i,2) - T(i,0);
	}
	
	X = A + ABC*transpose(Pts);
	W = abs(det(ABC, ABC.rows()))*Ws;
	X = transpose(X);
	
}

Vector switch_prod(Matrix X)
{
	Vector SP(X.rows());

	for (int i = 0; i < X.rows()-1; i++)
		SP(i) = X(i,0)*X(i+1,1) - X(i+1,0)*X(i,1);
	SP(X.rows()-1) = X(X.rows()-1,0)*X(0,1) - X(0,0)*X(X.rows()-1,1);

	return SP;
}

double F(double(*f)(int ind, double x, double y, Vector Xe, double he), double x, double y, int indX, int indY, Vector Xe, double he)
{
	return f(indX, x, y, Xe, he)*f(indY, x, y, Xe, he);
}

double sum(Vector X)
{
	double s = 0;

	for (int i = 0; i < X.rows(); i++)
		s += X(i);

	return s;
}

double trace(Matrix A)
{
	double tr = 0;
	for (int i = 0; i < A.rows(); ++i)
		tr += A(i,i);

	return tr;
}

double max(Vector X)
{
	double m = X(0);

	for (int i = 1; i < X.rows(); ++i)
		if (m < X(i))
			m = X(i);

	return m;
}

double norm(Vector X)
{
	double a = 0;
	for (int i = 0; i < X.rows(); i++)
		a += X(i)*X(i);
	return pow(a,0.5);
}

Matrix transpose(Matrix B)
{
    Matrix A(B.cols(), B.rows());
    for (int i = 0; i < B.cols(); i++)
        for (int j = 0; j < B.rows(); j++)
            A(i,j) = B(j,i);

    return A;
}

Vector Prod_vec(Vector X, Vector Y)
{
	Vector P(X.rows());

	for (int i = 0; i < X.rows(); i++)
		P(i) = X(i)*Y(i);

	return P;
}

Vector Prod_vecMat(Vector X, Matrix M)
{
	Vector P(M.cols());

	for (int j = 0; j < M.cols(); j++)
		for (int i = 0; i < M.rows(); ++i)
			P(j) += X(i)*M(i,j);

	return P;
}

double Prodsc(Vector X, Vector Y)
{
	double sc = 0;
	for (int i = 0; i < X.rows(); ++i)
		sc += X(i)*Y(i);

	return sc;
}

double integral_fct(double(*f)(int ind, double x, double y, Vector Xe, double he),  int indX, int indY, Matrix Q, Vector Xe, double he, int prec)
{
	double Integ = 0;
	Matrix T(2,3);
	T(0,2) = Xe(0); T(1,2) = Xe(1);
	Matrix X; Vector W;
	double FF = 0;
	for (int i = 0; i < Q.rows()-1; i++)
	{
		T(0,0) = Q(i,0); T(1,0) = Q(i,1);
		T(0,1) = Q(i+1,0); T(1,1) = Q(i+1,1);
		Quad2D(prec, T, W, X);
        
        FF = F(f, X(0,0), X(0,1), indX, indY, Xe, he);
        
		Integ += sum(W*FF);
	}
	T(0,0) = Q(0,0); T(1,0) = Q(0,1);
	T(0,1) = Q(Q.rows()-1,0); T(1,1) = Q(Q.rows()-1,1);
	
	Quad2D(prec, T, W, X);
	
    FF = F(f, X(0,0), X(0,1), indX, indY, Xe, he);
    
	Integ += sum(W*FF);

	return Integ;
}