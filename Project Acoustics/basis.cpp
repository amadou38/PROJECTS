#include "headers.hpp"

extern int myRank;
extern int nbTasks;

double basis(int ind, double x, double y, Vector Xe, double he)
{
	Vector p(3), Coord(2);
	Coord(0) = x; Coord(1) = y; 
	p(0) = (x - Xe(0))/he;
	p(1) = (y - Xe(1))/he;
	p(2) = 1;
	for (int i = 0; i < 2; ++i)
		p(i) = p(i) + 0/*integral_fct(basis, i, 2, Coord, Xe, he, 1)*/;
		
	return p(ind);
}

double area(Matrix Elt)
{
	return 0.5*abs(sum(switch_prod(Elt)));
}

Vector centroid(Matrix Elt)
{
	Vector Xe(Elt.cols());
	Vector X(Elt.rows());

	for (int j = 0; j < Elt.cols(); j++)
	{
		for (int i = 0; i < Elt.rows(); i++)
		{
			if (i == Elt.rows()-1)
				X(i) = Elt(i,j) + Elt(0,j);
			else
				X(i) = Elt(i,j) + Elt(i+1,j);
		}
		Xe(j) = sum(Prod_vec(X,switch_prod(Elt)));
		Xe(j) = Xe(j)/(6*0.5*sum(switch_prod(Elt)));
	}

	return Xe;
}

double diameter(Matrix Elt)
{
	double he = 0;
	Vector X(Elt.cols()), Y(Elt.cols());
	for (int i = 0; i < Elt.rows()-1; i++)
	{
		for (int k = 0; k < Elt.cols(); k++)
			X(k) = Elt(i,k);
		for (int j = i+1; j < Elt.rows(); j++)
		{
			for (int k = 0; k < Elt.cols(); k++)
				Y(k) = Elt(j,k);
			he = max(he, norm(X - Y));
		}
	}
return he;
}