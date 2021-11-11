#include "headers.hpp"

double basis(IntVector Ind, double x, double y, Vector Xe, double he)
{
	Matrix p(2,6);
	 
	p(0,0) = 1;					p(1,0) = 0;
	p(0,1) = 0;					p(1,1) = 1;
	p(0,2) = -(y - Xe(1))/he;	p(1,2) = (x - Xe(0))/he;
	p(0,3) = (y - Xe(1))/he;	p(1,3) = (x - Xe(0))/he;
	p(0,4) = (x - Xe(0))/he;	p(1,4) = 0;
	p(0,5) = 0;					p(1,5) = (y - Xe(1))/he;
	
	return p(Ind(0),Ind(1));
}

double area(Matrix Elt)
{
	return 0.5*fabs(sum(switch_prod(Elt)));
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