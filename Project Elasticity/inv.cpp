#include "headers.hpp"

extern int myRank;
extern int nbTasks;

double getDeterminant(Matrix vect) {

    int dimension = vect.rows();

    if(dimension == 0) {
        return 1;
    }

    if(dimension == 1) {
        return vect(0,0);
    }

    //Formula for 2x2-matrix
    if(dimension == 2) {
        return vect(0,0) * vect(1,1) - vect(0,1) * vect(1,0);
    }

    double result = 0;
    int sign = 1;
    for(int i = 0; i < dimension; i++) {

        //Submatrix
        Matrix subVect(dimension - 1, dimension - 1);
        for(int m = 1; m < dimension; m++) {
            int z = 0;
            for(int n = 0; n < dimension; n++) {
                if(n != i) {
                    subVect(m-1,z) = vect(m,n);
                    z++;
                }
            }
        }

        //recursive call
        result = result + sign * vect(0,i) * getDeterminant(subVect);
        sign = -sign;
    }

    return result;
}

Matrix getTranspose(Matrix matrix1) {

    Matrix solution(matrix1.rows(), matrix1.cols());

    //Filling solution-matrix
    for(int i = 0; i < matrix1.rows(); i++) {
        for(int j = 0; j < matrix1.cols(); j++) {
            solution(j,i) = matrix1(i,j);
        }
    }
    return solution;
}

Matrix getCofactor(Matrix vect) {
    
    Matrix solution(vect.rows(),vect.cols());
    Matrix subVect(vect.rows() - 1, vect.cols() - 1);

    for(int i = 0; i < vect.rows(); i++) {
        for(int j = 0; j < vect.cols(); j++) {

            int p = 0;
            for(int x = 0; x < vect.rows(); x++) {
                if(x == i) {
                    continue;
                }
                int q = 0;

                for(int y = 0; y < vect.cols(); y++) {
                    if(y == j) {
                        continue;
                    }

                    subVect(p,q) = vect(x,y);
                    q++;
                }
                p++;
            }
            solution(i,j) = pow(-1, i + j) * getDeterminant(subVect);
        }
    }
    return solution;
}
Matrix inv(Matrix vect) {
    
    double d = 1.0/getDeterminant(vect);
    Matrix solution(vect.rows(), vect.cols());

    for(int i = 0; i < vect.rows(); i++) {
        for(int j = 0; j < vect.cols(); j++) {
            solution(i,j) = vect(i,j); 
        }
    }

    solution = getTranspose(getCofactor(solution));

    for(int i = 0; i < vect.rows(); i++) {
        for(int j = 0; j < vect.cols(); j++) {
            solution(i,j) *= d;
        }
    }

    return solution;
}