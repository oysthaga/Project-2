#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <armadillo>

// Create a tridiagonal matrix tridiag(a,d,e) of size n*n, 
// from scalar input a, d, and e. That is, create a matrix where
// - all n-1 elements on the subdiagonal have value a
// - all n elements on the diagonal have value d
// - all n-1 elements on the superdiagonal have value e
arma::mat create_tridiagonal(int n, double a, double d, double e)
{
    // Start from identity matrix
    arma::mat A = arma::mat(n, n, arma::fill::eye);

    // Fill the first row (row index 0), e.g.
    A(0,0) = d;
    A(0,1) = e;

    // Loop that fills rows 2 to n-1 (row indices 1 to n-2)
    for (int i = 1; i <= n-2; i++)
    {
    A(i,i) = d;
    A(i,i-1) = a;
    A(i, i+1) = e;
    }
    // Fill last row (row index n-1)
    A(n-1, n-1) = d;
    A(n-1, n-2) = a;

    return A;
}

// Create a symmetric tridiagonal matrix tridiag(a,d,a) of size n*n
// from scalar input a and d.
arma::mat create_symmetric_tridiagonal(int n, double a, double d)
{
    // Call create_tridiagonal and return the result
    return create_tridiagonal(n, a, d, a);
}

int main()
{
    int N = 6; 
    int n = N+1;
    double h = 1./n;
    double a = -1./pow(h,2);
    double d = 2./pow(h,2);
    arma::mat A = create_symmetric_tridiagonal(N, a, d);

    // Test that the matrix is correct 
    // by saving in textfile. 
    A.save("A.txt", arma::raw_ascii);

    // Numeric eigenvalues and eigenvectors.
    arma::vec eigval;
    arma::mat eigvec;
    eig_sym(eigval, eigvec, A); 
    eigval = arma::normalise(eigval);
    eigvec = arma::normalise(eigvec);
    eigval.save("eigval.txt", arma::raw_ascii);
    eigvec.save("eigvec.txt", arma::raw_ascii);


    // Analytic eigenvalues and eigenvectors.
    arma::vec lam = arma::vec(N);
    arma::mat v = arma::mat(N,N);
    for (int i = 0; i <= N-1; i++)
    {
        lam[i] = d + 2*a*cos( ((i+1)*M_PI)/(N+1) );
        for (int j = 0; j <= N-1; j++)
        {
            v(i, j) = sin ( ((i+1)*(j+1)*M_PI) / (N+1) );
        }
    }

    lam = arma::normalise(lam);
    v = arma::normalise(v);

    lam.save("lam.txt", arma::raw_ascii);
    v.save("v.txt", arma::raw_ascii);
/*
    double tol = pow(10, -4);
    for (int i = 0; i <= N-1; i++)
    {
        if abs(lam[i]-eigval[i])>tol
        {
            for (int j = 0; j <= N-1; j++)
            {
                if abs(v(i,j)-eigvec(i,j))>tol
                {
                    ...
                }
            }
    }
*/
    return 0;
}