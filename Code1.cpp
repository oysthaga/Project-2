#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <armadillo>
#include <cassert>

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


// A function that finds the max off-diag element of a symmetric matrix A.
// - The matrix indices of the max element are returned by writing to the  
//   int references k and l (row and column, respectively)
// - The value of the max element A(k,l) is returned as the function
//   return value
double max_offdiag_symmetric(const arma::mat& A, int& k, int& l)
{
    // Get size of the matrix A. Use e.g. A.n_rows, see the Armadillo documentation
    int size = A.n_rows;

    // Possible consistency checks:
    // Check that A is square and larger than 1x1. Here you can for instance use A.is_square(), 
    // see the Armadillo documentation.
    // 
    // The standard function 'assert' from <assert.h> can be useful for quick checks like this 
    // during the code development phase. Use it like this: assert(some condition),
    // e.g assert(a==b). If the condition evaluates to false, the program is killed with 
    // an assertion error. More info: https://www.cplusplus.com/reference/cassert/assert/
    assert(A.is_square());

    // Initialize references k and l to the first off-diagonal element of A
    k = 0;
    l = 1; 

    // Initialize a double variable 'maxval' to A(k,l). We'll use this variable 
    // to keep track of the largest off-diag element.
    double maxval = std::abs(A(k,l));

    // Loop through all elements in the upper triangle of A (not including the diagonal)
    // When encountering a matrix element with larger absolute value than the current value of maxval,
    // update k, l and max accordingly.
    for (int i = 0; i <= size-2; i++)
    {
        for (int j = i+1; j <= size-1; j++)
        {
            if ( std::abs(A(i,j))  > maxval )
            {
                maxval = std::abs(A(i,j));
                k=i;
                l=j;

            }
        }
    }

    return maxval;
}

// Performs a single Jacobi rotation, to "rotate away"
// the off-diagonal element at A(k,l).
// - Assumes symmetric matrix, so we only consider k < l
// - Modifies the input matrices A and R
void jacobi_rotate(arma::mat& A, arma::mat& R, int k, int l)
{
    int N = A.n_rows;
    arma::mat A_mp1 = arma::mat(N,N);
    arma::mat R_mp1 = arma::mat(N,N);
    double t;

    double tau = (A(l,l)-A(k,k))/(2*A(k,l));
    if (tau>0)
    {
        t = -tau + sqrt( 1+ pow(tau,2) );
    }
    else
    {
        t = -tau - sqrt( 1+ pow(tau,2) );
    }
    double c = 1/sqrt( 1 + pow(t,2) );
    double s = c*t;

    A_mp1(k,k) = A(k,k)*pow(c,2) - 2*A(k,l)*c*s + A(l,l)*pow(s,2);
    A_mp1(l,l) = A(l,l)*pow(c,2) + 2*A(k,l)*c*s + A(k,k)*pow(s,2);
    A_mp1(k,l) = 0;
    A_mp1(l,k) = 0; 

    for (int i = 0; i <= N-1; i++)
    {
        if (i != k)
        {
            if (i != l) 
            {
                A_mp1(i,k) = A(i,k)*c - A(i,l)*s;
                A_mp1(k,i) = A_mp1(i,k); 
                A_mp1(i,l) = A(i,l)*c + A(i,k)*s; 
                A_mp1(l,i) = A_mp1(i,l); 
            }
        }
        R_mp1(i,k) = R(i,k)*c - R(i,l)*s;
        R_mp1(i,l) = R(i,l)*c + R(i,k)*s; 
    }
    A = A_mp1;
    R = R_mp1;
}

// Jacobi method eigensolver:
// - Runs jacobo_rotate until max off-diagonal element < eps
// - Writes the eigenvalues as entries in the vector "eigenvalues"
// - Writes the eigenvectors as columns in the matrix "eigenvectors"
//   (The returned eigenvalues and eigenvectors are sorted using arma::sort_index)
// - Stops if the number of iterations reaches "maxiter"
// - Writes the number of iterations to the integer "iterations"
// - Sets the bool reference "converged" to true if convergence was reached before hitting maxiter
void jacobi_eigensolver(const arma::mat& A, double eps, arma::vec& eigenvalues, arma::mat& eigenvectors, 
                        const int maxiter, int& iterations, bool& converged)
{
    arma::mat A_m = A;
    int N = A.n_rows;
    arma::mat R = arma::eye(N,N);
    int k;
    int l;
    double Amax;
    while ( Amax > eps )
    {
        // run jacobi_rotate
        Amax = max_offdiag_symmetric(A_m, k, l);
        jacobi_rotate(A_m, R, k, l);
        
        
        iterations += 1;
        if (iterations = maxiter)
        {
            break; 
        }
        if (Amax <= eps)
        {
            converged = true;
        }
    }
    arma::vec A_m_diag = A_m.diag();
    //arma::uvec val_indices = arma::sort_index(A_m_diag);

    //arma::mat EigValsSorted = A_m_diag.sort(val_indices);
    //arma::mat EigValsSorted = sort(A_m_diag, val_indices);
    //eigenvalues = arma::sort(A_m_diag);
    //eigenvectors = arma::sort(R);
    eigenvalues = A_m_diag;
    eigenvectors = R;
/*
    arma::vec diag(N);

    for(int i=0, i<N, i++)
    {
        diag(i,i) = A_m(i, i)
    }
*/
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


    int k, l; 
/*
    arma::mat A_test = "1., 0., 0., 0.5; 0., 1., -0.7, 0.; 0., -0.7, 1., 0.; 0.5, 0., 0., 1.;";
    double test_max = max_offdiag_symmetric(A_test, k, l);
    std::cout << test_max;
*/
    arma::vec eigenvalues;
    arma::mat eigenvectors;
    int iterations;
    bool converged;
    double epsilon = pow(10, -12);
    int max_iterations = 100000;
    jacobi_eigensolver(A, epsilon, eigenvalues, eigenvectors, max_iterations, iterations, converged);

    std::cout << eigenvalues;
    std::cout << eigenvectors;
    std::cout << converged;


    return 0;
}