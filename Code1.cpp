#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <armadillo>
#include <cassert>

// Create a tridiagonal matrix tridiag(a,d,e) of size n*n, 
// from scalar input a, d, and e. 
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
    // Get size of the matrix A. 
    int size = A.n_rows;


    // Check that A is square and larger than 1x1. 
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
    double t;
    arma::mat a_m = A; // Copy
    arma::mat r_m = R; // Copy
    

    double tau = (A(l,l)-A(k,k))/(2*A(k,l));
    if (tau>0)
    {
        t = 1/(tau+sqrt(1+pow(tau,2)));
    }
    else
    {
        t = -1/(-tau+sqrt(1+pow(tau,2)));
    }
    double c = 1/sqrt( 1 + pow(t,2) );
    double s = c*t;

    A(k,k) = A(k,k)*pow(c,2) - 2*A(k,l)*c*s + A(l,l)*pow(s,2);
    A(l,l) = A(l,l)*pow(c,2) + 2*A(k,l)*c*s + A(k,k)*pow(s,2);
    A(k,l) = 0;
    A(l,k) = 0; 

    for (int i = 0; i <= N-1; i++)
    {
        if (i != k)
        {
            if (i != l) 
            {
                A(i,k) = a_m(i,k)*c - a_m(i,l)*s;
                A(k,i) = A(i,k); 
                A(i,l) = a_m(i,l)*c + a_m(i,k)*s; 
                A(l,i) = A(i,l); 
            }
        }
        R(i,k) = r_m(i,k)*c - r_m(i,l)*s;
        R(i,l) = r_m(i,l)*c + r_m(i,k)*s; 

    }
}

// Jacobi method eigensolver:
// - Runs jacobo_rotate until max off-diagonal element < eps
// - Writes the eigenvalues as entries in the vector "eigenvalues"
// - Writes the eigenvectors as columns in the matrix "eigenvectors"
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
    
    double Amax = max_offdiag_symmetric(A_m, k, l);
    iterations = 0;
    while ( Amax > eps ) // Rotate until max of-diagonal value is close to 0.
    {
        // run jacobi_rotate
        Amax = max_offdiag_symmetric(A_m, k, l); // Maximum off-diagonal value
        jacobi_rotate(A_m, R, k, l); // One rotation
        
        
        iterations += 1;
        if (iterations == maxiter) 
        {
            std::cout << "\nBroke\n";
            converged = false;
            break; 
        }
        if (Amax <= eps)
        {
            converged = true;
        }
    }
    arma::vec A_m_diag = A_m.diag();
    eigenvalues = arma::normalise(A_m_diag);
    eigenvectors = arma::normalise(R);
    arma::uvec indices = arma::sort_index(eigenvalues);
    eigenvalues = eigenvalues(indices);
    eigenvectors = eigenvectors.cols(indices);

}




int main(int argc, char* argv[])
{
  // Check number of command-line arguments
    if (argc != 2)  // Expects 1 command-line arguments
    {
        // Get the name of the executable file
        std::string executable_name = argv[0];
        std::cerr << "Error: Wrong number of input arguments." << std::endl;
        std::cerr << "Usage: " << executable_name << " <some integer> " << std::endl;
        // Exit program with non-zero return code to indicate a problem
        return 1;   
    }

    int N = atoi(argv[1]);    // atoi converts argv[2] to an integer
    int n = N+1;              // Number of steps
    double h = 1./n;          // Step length
    double a = -1./pow(h,2);  // Diagonal
    double d = 2./pow(h,2);   // Sub-/superdiagonal
    arma::mat A = create_symmetric_tridiagonal(N, a, d);


    // Numeric eigenvalues and eigenvectors.
    arma::vec eigval;
    arma::mat eigvec;
    eig_sym(eigval, eigvec, A); 
    eigval = arma::normalise(eigval);
    eigvec = arma::normalise(eigvec);


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

    // Prints the matrix and the eigenvalues and -vectors. 
    std::cout << "\nProblem 2\n";
    std::cout << "\nA\n";
    std::cout << A;
    std::cout << "\nAnalytic eigenvalues\n";
    std::cout << lam;
    std::cout << "\nNumeric eigenvalues\n";
    std::cout << eigval;
    std::cout << "\nAnalytic eigenvectors\n";
    std::cout << v;
    std::cout << "\nNumeric eigenvectors\n";
    std::cout << eigvec;


    // Prints the maximum off-diagonal value of a test-matrix. 
    int k, l; 
    arma::mat A_test = "1., 0., 0., 0.5; 0., 1., -0.7, 0.; 0., -0.7, 1., 0.; 0.5, 0., 0., 1.;";
    double test_max = max_offdiag_symmetric(A_test, k, l);
    std::cout << "\nProblem 3\n";
    std::cout << "Max off-diagonal-value: ";
    std::cout << test_max;
    std::cout << "\n";


    // Runs Jacobi eigensolver
    arma::vec eigenvalues;
    arma::mat eigenvectors;
    int iterations;
    bool converged;
    double epsilon = pow(10, -14);
    int max_iterations = 100000;
    jacobi_eigensolver(A, epsilon, eigenvalues, eigenvectors, max_iterations, iterations, converged);

    eigenvalues.save("eigenvalues.txt", arma::raw_ascii);
    eigenvectors.save("eigenvectors.txt", arma::raw_ascii);
    int k1, l1; 

    std::cout << "\nProblem 4\n";
    std::cout << "\nA\n";
    std::cout << A;
    std::cout << "\nMax off-diagonal element\n";
    std::cout << max_offdiag_symmetric(A,k1,l1);
    std::cout << "\nAnalytic eigenvalues\n";
    std::cout << lam;
    std::cout << "\nAnalytic eigenvectors\n";
    std::cout << v;
    std::cout << "\nNumeric eigenvalues\n";
    std::cout << eigenvalues;
    std::cout << "\nNumeric eigenvectors\n";
    std::cout << eigenvectors;
    if (converged == true)
    {
        std::cout << "\nConverged\n";
    }
    else
    {
        {
        std::cout << "\nNot converged\n";
    }
    }
    std::cout << "\niterations\n";
    std::cout << iterations;
    std::cout <<"\n";

    return 0;
}