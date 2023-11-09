#include "Jacobi.hpp"


double max_offdiag(const arma::mat& A, int& k, int& l)
{
    int Len = A.n_rows;
    //copies the absolute values of the upper diagonal part of the matrix
    arma::mat B = trimatu(abs(A),1); 
    arma::uword i = B.index_max(); //linear index of max value


    l = static_cast<int>(i/Len); //Makes sure we do integer division and truncate the result.
    
    k = i - (Len * (l));

    return A(i);
}


void Jacobi_rotate(arma::mat& A, arma::mat& R, int k, int l, int& Len)
{
    //We calculate the different variables neccessary for doing the rotation (this is just math, dont mind it)
    double tau = (A(l,l) - A(k,k))/(2*A(k,l));

    double t;
    
    if (tau > 0)
    {
        t = 1/(tau + std::sqrt(1+ std::pow(tau, 2)));
    }

    else
    {   
        t = -1/(-tau + std::sqrt(1+ std::pow(tau, 2)));
    }

    // std::cout <<t <<std::endl;

    double c = 1/(std::sqrt(1+ std::pow(t,2)));
    double s = c*t;
    
    //we update the A and R matrix
    double Akk_temp = A(k, k); //a temporary variable because we rewrite A(k, k) and the use it
    A(k, k) = A(k, k)*c*c - 2*A(k, l)*c*s + A(l, l)*s*s;
    A(l, l) = A(l, l)*c*c + 2*A(k, l)*c*s + Akk_temp*s*s;
    A(k, l) = 0;
    A(l, k) = 0;

    for (int i = 0; i<Len; i++)
    {
        if (i == k || i == l )
        {
            ;
        }
        else 
        {
            double Ak_temp = A(i, k);
            double Al_temp = A(i, l);
            A(i, k) = Ak_temp*c - Al_temp*s;
            A(i, l) = Al_temp*c + Ak_temp*s;
            A(k, i) = A(i, k);
            A(l, i) = A(i, l);
        }

        double Rk_temp = R(i, k);
        double Rl_temp = R(i, l);
        R(i, k) = Rk_temp*c - Rl_temp*s;
        R(i, l) = Rl_temp*c + Rk_temp*s;
    
    }
}




void Jacobi_eigensolver(const arma::mat& A, double eps, arma::vec& eigenvalues, arma::mat& 
    eigenvectors, const int maxiter, int& iterations, bool& converged)
{
    //determine size of A and  test if square
    assert(true == A.is_square()); 
    
    int Len = A.n_rows;

    //Makes sure we don't make changes to A, but costs us memory
    arma::mat A_update = A; 

    //initialize the k and l vars.
    int k; int l;


    //start R as the identity matrix
    arma::mat  R(Len, Len, arma::fill::eye);



    //just starting a_max with something more that eps
    double a_max = eps + 1.;

    while(fabs(a_max) > eps)
    {
        //find max off-diag. element and its indexes
        a_max = max_offdiag(A_update, k, l);
        // std::cout << a_max << std::endl;

        if (k == l)
        {
            std::cout <<"Gave diag as a_max" << std::endl;
            break;
        }

        //Update A_update matrix and R matrix
        Jacobi_rotate(A_update, R, k, l, Len);

        // A_update.print();
        // std::cout << std::endl;


        //counts iterations
        iterations++;

        if (iterations > maxiter)
        {
            std::cout<< "We have reached our max number of iterations" <<std::endl;
            break;
        }
    }

    if (a_max < eps)
    {
        converged = true;
    }

    eigenvalues = A_update.diag();
    eigenvectors = R;

}


int main()
{
    int N = 5;

    //Generate random matrix
    arma::mat A(N,N, arma::fill::randn);
    A.print();
    std::cout << std::endl;

    //Declaring variables for the eigen-solver to update
    arma::vec eigenvalues;
    arma::mat eigenvectors;
    int iterations = 0;
    bool converged = false;


    //Finding the eigen values and eigen vectors of A
    Jacobi_eigensolver(A, 1e-8, eigenvalues, eigenvectors, 1e8, iterations, converged);

    std::cout << "Number of iterations:" << iterations << std::endl;
    std::cout << std::endl;

    std::cout << "Eigenvectors:" << std::endl;

    eigenvectors.print();
    std::cout << std::endl;

    std::cout << "Eigenvalues::" << std::endl;
    eigenvalues.print();
}