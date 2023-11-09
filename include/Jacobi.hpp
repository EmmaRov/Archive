#include <iostream>
#include <assert.h>
#include <armadillo>
#include <cmath>
#include <string>
#include <fstream>
#include <iomanip>


double max_offdiag(const arma::mat& A, int& k, int& l);

void Jacobi_rotate(arma::mat& A, arma::mat& R, int k, int l, int& Len);

void Jacobi_eigensolver(const arma::mat& A, double eps, arma::vec& eigenvalues, arma::mat& 
    eigenvectors, const int maxiter, int& iterations, bool& converged);