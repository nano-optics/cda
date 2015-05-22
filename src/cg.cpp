#include <RcppArmadillo.h>
#include <iostream>

using namespace Rcpp ;
using namespace RcppArmadillo ;
using namespace arma;
using namespace std;


//
// Solves a linear system iteratively via conjugate-gradient
//
// A is a complex matrix
// y is the right-hand side
// x0 is the initial guess
// nmax is the max number of iterations
// tol is the tolerance
arma::cx_colvec cg_solve(const arma::cx_mat& A, 
			      const arma::cx_vec& y,	
            const arma::cx_vec& x0,		
            const double nmax, 	
            const double tol)
            {
            
            int n = 0; // counter
            double rel_error = 1e4; // large value
            const int N = y.n_elem;
            
            // temporary variables
            double alpha_i, beta_i;
            const arma::cx_mat B = A.t(); // Hermitian transpose of A
            arma::cx_colvec z(N), g_i(N), p_i(N), w_i(N), v_i(N), x_i(N);   
            arma::cx_colvec g_ip1(N), p_ip1(N), w_ip1(N), v_ip1(N), x_ip1(N);   
            
            // initial step

            z = B * y;
            g_i = z - B * A * x0;
            p_i = g_i;
            w_i = A * x0;
            v_i = A * p_i;
            x_i = x0;

            while ((n < nmax) && (rel_error > tol)){ 
            
            alpha_i = real(cdot(g_i, g_i) / cdot(v_i, v_i));
            // std::cout << alpha_i << std::endl;
            x_ip1 = x_i + alpha_i * p_i;
            w_ip1 = w_i + alpha_i * v_i;
            g_ip1 = z - B * w_ip1;
            beta_i = real(cdot(g_ip1,g_ip1) / cdot(g_i,g_i));
            // updating
            p_i = g_ip1 + beta_i * p_i;
            v_i = A * p_i;
            rel_error = norm(x_ip1 - x_i) / norm(x_ip1) ;
            x_i = x_ip1;
            w_i = A * x_i;
            g_i = g_ip1;
            n++;
            }; // terminate while
            
            return x_i;
            
            
            } 

