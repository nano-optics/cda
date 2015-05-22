set.seed(123)
n <- 1000
A = matrix(rnorm(n^2) + rnorm(n^2)*1i, n)
y <- rnorm(n) + rnorm(n)*1i

norm <- function(x){
  c(Re(crossprod(Conj(x), x)))
}

## solve Ax=y iteratively
cg <- function(A, y, x_0 = NULL, tol=1e-3, nmax=20,
               variant = FALSE, verbose = TRUE){
  
  n <- 0
  rel.error <- Inf
  
  ## init
  if(is.null(x_0))
    x_0 <- 0*y
  B <- t(Conj(A))
  z <- B %*% y
  g_i <- z - B %*% A %*% x_0
  p_i <- g_i
  w_i <- A %*% x_0
  v_i <- A %*% p_i
  x_i <- x_0
  while(n < nmax && rel.error > tol){
    
    alpha_i <- norm(g_i) / norm(v_i)
    #print(alpha_i)
    x_ip1 <- x_i + alpha_i * p_i
    w_ip1 <- w_i + alpha_i * v_i
    g_ip1 <- z - B %*% w_ip1
    beta_i <- norm(g_ip1) / norm(g_i)
    # updating
    p_i <- g_ip1 + beta_i * p_i
    if(variant) v_i <- A %*% g_ip1 + beta_i * v_i else
      v_i <- A %*% p_i
    rel.error <- norm(x_ip1 - x_i) / norm(x_ip1) 
    x_i <- x_ip1 
    if(variant) w_i <- w_ip1 else
      w_i <- A %*% x_i
    g_i <- g_ip1 
    
    if(verbose)
      message("n: ", n, ", rel.error: ", rel.error)
    n <- n + 1
    
  }
  
  c(x_i)
  
}


library(Rcpp)
cppFunction(depends = "RcppArmadillo",
            'arma::cx_colvec cg_solve(const arma::cx_mat& A, 
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
            
            
            } ')


set.seed(123)
n <- 2000
A = matrix(rnorm(n^2) + rnorm(n^2)*1i, n)
y <- rnorm(n) + rnorm(n)*1i

test <- solve(A, y)
test1 <- cg(A, y, nmax=10, verbose = FALSE)
test2 <- cg_solve(A, y, 0*y, 10, 1e-2)
test2
max(Mod(test1 - test2))
# 
# microbenchmark(blind = cg_solve(A, y, 0*y, 10, 1e-2),
#                 guess = cg_solve(A, y, test, 10, 1e-2))

# Unit: milliseconds
# expr      min       lq     mean   median       uq      max neval cld
# blind 445.8212 457.1995 464.8376 464.0516 471.3204 497.2138   100   b
# guess 246.6692 253.1206 257.0085 256.8858 259.9855 270.0985   100  a 
library(microbenchmark)
# 
# microbenchmark(bf = solve(A, y),
#                cg = cg(A, y, nmax=10, verbose = FALSE),
#                cpp = cg_solve(A, y, 0*y, 10, 1e-2))
# 
# Unit: milliseconds
# expr       min        lq      mean    median        uq       max neval cld
# bf  617.5702  631.0103  654.0801  639.8153  653.4558  771.4092   100  b 
# cg 1948.7679 1977.7604 2014.2302 1993.6667 2026.0109 2140.7054   100   c
# cpp  107.7403  113.9244  116.4186  116.3641  118.9065  125.7930   100 a  
