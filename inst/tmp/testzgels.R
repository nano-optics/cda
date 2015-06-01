
library(Rcpp)
cppFunction(depends = "RcppArmadillo",
            '
            #define ARMA_NO_DEBUG
            arma::colvec main_call(const arma::mat& A, const arma::colvec& b)
            {
            const int n=A.n_rows;
            arma::colvec x(n);
            x = solve(A, b);
            return x;
            } 
            ')


set.seed(123)
n <- 10
A = matrix(rnorm(n^2), n)
b <- rnorm(n) 

test <- solve(A, b)
test2 <- main_call(A, b)
