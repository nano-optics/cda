// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <iostream>

using namespace Rcpp ;
using namespace RcppArmadillo ;
using namespace arma;
using namespace std;

/*
GOAL
"hard-core" process: generate N points on a sphere no closer than exclusion zone


STRATEGY
- generate N + k points at random, store col-wise in matrix s
- test all pairwise dists, store results in bool vector indices
- if more than N good points, return first N points
- if n missing, generate n + k new random points and iterate
- terminate with warning if more than nmax iterations,
- always return N points
*/


//' @export
// [[Rcpp::export]]
arma::mat sample_random(const int N){
  arma::rowvec phi = arma::rowvec(N), theta=phi;
  arma::mat result = arma::mat(3, N);
  phi = datum::pi*(2*arma::randu<arma::rowvec>(N) - 1.0);
  theta = acos(2*arma::randu<arma::rowvec>(N) - 1.0);
  result.row(0) = sin(theta)%cos(phi);
  result.row(1) = sin(theta)%sin(phi);
  result.row(2) = cos(theta);
  return(result);
}

//' @export
// [[Rcpp::export]]
arma::mat sample_hc(const int N,
  const double exclusion,
  const int maxiter=200,
  const int k=30){

    arma::mat result = arma::mat(N,3);
    arma::mat s = arma::mat(N+k,3);
    arma::uvec indices = arma::uvec(N+k);

    // initial step: random points
    s = sample_random(N+k);
    indices.ones(); // assume good until shown otherwise

    int ii, jj;
    double dist;
    bool test;
    for(ii=0;ii<N+k;ii++){// loop over all points
      for(jj=ii+1;jj<N+k;jj++){//loop over all other points
        dist = norm(s.col(ii) - s.col(jj), 2);
        test = (dist < exclusion);
        if(test) {
          indices(ii) = 0;
          break; // this ii point was actually bad
        }
      }
    }
    bool todo = (sum(indices) < N);
    // if more than N, we're done
    uvec pick = find(indices);
    if(!todo){
      // std::cout << "here" << std::endl;
      return(s.cols(pick(span(0,N-1))));
    }

    // otherwise, replace bad elements and try again
    uvec bad;
    int p;
    int iter =0;
    while(todo){

      bad = find(indices == 0);
      p = bad.n_elem;
      // std::cout << "iter: " << iter << ", p: " << p << std::endl;
      s.cols(bad) = sample_random(p);
      for(ii=0; ii<p; ii++){ // points to test
        indices(bad(ii)) = 1; // assume it is good until shown otherwise
        for(jj=0; jj<N+k; jj++){// loop over all other points
          if(bad(ii) == jj) continue;
          dist = norm(s.col(bad(ii)) - s.col(jj), 2);
          test = (dist < exclusion);
          if(test) {
            indices(bad(ii)) = 0;
            break; // this ii point was actually bad
          }
        }
      }
      // if more than N, we're done
      pick = find(indices);
      if(pick.n_elem >= N){
        // std::cout << "here2" << std::endl;
        return(s.cols(pick(span(0,N-1))));
      }
      if(iter >= maxiter){
        // should warn somewhow
        Rcpp::Rcout << "warning: max iter reached" << std::endl;
        break;
      }
      todo = (pick.n_elem < N) && (iter < maxiter);
      iter++;
    }
    // if we get to here, this won't be a good sample
    // but for consistency always return N points regardless
    return(s.cols(span(0,N-1)));
  }
