// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <iostream>

using namespace Rcpp ;
using namespace RcppArmadillo ;
using namespace arma;
using namespace std;

//' @title Generate a random sample of points on the unit sphere
//' @description Random sample
//' @describeIn  sample_random random sample
//' @return 3xN matrix
//' @family sample
//' @examples 
//' sample_random(10)
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

//' @title Generate a sample of points on the unit sphere
//' @description Random sample with minimum exlusion zone ("hard-core process")
//' @param exclusion minimum distance allowed between points
//' @param maxiter maximum number of iterations
//' @param k number of extra new points to try at each iteration
//' @return 3xN matrix
//' @describeIn sample_random random sample with exclusion zone
//' @family sample
//' @examples 
//' sample_hc(10)
//' @export
// [[Rcpp::export]]
arma::mat sample_hc(const int N,
  const double exclusion=0.1,
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
      return(s.cols(pick(span(0,N-1))));
    }

    // otherwise, replace bad elements and try again
    uvec bad;
    int p;
    int iter =0;
    while(todo){

      bad = find(indices == 0);
      p = bad.n_elem;
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
        return(s.cols(pick(span(0,N-1))));
      }
      if(iter >= maxiter){
        // should warn somewhow
        Rcpp::warning("max iter reached");
        break;
      }
      todo = (pick.n_elem < N) && (iter < maxiter);
      iter++;
    }
    // if we get to here, this won't be a good sample
    // but for consistency always return N points regardless
    return(s.cols(span(0,N-1)));
  }
