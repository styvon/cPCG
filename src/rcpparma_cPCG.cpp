// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#include <omp.h>
#include "RcppArmadillo.h"
#include <math.h> 

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;
using namespace std;
using namespace arma;


// incomplete cholesky factorization
// [[Rcpp::export]]
arma::mat icc(arma::mat A){
  int N = A.n_cols ;
  arma::mat temp = A;
  for(int k = 0; k < N; k++){
    temp(k,k) = sqrt(temp(k,k));
    for(int i = k + 1; i < N; i++){
      if(temp(i,k) != 0){
        temp(i,k) = temp(i,k)/temp(k,k);
      }
    }
    for(int j = k + 1; j < N; j++){
      for(int i= j; i < N; i++){
        if(temp(i,j) != 0){
          temp(i,j) = temp(i,j) - temp(i,k)*temp(j,k);
        }
      }
    }
  }
  
  for(int i = 0; i<N; i++){
    for(int j = i+1; j<N; j++){
      temp(i,j) = 0;
    }
  }
  
  return temp;
}


//' Conjugate gradient method
//'
//' Conjugate gradient method for solving system of linear equations Ax = b,
//' where A is symmetric and positive definite.
//'
//' @title Solve for x in Ax = b using conjugate gradient method.
//' @param A matrix, symmetric and positive definite.
//' @param b vector, with same dimension as number of rows of A.
//' @param tol numeric, threshold for convergence, default is \code{1e-6}.
//' @param maxIter numeric, maximum iteration, default is \code{1000}.
//' @return A vector representing solution x.
//' @examples
//' \dontrun{
//' test_A <- matrix(c(4,1,1,3), ncol = 2)
//' test_b <- matrix(1:2, ncol = 1)
//' cgsolve(test_A, test_b, 1e-6, 1000)
//' }
// [[Rcpp::export]]
arma::vec cgsolve(arma::mat A, arma::vec b, float tol = 1e-6, int maxIter = 1000) {
  /* Function for solving linear equations Ax = b using conjugate gradient
  !!!todo: preconditioning,c++
    Input:
    A: matrix.
  b: vector
  Output
  x: vector
  */
    // get number of columns of A
  int C = A.n_cols ;
  int R = A.n_rows ;
  // initiate solution x as zeros
  arma::vec x(C) ;
  x.zeros() ; 
  arma::vec oneVec(C);
  oneVec.ones() ;
  
  arma::vec r = b - A * x;
  arma::vec p = r;
  double rs_old = as_scalar( r.t() * r );
  double rs_new=1;
  arma::vec rs_ratio(1);
  
  arma::vec Ap(R);
  double alpha;
  // vector version of alpha
  arma::vec alphaVec(1);
  
  for(int iter = 0; (iter < maxIter) && (rs_new > tol); iter++){
    Ap = A * p;
    alpha = rs_old / as_scalar(p.t() * Ap);
    alphaVec.fill(alpha); 
    x += (oneVec * alphaVec) % p;
    r -= (oneVec * alphaVec) % Ap;
    rs_new = as_scalar( r.t() * r );
    rs_ratio.fill(rs_new / rs_old); 
    
    p = r + (oneVec * rs_ratio) % p;
    rs_old = rs_new;
    if (iter >= maxIter){
      Rcout << "cg did not converge." << endl;
    }
  }
  
  return x;
  
} 


//' Preconditioned conjugate gradient method
//'
//' Preconditioned conjugate gradient method for solving system of linear equations Ax = b,
//' where A is symmetric and positive definite.
//'
//' @title Solve for x in Ax = b using preconditioned conjugate gradient method.
//' @param A matrix, symmetric and positive definite.
//' @param b vector, with same dimension as number of rows of A.
//' @param preconditioner string, method for preconditioning: \code{"Jacobi"} (default), \code{"SSOR"}, or \code{"ICC"}.
//' @param tol numeric, threshold for convergence, default is \code{1e-6}.
//' @param maxIter numeric, maximum iteration, default is \code{1000}.
//' @return A vector representing solution x.
//' @examples
//' \dontrun{
//' test_A <- matrix(c(4,1,1,3), ncol = 2)
//' test_b <- matrix(1:2, ncol = 1)
//' pcgsolve(test_A, test_b, "ICC")
//' }
// [[Rcpp::export]]
arma::vec pcgsolve(arma::mat A, arma::vec b, std::string preconditioner = "Jacobi", float tol = 1e-6, int maxIter = 1000) {
  /* Function for solving linear equations Ax = b using preconditioned conjugate gradient
  Input:
  A: matrix.
  b: vector
  preconditioner: string, type of preconditioner
  Output
  x: vector
  */
  // get number of columns of A
  int C = A.n_cols ;
  int R = A.n_rows ;
  
  // get preconditioner M
  arma::mat M;
  if (preconditioner == "Jacobi"){
    M = arma::diagmat(A);
  } else if(preconditioner == "SSOR"){
    arma::mat D = arma::diagmat(A);
    arma::mat L = arma::trimatl(A);
    M = (D+L) * D.i() * (D+L).t();
  } else if(preconditioner == "ICC"){
    M = icc(A);
  }

  
  // initiate solution x as zeros
  arma::vec x(C) ;
  x.zeros() ; 
  
  arma::vec oneVec(C);
  oneVec.ones() ;
  
  arma::vec r = b - A * x;
  arma::mat Minv = M.i();
  arma::vec z = Minv * r;
  arma::vec p = z;
  double rz_old = as_scalar( r.t() * z );
  double rz_new=1;
  arma::vec rz_ratio(1);
  
  arma::vec Ap(R);
  double alpha;
  // vector version of alpha
  arma::vec alphaVec(1);
  
  for(int iter = 0; (iter < maxIter) && (rz_new > tol); iter++){
    Ap = A * p;
    alpha = rz_old / as_scalar(p.t() * Ap);
    alphaVec.fill(alpha); 
    x += (oneVec * alphaVec) % p;
    r -= (oneVec * alphaVec) % Ap;
    z = Minv * r;
    rz_new = as_scalar( z.t() * r );
    rz_ratio.fill(rz_new / rz_old); 
    
    p = z + (oneVec * rz_ratio) % p;
    rz_old = rz_new;
    if (iter >= maxIter){
      Rcout << "pcg did not converge." << endl;
    }
  }
  
  return x;
  
} 


void matvecmult_sparse(const sp_mat& A, const vec *x, vec *outvec, int nrows){
  int i, ID, nt;
  double rowsum;

  // omp_set_num_threads(nthreads);


  #pragma omp parallel 
  {
    // nt=omp_get_num_threads();
    // printf("OMP with %d threads\n", nt);
    #pragma omp for private(i,rowsum,ID) nowait
    for(i=0; i<nrows; i++){
      rowsum = 0.0;
      arma::sp_mat::const_row_iterator it = A.begin_row(i);
      arma::sp_mat::const_row_iterator it_end = A.end_row(i);
      while(it != it_end){
        rowsum += (*it) * (*x)( it.col() )  ;
        // cout << (*it) << endl;
        ++it;
      }
      (*outvec)(i) = rowsum;

    // ID=omp_get_thread_num();
    // if(i% 30 ==0)
    //  printf("i: %d, OMPthread(%d)\n", i, ID);
    }
  }

}


// test output
vec getMatvecmult_sparse(const sp_mat& A, const vec& xvec){
  int nrows = A.n_rows;

  const vec* px = &xvec;
  vec outvec(nrows);
  vec *poutvec = &outvec;

  matvecmult_sparse(A, px, poutvec, nrows);

  return outvec;
} 





//' Calculation of sparse matrix inverse and vector product using conjugate gradient descent
//'
//' This Rcpp function utilizes OpenMP and calculates inverse of a sparse matrix and product with a vector 
//' @param A          A (n x n) matrix.
//' @param b          A n-vector.
//' @param tol        Tolerance for conjugate gradient algorithm.
//' @param maxIter    Max number of iterations for conjugate gradient algorithm.
//' @param nThreads   Number of threads for OMP.
//' @param nSubThreads   Number of threads for OMP nested tasks.
//' @return      A n-vector of A_inv b
// [[Rcpp::export]]
arma::vec cgsolve_sparseOMP(const arma::sp_mat& A, const arma::vec& b, float tol= 1e-6, int maxIter = 1000, int nThreads=1, int nSubThreads=1){
  // get number of columns of A
  int C = A.n_cols ;
  int R = A.n_rows ;

  // initiate solution x as zeros
  arma::vec x(C) ;
  x.zeros() ; 
  arma::vec oneVec(C);
  oneVec.ones() ;

  arma::vec r = b - A * x; 
  arma::vec p = r;
  double rs_old = as_scalar( r.t() * r );
  double rs_new=1;
  arma::vec rs_ratio(1);

  arma::vec Ap(R);
  double pAp;
  double alpha;
  // vector version of alpha
  arma::vec alphaVec(1);

  //pointers
  arma::vec* p_Ap = &Ap;   //pointer to Ap
  const arma::vec* p_p = &p; //pointer to p

  int iter=0, abort=0;

  omp_set_num_threads(nThreads);
  omp_set_nested(true);
  
#pragma omp parallel default(shared) private(iter)
{
  while(!abort){  
    matvecmult_sparse(A, p_p , p_Ap, R); // get Ap


#pragma omp single
    pAp = 0.0;

#pragma omp parallel for reduction(+:pAp) num_threads(nSubThreads)
    for(int i=0; i <  R; i++)
      pAp += p(i)* (*p_Ap)(i) ;

#pragma omp single
{
    alpha = rs_old / pAp;
    alphaVec.fill(alpha); 
    x += (oneVec * alphaVec) % p;
      r -= (oneVec * alphaVec) % Ap;

      rs_new = 0.0;
}

#pragma omp parallel for reduction(+:rs_new) num_threads(nSubThreads)
    for(int i=0; i <  R; i++)
      rs_new += r(i)* r(i) ;

#pragma omp single
{
    rs_ratio.fill(rs_new / rs_old);
    p = r + (oneVec * rs_ratio) % p;
    rs_old = rs_new;
    if(rs_new <= tol)
      abort=1;
}

  }// end of while(!abort)
  iter++;
}//end of omp parallel outer
  return x;
}




///*** R
// #vecA = rnorm(100000)
// #sumsq_parallel(vecA,3)


// library('Matrix')
// #n <- 1000
// #M<- Matrix(0, nrow = n, ncol = n, sparse = TRUE)
// #M <- M + diag(n)
// M <- Matrix(c(4,1,1,3), nrow=2, ncol=2, sparse = TRUE)
// x <- 1:2
// #x <- 1:1000
// #out <- getMatvecmult_sparse(M,x, 3)
// out <- cgsolve_sparseOMP(M,x, nThreads=3, nSubThreads=3)
// out
//*/