// [[Rcpp::depends(RcppEigen)]]
 
#include <RcppEigen.h>
  
// [[Rcpp::export]]
SEXP eigenMatMult(Eigen::MatrixXd A, Eigen::MatrixXd B){
    Eigen::MatrixXd C = A * B;

    return Rcpp::wrap(C);
}


// [[Rcpp::export]]
SEXP eigenMatMult3(Eigen::MatrixXd A, Eigen::MatrixXd B, Eigen::MatrixXd C){
    Eigen::MatrixXd D = A * B * C;
    
    return Rcpp::wrap(D);
}


/* // [[Rcpp::export]]
SEXP eigenMatMult1(Eigen::MatrixXd A ){
    Eigen::MatrixXd C = A.adjoint() * A;

    return Rcpp::wrap(C);
}  */

// [[Rcpp::export]]
SEXP eigenMapMatMult(const Eigen::Map<Eigen::MatrixXd> A, Eigen::Map<Eigen::MatrixXd> B){
    Eigen::MatrixXd C = A * B;

    return Rcpp::wrap(C);
}

// [[Rcpp::export]]
SEXP eigenMapMatMult3(const Eigen::Map<Eigen::MatrixXd> A, Eigen::Map<Eigen::MatrixXd> B, Eigen::Map<Eigen::MatrixXd> C){
    Eigen::MatrixXd D = A * B * C;
    
    return Rcpp::wrap(D);
}


// [[Rcpp::export]]
SEXP eigenMapMatMult4(const Eigen::Map<Eigen::MatrixXd> A, Eigen::Map<Eigen::MatrixXd> B, Eigen::Map<Eigen::MatrixXd> C, Eigen::Map<Eigen::MatrixXd> D){
    Eigen::MatrixXd E = A * B * C * D;
    
    return Rcpp::wrap(E);
}


// [[Rcpp::export]]
SEXP eigenMapMatMult1(const Eigen::Map<Eigen::MatrixXd> A, Eigen::Map<Eigen::MatrixXd> B){
    Eigen::MatrixXd C = A.adjoint() * B;

    return Rcpp::wrap(C);
} 

/* // [[Rcpp::export]]
SEXP eigenMapMatMult1(const Eigen::Map<Eigen::MatrixXd> A){
    Eigen::MatrixXd C =A.adjoint() * A;

    return Rcpp::wrap(C);
}  */

// [[Rcpp::export]]
SEXP eigenMapMatMult2(const Eigen::Map<Eigen::MatrixXd> A){
	  int m(A.cols());
    Eigen::MatrixXd C(m, m);
	C.setZero();
	C+=C.selfadjointView<Eigen::Lower>().rankUpdate(A.adjoint()); 
    
    return Rcpp::wrap(C);
} 
 
 
