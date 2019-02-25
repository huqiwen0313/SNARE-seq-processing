//cpp function to generate single-cell count matrix - sparse matrix structure
//Qiwen Hu 2018

#include <RcppEigen.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]

Eigen::MatrixXd ccwindow_n_tags_barcodes(NumericVector pos, NumericVector barcode, int nbarcode,
                                        NumericVector wpos, double window_half_size){
// pos - vector contains sorted position of each reads
// barcode - vector contain the indexed barcode information
// nbarcode - number of total barcodes
// wpos - centered peak positions
// window_half_size - half sliding window for counting reads
  int nw = wpos.size();
  int n = pos.size();
  double whs = window_half_size;
  
  //output matrix - sparse matrix
  Eigen::SparseMatrix<double> count(nw, nbarcode);
  
  // current array start/end indecies
  int cs = 0; int ce = 0;
  
  for(int i=0; i<nw; i++){
    // increment ce to add windows that are overlapped
    double ep = wpos[i] + whs;
    while(ep > pos[ce] & ce<n) { ce++; }
    
    //increment cs to drop windows that have already ended
    double sp = wpos[i] - whs;
    while(pos[cs] < sp) { cs++ ; }
    
    for(int j=cs; j<ce; j++){
      double value = count.coeff(i, barcode[j] - 1) + 1;
      count.insert(i, barcode[j] - 1) = value;
    }
  }
  return(Eigen::MatrixXd(count));
}
