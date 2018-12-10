#include <Rcpp.h>
#include "modeFinder_types.h" // My header
using namespace Rcpp;


//' Compute the univariate empirical mode via Bernstein polynomials.
//' 
//' 
//'  
//' 
//' @param data A univariate sample data.
//' @param threshold A number of subsample size used in Anderson-Darling tests, and default is 100L.
//' @param smooth A indicator of if performing Bernstein polynomials smoothing or not; If TRUE (default), adapt Bernstein polynomials smoothing.
//' @param smooth_option A name of method of finding maximum based on the Bernstein polynomials density, and "optimalOfPdf" (default) indicates finding optimum on the smoothing density.
//' @param fix_lower A value of the lower bound of the the sample data; If NULL (default), set min{0, the minimum of the sample data}; Set a particualr value otherwise.
//' @param fix_upper A value of the upper bound of the the sample data; If NULL (default), set max{1, the maximum of the sample data}; Set a particualr value otherwise.
//' @param m_degree A number of degree of Bernstein polynomials; If NA (default), set the one with maximal p-value that less than 0.95; Set a particualr value otherwise.
//' @param m_knots A number of knots in finding the maximum of Bernstein polynomials density; If NA (default), set n / log(n); Set a particualr value otherwise.
//' @return The empirical mode of \code{data}, the degree of Bernstein polynomials used and the corresponding p-value of the Anderson-Darling test based on the subsample and Bernstein polynomials.
//' @examples
//' m_threshold = 1e3
//' n = 1e6
//' set.seed(1)
//' data = rbeta(n, 2, 5)
//' emp_mode(data, threshold = m_threshold)
//' @export
// [[Rcpp::export]]
Rcpp::List emp_mode(std::vector<double> data, 
                    Rcpp::Nullable<int> threshold_ = R_NilValue, 
                             bool smooth = true, 
                             String smooth_option = "optimalOfPdf", 
                             Rcpp::Nullable<double> fix_lower = R_NilValue, 
                             Rcpp::Nullable<double> fix_upper = R_NilValue, 
                             Rcpp::Nullable<int> m_degree = R_NilValue, 
                             Rcpp::Nullable<int> m_knots = R_NilValue)
{
  if(smooth){
    if(smooth_option == "optimalOfPdf"){
// Load data ---------------------------
      emp_mode_class * new_class = new emp_mode_class(data);
// Set subsample size ---------------------------
      int threshold = threshold_.isNull() ? data.size() : Rcpp::as<int>(threshold_) ;
      new_class->set_threshold(threshold);
// Find emprirical mode ---------------------------
      Rcpp::List res = new_class->emp_mode_c(fix_lower.isNull(), 
                                  fix_lower.isNull() ? 0 : Rcpp::as<double>(fix_lower), 
                                  fix_upper.isNull(), 
                                  fix_upper.isNull() ? 0 : Rcpp::as<double>(fix_upper));
// Return results ---------------------------
        return res;
    }
  }
  return List();
}




