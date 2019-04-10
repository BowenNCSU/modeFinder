#include <Rcpp.h>
#include "BP.h" // My header
#include "empirical.h" // My header
#include "kernel.h" // My header
#include "set_bound.h" // My header
using namespace Rcpp;
// Create emp_mode class and a method which returns a List and Brent within the class
// Adopt pointer when updating coefficient maps


//' Compute the univariate empirical mode via Bernstein polynomials.
//' 
//' 
//'  
//' 
//' @param data An univariate sample data.
//' @param fix_lower A lower bound of the support; If \code{NULL} (default), set min{0, the minimum of the sample data}; Set a particualr value otherwise.
//' @param fix_upper An upper bound of the support; If \code{NULL} (default), set max{1, the maximum of the sample data}; Set a particualr value otherwise.
//' @param smooth An indicator of if performing Bernstein polynomials smoothing or not; If \code{TRUE} (default), adopt Bernstein polynomials density function; If \code{FALSE}, adopt empirical density function without smoothing.
//' @param smooth_option If \code{smooth = TRUE}, an option of smoothing method for empirical density function; \code{"Bernstein"} (default) indicates Bernstein polynomials density estimate, and \code{"kernel"} the kernel density estimate.
//' @param smooth_density If \code{smooth = TRUE}, an indicator of if output the density function or not; If \code{FALSE} (default), only the mode estimate is given; If \code{TRUE}, the smooth density function is also included in the output.
//' @param density_points If \code{smooth = TRUE}, a number of equally spaced points for the output density function or not, default is 512.
//' @param m_degree If \code{smooth = TRUE} and \code{smooth_option = "Bernstein"}, a degree of Bernstein polynomials; If \code{NA} (default), pick the one with maximal p-value until exceeding 0.95; Set a particualr value otherwise.
//' @param sample_size If \code{smooth = TRUE} and \code{smooth_option = "Bernstein"}, a subsample size used for Anderson-Darling criterion, and default is 1e+03.
//' @param m_knots If \code{smooth = TRUE} and \code{smooth_option = "kernel"}, a number of knots in finding the maximum of empirical density function; If \code{NA} (default), set \code{m_knots} = sample size of \code{data}; Set a particualr value otherwise.
//' @return An empirical mode of \code{data}, a degree of Bernstein polynomials adopted with the maximum p-value and the corresponding maximum p-value of the Anderson-Darling criterion based on the subsample and Bernstein polynomials; \code{NA} if not applicable.
//' @examples
//' n = 1e6
//' set.seed(1)
//' data = rbeta(n, 2, 5)
//' emp_mode(data)
//' @export
// [[Rcpp::export]]
Rcpp::List emp_mode(std::vector<double> data, 
                    Rcpp::Nullable<double> fix_lower = R_NilValue, 
                    Rcpp::Nullable<double> fix_upper = R_NilValue, 
                    bool smooth = true, 
                    String smooth_option = "Bernstein", 
                    bool smooth_density = false,
                    int density_points = 512L,
                    Rcpp::Nullable<int> m_degree = R_NilValue, 
                    int sample_size = 1000L, 
                    Rcpp::Nullable<int> m_knots = R_NilValue){
  double lower_ptr;
  double upper_ptr;
  set_bound(data, fix_lower, fix_upper, & lower_ptr, & upper_ptr);
  
  if(smooth){ // Default: BP smooth
    if(smooth_option == "Bernstein"){
// Load data ---------------------------
      emp_mode_class * new_class = new emp_mode_class(data);
// Set subsample size ---------------------------
      new_class->set_sample_size(sample_size);
// Set degree ---------------------------
      if(m_degree.isNull()){
        new_class->set_do_choose_m(true);
      }else{
        new_class->set_m(Rcpp::as<int>(m_degree));
        new_class->set_do_choose_m(false);
      }
// Set density output
      new_class->set_do_density(smooth_density);
      new_class->set_n_points(density_points);
// Find emprirical mode ---------------------------
      Rcpp::List res = new_class->emp_mode_c(lower_ptr, upper_ptr);
// Return results ---------------------------
        return res;
    }else if(smooth_option == "kernel"){ // Kernel density
      // Set density output
      Rcpp::List res = mode2(data, smooth_density, density_points);
      return res;
    }
  }else{ // Empirical approach
    int knots = m_knots.isNull() ? data.size() : Rcpp::as<int>(m_knots);

    Rcpp::List res = empirical(data, knots, lower_ptr, upper_ptr);
// Return results ---------------------------
    return res;
  }
  return List();
}





