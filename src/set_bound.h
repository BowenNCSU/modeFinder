#include <vector>
#include <Rcpp.h>
using namespace std;

// Calcualte the left and right boundary for the "emp_mode"
//
//

void set_bound(std::vector<double> & data, 
          Rcpp::Nullable<double> fix_lower, 
          Rcpp::Nullable<double> fix_upper, 
          double * lower,
          double * upper){
std::sort(data.begin(), data.end());
double min_data = *data.begin();
double max_data = data.back();
double delta = 1 / (double) data.size(); 

if(fix_lower.isNull()){
  *lower = min_data < 0 ? min_data - delta : 0;
}else{
  *lower = Rcpp::as<double>(fix_lower);
}

if(fix_upper.isNull()){
  *upper = max_data > 1 ? max_data + delta : 1;
}else{
  *upper = Rcpp::as<double>(fix_upper);
}

}
