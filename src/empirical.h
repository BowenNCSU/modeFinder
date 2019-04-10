#include <vector>
#include <Rcpp.h>
using namespace std;

Rcpp::List empirical(std::vector<double> & data, int & knots, double & lower, double & upper){
  Rcpp::List res;
  int n = data.size();
  double max_diff = -1;
  int max_left_bound = -1;
  int cum_index = 0;
  double criterion = 0 * (upper - lower) + lower; 
  while(cum_index < n && data[cum_index] <= criterion){
    cum_index ++ ;
  }
  double Fn_value_0 = cum_index / (double) n;
  
  for(int i = 0; i < knots; i++) {
    criterion = (i + 1) / (double) knots * (upper - lower) + lower;
    while(cum_index < n && data[cum_index] <= criterion){
      cum_index ++ ;
    }
    double Fn_value = cum_index / (double) n;
    double diff = Fn_value - Fn_value_0;
    
    if(diff > max_diff){
      max_left_bound = i;
      max_diff = diff;
    }
    Fn_value_0 = Fn_value;
  }
  
  res["number of knots"] = knots;
  res["mode"] = (2 * max_left_bound + 1) / (double) knots / 2.0 * (upper - lower) + lower;
  Rcpp::Rcout << "The mode estimate is choosen as the mid point of interval with maximum difference among knots equally-spaced intervals." << std::endl;
  
  return(res);
}



