#include <vector>
#include <Rcpp.h>
using namespace std;

Rcpp::Environment base("package:stats"); 
Rcpp::Function density_r = base["density"];

Rcpp::List mode2(std::vector<double> & data, 
                 bool & smooth_density,
                 int & density_points){
  
  Rcpp::List fit = density_r(Rcpp::_["x"] = data, Rcpp::_["n"] = density_points); // Default setting
  double max_y = Rcpp::as<Rcpp::NumericVector>(fit["y"])[0];
  double max_x = Rcpp::as<Rcpp::NumericVector>(fit["x"])[0];
  int data_size = data.size();
  for(int i = 0; i < data_size; i++){
    if(Rcpp::as<Rcpp::NumericVector>(fit["y"])[i] > max_y){
      max_y = Rcpp::as<Rcpp::NumericVector>(fit["y"])[i];
      max_x = Rcpp::as<Rcpp::NumericVector>(fit["x"])[i];
    }
  }
  Rcpp::List res;
  res["mode"] = max_x;
  if(smooth_density){
    res["y"] = Rcpp::as<Rcpp::NumericVector>(fit["y"]);
    res["x"] = Rcpp::as<Rcpp::NumericVector>(fit["x"]);
  }
  return res;
}






