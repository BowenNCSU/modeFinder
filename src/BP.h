#include <vector>
#include <map>
#include <Rcpp.h>
#include "lbfgsb.cpp"
#include "AnDarl.cpp"
using namespace std;



class emp_mode_class {
public:
  emp_mode_class(std::vector<double> & data_){
    data = data_;
    n = data_.size();
    // std::sort(data.begin(), data.end());
    subsample_size = 1000;
    do_choose_m = true;
    do_density = false;
    n_points = 512;
  }


  
  void set_m(int m_){
    m = m_;
  }
  
  void set_do_choose_m(bool TorF){
    do_choose_m = TorF;
  }
  
  void set_do_density(bool TorF){
    do_density = TorF;
  }  
  
  void set_n_points(int np){
    n_points = np;
  }
  
  void set_sample_size(int t_){
    subsample_size = t_;
  }
  
  double mix_cdf_c(double u) {
    if(u < 0){
      return 0;
    }else if(u >= 1){
      return 1;
    }else if(u == 0){
      int cum_index = 0;
      while(cum_index < n && data[cum_index] <= 0){
        cum_index ++ ;
      }
      double Fn_value = cum_index / (double) n;
      return Fn_value;
    }else if ( coef_map_cdf_bern.find(m) == coef_map_cdf_bern.end() ) {
      // not found
      // coef.vec.cdf.bern <- Fn(0 : m / m) * dbeta(u, 0 : m + 1, m - (0 : m) + 1) / (m + 1)
      // coef.vec.cdf.bern <- log(Fn(0 : m / m) / (m + 1)) + dbeta_r(u, 0 : m + 1, m - (0 : m) + 1, log = TRUE)
      std::vector<double> coef_vec_cdf_bern(m + 1);
      std::vector<double> Fn_vec_cdf_bern(m + 1);
      double summation = 0;
      int cum_index = 0;
      double dbeta_log = m * log(1 - u);
      for(int i = 0; i < (m + 1); i++) {
        // double dbeta_log = R::dbeta(u, i + 1, m - i + 1, 1);
        if(i > 0){
          dbeta_log += log((m - i + 1) / (double) i * u / (1 - u));
        }
        double criterion = i / (double) m;
        while(cum_index < n && data[cum_index] <= criterion){
          cum_index ++ ;
        }
        double Fn_value = cum_index / (double) n; // / (double) (m + 1);
        Fn_vec_cdf_bern[i] = Fn_value;
        coef_vec_cdf_bern[i] = log(Fn_value) + dbeta_log;
        summation += Fn_value * exp(dbeta_log);
      }
      // m.vec.cdf.bern <<- c(m.vec.cdf.bern, m)
      coef_map_cdf_bern[m] = coef_vec_cdf_bern;
      Fn_map_cdf_bern[m] = Fn_vec_cdf_bern;
      u_map_cdf_bern[m] = u;
      // return(sum(Fn(0 : m / m) * dbeta(u, 0 : m + 1, m - (0 : m) + 1)) / (m + 1)) # 1 / (m+2) 
      return summation;
    } else {
      // found
      double & u_post = u_map_cdf_bern[m];
      std::vector<double> & coef_vec_cdf_bern = coef_map_cdf_bern[m];
      // sum.cdf.bern <- crossprod(coef.ls.cdf.bern[[which(m.vec.cdf.bern == m)]], 
      // (u / u.post) ^ (0 : m) * ((1 - u) / (1 - u.post)) ^ (m : 0))
      double summation = 0;
      for(int i = 0; i < (m + 1); i++){
        summation += exp(coef_vec_cdf_bern[i] + 
          log(u / u_post) * i + log((1 - u) / (1 - u_post)) * (m - i));
      }

      if(isnan(summation) || isinf(summation) || summation < 0 || summation > 1){
        // coef.vec.cdf.bern <- Fn(0 : m / m) * dbeta(u, 0 : m + 1, m - (0 : m) + 1) / (m + 1)
        // coef.vec.cdf.bern <- log(Fn(0 : m / m) / (m + 1)) + dbeta(u, 0 : m + 1, m - (0 : m) + 1, log = TRUE)
        // u.vec.cdf.bern[which(m.vec.cdf.bern == m)] <<- u
        // coef.ls.cdf.bern[[which(m.vec.cdf.bern == m)]] <<- coef.vec.cdf.bern
        // return(sum(exp(coef.vec.cdf.bern)))
        summation = 0;
        std::vector<double> & Fn_vec_cdf_bern = Fn_map_cdf_bern[m];
        double dbeta_log = m * log(1 - u);
        for(int i = 0; i < (m + 1); i++){
          // double dbeta_log = R::dbeta(u, i + 1, m - i + 1, 1);
          if(i > 0){
            dbeta_log += log((m - i + 1) / (double) i * u / (1 - u));
          }
          double Fn_value = Fn_vec_cdf_bern[i];
          coef_vec_cdf_bern[i] = log(Fn_value) + dbeta_log;
          summation += Fn_value * exp(dbeta_log);
        }        
        // coef_map_cdf_bern[m] = coef_vec_cdf_bern;
        u_post = u;
        // return(sum(Fn(0 : m / m) * dbeta(u, 0 : m + 1, m - (0 : m) + 1)) / (m + 1)) # 1 / (m+2) 
        return summation;          
        
      }else{
        // return(sum.cdf.bern)
        return summation;
      }
    }
  }
  
  std::vector<double> cdf_bern_c(std::vector<double> & u){
    int n = u.size();
    std::vector<double> out(n);
    for(int i = 0; i < n; i++){
      out[i] = mix_cdf_c(u[i]);
    }
    return out;
  }
  
  
  double mix_pdf_c(double u) {
    if(u * (1 - u) < 0){
      return 0;
    }else if (u == 0){
      int cum_index = 0;
      while(cum_index < n && data[cum_index] <= 0){
        cum_index ++ ;
      }
      double Fn_value_0 = cum_index / (double) n;
      while(cum_index < n && data[cum_index] <= 1 / (double) m){
        cum_index ++ ;
      }
      double Fn_value_1 = cum_index / (double) n;
      return (Fn_value_1 - Fn_value_0) * m;
    }else if (u == 1){
      int cum_index = n;
      while(cum_index > 0 && data[cum_index - 1] > 1 - 1 / (double) m){
        cum_index -- ;
      }
      double Fn_value = cum_index / (double) n;
      return (1 - Fn_value) * m;
    }else if ( coef_map_pdf_bern.find(m) == coef_map_pdf_bern.end() ) {
      // not found
      // log(Fn((1 : m) / m) - Fn((0 : (m - 1)) / m)) + dbeta(u, 1 : m, m - (1 : m) + 1, log = TRUE)
      std::vector<double> coef_vec_pdf_bern(m);
      std::vector<double> Fn_vec_pdf_bern(m);
      double summation = 0;
      int cum_index = 0;
      double criterion = 0; // m
      while(cum_index < n && data[cum_index] <= criterion){
        cum_index ++ ;
      }
      double Fn_value_0 = cum_index / (double) n;
      double dbeta_log = log(m) + (m - 1) * log(1 - u);
      
      for(int i = 0; i < m; i++) {
        // double dbeta_log = R::dbeta(u, i + 1, m - i, 1);
        if(i > 0){
          dbeta_log += log((m - i) / (double) i * u / (1 - u));
        }
        criterion = (i + 1) / (double) m;
        while(cum_index < n && data[cum_index] <= criterion){
          cum_index ++ ;
        }
        double Fn_value = cum_index / (double) n;
        Fn_vec_pdf_bern[i] = Fn_value - Fn_value_0;
        coef_vec_pdf_bern[i] = log(Fn_value - Fn_value_0) + dbeta_log;
        summation += (Fn_value - Fn_value_0) * exp(dbeta_log);
        Fn_value_0 = Fn_value;
      }
      // m.vec.pdf.bern <<- c(m.vec.pdf.bern, m)
      coef_map_pdf_bern[m] = coef_vec_pdf_bern;
      Fn_map_pdf_bern[m] = Fn_vec_pdf_bern;
      u_map_pdf_bern[m] = u;
      // return(sum(Fn(0 : m / m) * dbeta(u, 0 : m + 1, m - (0 : m) + 1)) / (m + 1)) # 1 / (m+2) 
      return summation;
    } else {
      // found
      double & u_post = u_map_pdf_bern[m];
      std::vector<double> & coef_vec_pdf_bern = coef_map_pdf_bern[m];
      // sum.pdf.bern <- crossprod(coef.ls.pdf.bern[[which(m.vec.pdf.bern == m)]], 
      // (u / u.post) ^ (0 : m) * ((1 - u) / (1 - u.post)) ^ (m : 0))
      double summation = 0;
      for(int i = 0; i < m; i++){
        summation += exp(coef_vec_pdf_bern[i] + 
          log(u / u_post) * i + log((1 - u) / (1 - u_post)) * (m - 1 - i));
      }
      
      if(isnan(summation) || summation < 0){
        // coef.vec.pdf.bern <- Fn(0 : m / m) * dbeta(u, 0 : m + 1, m - (0 : m) + 1) / (m + 1)
        // coef.vec.pdf.bern <- log(Fn(0 : m / m) / (m + 1)) + dbeta(u, 0 : m + 1, m - (0 : m) + 1, log = TRUE)
        // u.vec.pdf.bern[which(m.vec.pdf.bern == m)] <<- u
        // coef.ls.pdf.bern[[which(m.vec.pdf.bern == m)]] <<- coef.vec.pdf.bern
        // return(sum(exp(coef.vec.pdf.bern)))
        summation = 0;
        std::vector<double> & Fn_vec_pdf_bern = Fn_map_pdf_bern[m];
        double dbeta_log = log(m) + (m - 1) * log(1 - u);
        
        for(int i = 0; i < m; i++){
          // double dbeta_log = R::dbeta(u, i + 1, m - i, 1);
          if(i > 0){
            dbeta_log += log((m - i) / (double) i * u / (1 - u));
          }
          double Fn_value = Fn_vec_pdf_bern[i];
          coef_vec_pdf_bern[i] = log(Fn_value) + dbeta_log;
          summation += Fn_value * exp(dbeta_log);
        }        
        // coef_map_pdf_bern[m] = coef_vec_pdf_bern;
        u_post = u;
        // return(sum(Fn(0 : m / m) * dbeta(u, 0 : m + 1, m - (0 : m) + 1)) / (m + 1)) # 1 / (m+2) 
        return summation;          
        
      }else{
        // return(sum.pdf.bern)
        return summation;
      }
    }
  }
  
  std::vector<double> pdf_bern_c(std::vector<double> & u){
    int n = u.size();
    std::vector<double> out(n);
    for(int i = 0; i < n; i++){
      out[i] = mix_pdf_c(u[i]);
    }
    return out;
  }
  
  
  double mix_dpdf_c(double u) {
      // Obtain environment containing function
      // Rcpp::Environment base("package:stats"); 
      // Make function callable from C++
      // Rcpp::Function dbeta_r = base["dbeta"];   
      // Rcpp::Function sum_r = base["sum"];   
      
      if(u * (1 - u) < 0){
        return 0;
      }else if (u == 0){
        int cum_index = 0;
        while(cum_index < n && data[cum_index] <= 0){
          cum_index ++ ;
        }
        double Fn_value_0 = cum_index / (double) n;
        while(cum_index < n && data[cum_index] <= 1 / (double) m){
          cum_index ++ ;
        }
        double Fn_value_1 = cum_index / (double) n;
        while(cum_index < n && data[cum_index] <= 2 / (double) m){
          cum_index ++ ;
        }
        double Fn_value_2 = cum_index / (double) n;
        return (Fn_value_2 - 2 * Fn_value_1 + Fn_value_0) * m * (m - 1);
      }else if (u == 1){
        int cum_index = n;
        while(cum_index > 0 && data[cum_index - 1] > 1 - 1 / (double) m){
          cum_index -- ;
        }
        double Fn_value = cum_index / (double) n;
        while(cum_index > 0 && data[cum_index - 1] > 1 - 2 / (double) m){
          cum_index -- ;
        }
        double Fn_value_ = cum_index / (double) n;
        return (1 - 2 * Fn_value + Fn_value_) * m * (m - 1);
      }else if ( coef_map_dpdf_bern.find(m) == coef_map_dpdf_bern.end() ) {
        // not found
        // (Fn(2 : m / m) - 2 * Fn((2 : m - 1) / m) + Fn((2 : m - 2) / m)) * m * dbeta(u, 2 : m - 1, m - (2 : m) + 1)
        std::vector<double> coef_vec_dpdf_bern(m - 1);
        std::vector<double> Fn_vec_dpdf_bern(m - 1);
        std::vector<bool> sign_vec_dpdf_bern(m - 1);
        double summation = 0;
        int cum_index = 0;
        double criterion = 0; // m
        while(cum_index < n && data[cum_index] <= criterion){
          cum_index ++ ;
        }
        double Fn_value_0 = cum_index / (double) n;
        
        criterion = 1 / (double) m;
        while(cum_index < n && data[cum_index] <= criterion){
          cum_index ++ ;
        }
        double Fn_value_1 = cum_index / (double) n;
        double dbeta_log = log(m * (m - 1)) + (m - 2) * log(1 - u);
        
        for(int i = 0; i < (m - 1); i++) {
          // double dbeta_log = R::dbeta(u, i + 1, m - i - 1, 1);
          if(i > 0){
            dbeta_log += log((m - i - 1) / (double) i * u / (1 - u));
          }
          criterion = (i + 2) / (double) m;
          while(cum_index < n && data[cum_index] <= criterion){
            cum_index ++ ;
          }
          double Fn_value = cum_index / (double) n;
          Fn_vec_dpdf_bern[i] = (Fn_value - 2 * Fn_value_1 + Fn_value_0); // * m;
          sign_vec_dpdf_bern[i] = Fn_value - 2 * Fn_value_1 + Fn_value_0 > 0;
          coef_vec_dpdf_bern[i] = log(std::abs(Fn_value - 2 * Fn_value_1 + Fn_value_0)) + dbeta_log; // * m;
          summation += (Fn_value - 2 * Fn_value_1 + Fn_value_0) * exp(dbeta_log); // * m;
          Fn_value_0 = Fn_value_1;
          Fn_value_1 = Fn_value;
        }
        // m.vec.dpdf.bern <<- c(m.vec.dpdf.bern, m)
        coef_map_dpdf_bern[m] = coef_vec_dpdf_bern;
        sign_map_dpdf_bern[m] = sign_vec_dpdf_bern;
        Fn_map_dpdf_bern[m] = Fn_vec_dpdf_bern;
        u_map_dpdf_bern[m] = u;
        // return(sum(Fn(0 : m / m) * dbeta(u, 0 : m + 1, m - (0 : m) + 1)) / (m + 1)) # 1 / (m+2) 
        return summation;
      } else {
        // found
        double & u_post = u_map_dpdf_bern[m];
        std::vector<double> & coef_vec_dpdf_bern = coef_map_dpdf_bern[m];
        std::vector<bool> & sign_vec_dpdf_bern = sign_map_dpdf_bern[m];
        // sum.dpdf.bern <- crossprod(coef.ls.dpdf.bern[[which(m.vec.dpdf.bern == m)]], 
        // (u / u.post) ^ (0 : m) * ((1 - u) / (1 - u.post)) ^ (m : 0))
        double summation = 0;
        for(int i = 0; i < (m - 1); i++){
          summation += (2 * sign_vec_dpdf_bern[i] - 1) * exp(coef_vec_dpdf_bern[i] + 
            log(u / u_post) * i + log((1 - u) / (1 - u_post)) * (m - 2 - i));
        }
        
        if(isnan(summation)){
          summation = 0;
          std::vector<double> & Fn_vec_dpdf_bern = Fn_map_dpdf_bern[m];
          double dbeta_log = log(m * (m - 1)) + (m - 2) * log(1 - u);
          
          for(int i = 0; i < (m - 1); i++){
            // double dbeta_log = R::dbeta(u, i + 1, m - i - 1, 1);
            if(i > 0){
              dbeta_log += log((m - i - 1) / (double) i * u / (1 - u));
            }
            double Fn_value = Fn_vec_dpdf_bern[i];
            coef_vec_dpdf_bern[i] = log(std::abs(Fn_value)) + dbeta_log;
            summation += Fn_value * exp(dbeta_log);
          }        
          // coef_map_dpdf_bern[m] = coef_vec_dpdf_bern;
          u_post = u;
          // return(sum(Fn(0 : m / m) * dbeta(u, 0 : m + 1, m - (0 : m) + 1)) / (m + 1)) # 1 / (m+2) 
          return summation;          
          
        }else{
          // return(sum.dpdf.bern)
          return summation;
        }
      }
    }
  
  
  std::vector<double> dpdf_bern_c(std::vector<double> & u){
    int n = u.size();
    std::vector<double> out(n);
    for(int i = 0; i < n; i++){
      out[i] = mix_dpdf_c(u[i]);
    }
    return out;
  }
  
  
  double fminfn(double * x){return -mix_pdf_c(*x); }
  
  void fmingr(double * x, double * res){res[0] = -mix_dpdf_c(*x); }
  
  
  //////////////////////////////
  // https://github.com/Microsoft/microsoft-r-open/blob/master/source/src/appl/optim.c
  double * vect(int n)
  {
    return (double *)R_alloc(n, sizeof(double));
  }
  
  void lbfgsb(double *x, double *l, double *u, int *nbd, char *msg, int *fncount, int *grcount, 
              double *Fmin, int *fail,
              int n = 1, int m = 5, double factr = 1e7, double pgtol = 0.,
              int maxit = 100, int trace = 0, int nREPORT = 10)
  {
    char task[60];
    double f, *g,  *wa;
    int tr = -1, iter = 0, *iwa, isave[21];
    isave[12] = 0; // -Wall
    tr = -1; 
    *fail = 0;
    g = vect(n);
    /* this needs to be zeroed for snd in mainlb to be zeroed */ 
    // std::vector<double> wa_vec (2*m*n+4*n+11*m*m+8*m);
    // wa = &wa_vec[0];
    // std::vector<int> iwa_vec (3*n);
    // iwa = &iwa_vec[0];
    wa = (double *) S_alloc(2*m*n+4*n+11*m*m+8*m, sizeof(double));
    iwa = (int *) R_alloc(3*n, sizeof(int));
    strcpy(task, "START");
    while(1) {
      setulb(n, m, x, l, u, nbd, &f, g, factr, &pgtol, wa, iwa, task,
             tr, isave);
      /*	Rprintf("in lbfgsb - %s\n", task);*/
      if (strncmp(task, "FG", 2) == 0) {
        f = fminfn(x);
        if (!R_FINITE(f))
          break; // error(("L-BFGS-B needs finite values of 'fn'"));
        fmingr(x, g);
      } else if (strncmp(task, "NEW_X", 5) == 0) {
        iter++;
        if(trace == 1 && (iter % nREPORT == 0)) {
          Rprintf("iter %4d value %f\n", iter, f);
        }
        if (iter > maxit) {
          *fail = 1;
          break;
        }
      } else if (strncmp(task, "WARN", 4) == 0) {
        *fail = 51;
        break;
      } else if (strncmp(task, "CONV", 4) == 0) {
        break;
      } else if (strncmp(task, "ERROR", 5) == 0) {
        *fail = 52;
        break;
      } else { /* some other condition that is not supposed to happen */
      *fail = 52;
        break;
      }
    }
    *Fmin = f;
    *fncount = *grcount = isave[12];
    strcpy(msg, task);
  }
  
  
  // https://github.com/SurajGupta/r-source/blob/master/src/library/stats/src/optimize.c
  /* Formerly in src/appl/fmim.c 
  fmin.f -- translated by f2c (version 19990503).
  R's  optimize() :   function	fmin(ax,bx,f,tol)
  =    ==========		~~~~~~~~~~~~~~~~~
  an approximation  x  to the point where  f  attains a minimum  on
  the interval  (ax,bx)  is determined.
  INPUT..
  ax    left endpoint of initial interval
  bx    right endpoint of initial interval
  f     function which evaluates  f(x, info)  for any  x
  in the interval  (ax,bx)
  tol   desired length of the interval of uncertainty of the final
  result ( >= 0.)
  OUTPUT..
  fmin  abcissa approximating the point where  f  attains a minimum */
  
  double Brent_fmin(double ax, double bx, double tol = pow(DBL_EPSILON, 0.25))
  {
    /*  c is the squared inverse of the golden ratio */
    const double c = (3. - sqrt(5.)) * .5;
    
    /* Local variables */
    double a, b, d, e, p, q, r, u, v, w, x;
    double t2, fu, fv, fw, fx, xm, eps, tol1, tol3;
    
    /*  eps is approximately the square root of the relative machine precision. */
    eps = DBL_EPSILON;
    tol1 = eps + 1.;/* the smallest 1.000... > 1 */
    eps = sqrt(eps);
    
    a = ax;
    b = bx;
    v = a + c * (b - a);
    w = v;
    x = v;
    
    d = 0.;/* -Wall */
    e = 0.;
    fx = fminfn(&x);
    fv = fx;
    fw = fx;
    tol3 = tol / 3.;
    
    /*  main loop starts here ----------------------------------- */
    
    for(;;) {
      xm = (a + b) * .5;
      tol1 = eps * fabs(x) + tol3;
      t2 = tol1 * 2.;
      
      /* check stopping criterion */
      
      if (fabs(x - xm) <= t2 - (b - a) * .5) break;
      p = 0.;
      q = 0.;
      r = 0.;
      if (fabs(e) > tol1) { /* fit parabola */
      
      r = (x - w) * (fx - fv);
        q = (x - v) * (fx - fw);
        p = (x - v) * q - (x - w) * r;
        q = (q - r) * 2.;
        if (q > 0.) p = -p; else q = -q;
        r = e;
        e = d;
      }
      
      if (fabs(p) >= fabs(q * .5 * r) ||
          p <= q * (a - x) || p >= q * (b - x)) { /* a golden-section step */
      
      if (x < xm) e = b - x; else e = a - x;
      d = c * e;
      }
      else { /* a parabolic-interpolation step */
      
      d = p / q;
        u = x + d;
        
        /* f must not be evaluated too close to ax or bx */
        
        if (u - a < t2 || b - u < t2) {
          d = tol1;
          if (x >= xm) d = -d;
        }
      }
      
      /* f must not be evaluated too close to x */
      
      if (fabs(d) >= tol1)
        u = x + d;
      else if (d > 0.)
        u = x + tol1;
      else
        u = x - tol1;
      
      fu = fminfn(&u);
      
      /*  update  a, b, v, w, and x */
      
      if (fu <= fx) {
        if (u < x) b = x; else a = x;
        v = w;    w = x;   x = u;
        fv = fw; fw = fx; fx = fu;
      } else {
        if (u < x) a = u; else b = u;
        if (fu <= fw || w == x) {
          v = w; fv = fw;
          w = u; fw = fu;
        } else if (fu <= fv || v == x || v == w) {
          v = u; fv = fu;
        }
      }
    }
    /* end of main loop */
    
    return x;
  }
  
  
  
  
  Rcpp::List emp_mode_c(double & fix_lower, double & fix_upper){
    Rcpp::List res;
    double min_data = *data.begin();
    // double max_data = data.back();
    
    if(data.size() == 1){
      res["optimal alpha"] = R_NilValue;
      res["optimal p-value"] = R_NilValue;
      res["mode"] = min_data;
      Rcpp::Rcout << "The sample size is 1 and no smoothing is applied." << std::endl;
      return res;
    }

    // double delta = 1 / (double) n; 
    lower = fix_lower;
    upper = fix_upper;
      
    for(int i = 0; i < n; i++){
      data[i] = (data[i] - lower) / (upper - lower);
    }
    
    
    
    std::vector<double> newdata;
    if(n > subsample_size){
      std::vector<double> newdata_(subsample_size);
      for(int i = 0 ; i < subsample_size ; i++){
        newdata_[i] = data[n * (i + 1) / (subsample_size + 1)];
      }
      newdata = newdata_;
    }else{
      newdata = data;
    }
    
    
    double max_pvalue = -1;
    double max_alpha = -1; 
    if(do_choose_m){
    for (double alpha = 0.3; alpha <= 0.95; alpha += 0.05) {
      // Rcpp::Rcout << alpha;
      set_m((int) ceil(pow(n, alpha)));
      std::vector<double> v = cdf_bern_c(newdata);
      double* x = &v[0];
      double pvalue_curr = 1. - ADtest(newdata.size(), x);
      
      // int len = newdata.size();
      // // ADtest
      // double t, z = 0;
      // for(int i = 0; i < len; i++)   {
      //   t = mix_cdf_c(newdata[i]) * (1. - mix_cdf_c(newdata[len - 1 - i]));
      //   z = z - (i + i + 1) * log(t);
      // }
      // double pvalue_curr = AD(len, -len + z / len);
      
      
      // double adstat = ADstat(newdata.size(), x);
      if(pvalue_curr > max_pvalue){
        max_pvalue = pvalue_curr;
        max_alpha = alpha;
      }
      if(max_pvalue > 0.95){
        break;
      }
    }
    set_m((int) ceil(pow(n, max_alpha)));
    }else{
      std::vector<double> v = cdf_bern_c(newdata);
      double* x = &v[0];
      double pvalue_curr = 1. - ADtest(newdata.size(), x);
      max_pvalue = pvalue_curr;
    }
    
    
    // gradient descent
    double max_x = -1;
    double max_density = -1;
    
    mix_pdf_c(0.5);
    double arg_emp;
    double max_diff = -1;
    for(int i = 0; i < m; ++i){
      if(Fn_map_pdf_bern[m][i] > max_diff){
        arg_emp = i / (double) m;
        max_diff = Fn_map_pdf_bern[m][i];
      }
    }
    
    double pdfpar_vec[] = {0, 0.5, 1, arg_emp};
    for(int i = 0; (unsigned) i < sizeof(pdfpar_vec) / sizeof(pdfpar_vec[0]); ++i){
    try{
      // https://github.com/Microsoft/microsoft-r-open/blob/master/source/src/library/stats/src/optim.c
      double dpar = pdfpar_vec[i]; // 0.5, 1
      double lower = 0;
      double upper = 1;
      int nbd = 2;
      double val = 0.0;
      int ifail = 0;
      char msg[60];
      int fncount = 0, grcount = 0;
      lbfgsb(&dpar, &lower, &upper, &nbd, msg, &fncount, &grcount, &val, &ifail); 
      if(-val > max_density){
        max_x = dpar;
        max_density = -val;
      }
    }catch(std::exception e){
      Rcpp::Rcout << e.what() << std::endl;
    }
    }
    
    double pdfilon_vec[] = {0, 0.1, 0.5};
    for(int i = 0; (unsigned) i < sizeof(pdfilon_vec) / sizeof(pdfilon_vec[0]); ++i){
      try{
        double dpar = Brent_fmin(0 - pdfilon_vec[i], 1 + pdfilon_vec[i]);
        double val = mix_pdf_c(dpar);
        if(val > max_density){
          max_x = dpar;
          max_density = val;
        }
      }catch(std::exception e){
        Rcpp::Rcout << e.what() << std::endl;
      }
    }
          
    std::vector<double> & Fn_vec_cdf_bern = Fn_map_cdf_bern[m];
    double prob_left = Fn_vec_cdf_bern[2];
    double prob_right = 1 - Fn_vec_cdf_bern[Fn_vec_cdf_bern.size() - 3];
    int cum_index = 0;
    while(cum_index < n && data[cum_index] <= max_x - 1 / (double) m){
      cum_index ++ ;
    }
    double Fn_value_0 = cum_index / (double) n;
    while(cum_index < n && data[cum_index] <= max_x + 1 / (double) m){
      cum_index ++ ;
    }
    double Fn_value_1 = cum_index / (double) n;
    double prob_mid = Fn_value_1 - Fn_value_0;
    
    if(max_pvalue > 0.05 || (prob_mid > prob_left && prob_mid > prob_right)){
      Rcpp::Rcout << "The mode estimate is choosen as the maximum point of Bernstein polynomials density estimate." << std::endl;
      res["mode"] = max_x * (upper - lower) + lower;
    }else if(prob_left > prob_right){
      Rcpp::Rcout << "The mode estimate is choosen as a boundary point and Anderson-Darling criterion is not satisfied with p-value less than 0.05." << std::endl;
      res["mode"] = 0 * (upper - lower) + lower;
    }else if(prob_left < prob_right){
      Rcpp::Rcout << "The mode estimate is choosen as a boundary point and Anderson-Darling criterion is not satisfied with p-value less than 0.05." << std::endl;
      res["mode"] = 1 * (upper - lower) + lower;
    }else{
      Rcpp::Rcout << "The mode estimate is choosen as a boundary point and Anderson-Darling criterion is not satisfied with p-value less than 0.05." << std::endl;
      res["mode"] = floor(rand() / double(RAND_MAX) * 2) * (upper - lower) + lower;
    }
    
    if(do_choose_m){
      res["optimal alpha"] = max_alpha;
      res["optimal p-value"] = max_pvalue;
    }else{
      res["given degree"] = m;
      res["p-value"] = max_pvalue;
    }
    
    if(do_density){
      std::vector<double> x_cord(n_points);
      for(int i = 0; i < n_points; i++){
        x_cord[i] = i / (double) n_points;
      }
      res["y"] = pdf_bern_c(x_cord);
      for(int i = 0; i < n_points; i++){
        x_cord[i] = i / (double) n_points * (upper - lower) + lower;
      }
      res["x"] = x_cord;
    }
    
  return res;
  }
  
  
// private:
  //////////////////////
  std::vector<double> data;
  int n;
  int m;
  int subsample_size;
  double lower, upper;
  bool do_choose_m;
  bool do_density;
  int n_points;

  std::map<int, std::vector<double> > coef_map_cdf_bern;
  std::map<int, std::vector<double> > Fn_map_cdf_bern;
  std::map<int, double> u_map_cdf_bern;
  
  std::map<int, std::vector<double> > coef_map_pdf_bern;
  std::map<int, std::vector<double> > Fn_map_pdf_bern;
  std::map<int, double> u_map_pdf_bern;
  
  std::map<int, std::vector<double> > coef_map_dpdf_bern;
  std::map<int, std::vector<bool> > sign_map_dpdf_bern;
  std::map<int, std::vector<double> > Fn_map_dpdf_bern;
  std::map<int, double> u_map_dpdf_bern;
};








