// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>





using namespace Rcpp;

//------------------------------------ Class Specific Functions -------------------------------------------//
class pch{

public:

  // data
  Rcpp::List curr_tables; 
  Rcpp::List hist_tables;
  arma::vec a0_vec;
  arma::vec n_intervals;
  int P; // number of columns of X
  int D; // total number of parameters
  
  

  // model definition
  //bool borrrow_control;
  bool dCurrent; // if false, then current likelihood contribution is zero
  bool shared_blh; // if false, then historical data has separate baseline hazard parameters
  
  // priors
  std::string prior_beta;
  std::string prior_lambda;
  arma::vec prior_beta_mu;
  arma::vec prior_beta_sigma;
  arma::vec prior_lambda_hp1;
  arma::vec prior_lambda_hp2;
  arma::vec prior_lambdah_hp1; // if baseline hazard parameters are not shared
  arma::vec prior_lambdah_hp2;
  
  
  arma::vec                lower_limits;
  arma::vec                upper_limits;
  arma::vec                slice_widths;
  int m;


  // public member functions;
  pch(Rcpp::List & curr_tables0, Rcpp::List & hist_tables0, arma::vec a0_vec0, 
      arma::vec & n_intervals0, bool & shared_blh0, int P0,
      std::string prior_beta0, arma::vec prior_beta_mu0, arma::vec prior_beta_sigma0,
      std::string prior_lambda0, arma::vec prior_lambda_hp10, arma::vec prior_lambda_hp20,
      arma::vec prior_lambdah_hp10, arma::vec prior_lambdah_hp20,
      arma::vec & lower_limits0, arma::vec & upper_limits0, arma::vec & slice_widths0, bool & dCurrent0);


  double logFC(arma::vec & parm0, const int & p);

};

pch::pch( Rcpp::List & curr_tables0, Rcpp::List & hist_tables0, arma::vec a0_vec0,
          arma::vec & n_intervals0, bool & shared_blh0, int P0,
          std::string prior_beta0, arma::vec prior_beta_mu0, arma::vec prior_beta_sigma0, 
          std::string prior_lambda0, arma::vec prior_lambda_hp10, arma::vec prior_lambda_hp20,
          arma::vec prior_lambdah_hp10, arma::vec prior_lambdah_hp20,
         arma::vec & lower_limits0, arma::vec & upper_limits0, arma::vec & slice_widths0, bool & dCurrent0)
{

  curr_tables = curr_tables0;
  hist_tables = hist_tables0;
  a0_vec = a0_vec0;
  
  n_intervals = n_intervals0;
  shared_blh = shared_blh0;
  P = P0;
  
  if(shared_blh0==TRUE){
    D = P0 + sum(n_intervals0);
  }else{
    D = P0 + 2*sum(n_intervals0);
  }

  
  dCurrent = dCurrent0;
  prior_beta = prior_beta0;
  prior_beta_mu = prior_beta_mu0;
  prior_beta_sigma = prior_beta_sigma0;
  prior_lambda = prior_lambda0;
  prior_lambda_hp1 = prior_lambda_hp10;
  prior_lambda_hp2 = prior_lambda_hp20;
  prior_lambdah_hp1 = prior_lambdah_hp10;
  prior_lambdah_hp2 = prior_lambdah_hp20;


  lower_limits = lower_limits0;
  upper_limits = upper_limits0;
  slice_widths = slice_widths0;

  m=10;
}

// Define the log likelihood 
double pch::logFC(arma::vec & parm0, const int & p)
{
  //Rcout << "parm0: " << parm0 << "\n";
  // extract regression parameters;
  int n_int = sum(n_intervals);
  arma::vec beta = parm0.subvec(0,P-1);
  //arma::vec beta = {1};
  arma::vec lambda = parm0.subvec(P,P+n_int-1);
  //arma::vec lambda = {1,2};
  arma::vec lambda_h;
  if(shared_blh==FALSE){
    lambda_h = parm0.subvec(P+n_int,D-1);
  }else{
    lambda_h = lambda;
  }
  //lambda_h = lambda;
  //parm0 = {1,1,2,1,2};
  //Rcout << "beta: " << beta << "\n";
  //Rcout << "p: " << p << "\n";
  //Rcout << "lambda: " << lambda << "\n";
  //Rcout << "lambda_h: " << lambda_h << "\n";
  
  // compute likelihood contribution;
  double ll = 0;
  
  if (dCurrent==TRUE){
    
    
    // current data likelihood;
    int Sn = curr_tables.size(); // number of strata
    for(int s = 0; s < Sn; s++){
      arma::mat tb = curr_tables[s];
      //Rcout << "tb: " << tb << "\n";



      //arma::mat rt = rt_data[s];

      arma::vec cumu = cumsum(n_intervals);
      int J = n_intervals[s];
      
      arma::vec lmd;
      if (s==0){
        lmd = lambda.subvec(0,J-1);
      }else{
        lmd = lambda.subvec(cumu[s-1],cumu[s-1]+J-1);
      }
      //Rcout << "lmd: " << lmd << "\n";
      for(int j = 0; j < J; j++){
        
        arma::mat tb_j = tb.rows(find(tb.col(P)==(j+1)));
        arma::mat X = tb_j.cols(0,P-1);
        arma::mat stats = tb_j.cols(arma::span(tb_j.n_cols-2, tb_j.n_cols-1));
        
        ll +=   sum(stats.col(1)) * log(lmd[j])
              - lmd[j] * sum(exp(X * beta) % stats.col(0))
              +  sum((X * beta) % stats.col(1)); 
              
       //Rcout << "ll: " << sum(stats.col(1)) * log(lmd[j])
       //                                - lmd[j] * sum(exp(X * beta) % stats.col(0))
       //                                +  sum((X * beta) % stats.col(1)) << "\n";
       //Rprintf("the value of v[%i] : %f \n", sum(stats.col(1)) * log(lmd[j])
      //          - lmd[j] * sum(exp(X * beta) % stats.col(0))
      //           +  sum((X * beta) % stats.col(1)));
      }
      
      
      
    }
    //Rcout << "ll: " << ll << "\n";
   
  }
  
  
  if(hist_tables.size() != 0){
    
    //Rcout << "lambda_h: " << lambda_h << "\n";
    for(int k = 0; k < hist_tables.size(); k++){
      Rcpp::List histdata = hist_tables[k];
      double a0 = a0_vec[k];
      //Rcout << "a0: " << a0 << "\n";
      int Sn = histdata.size(); // number of strata
      for(int s = 0; s < Sn; s++){
        arma::mat tb = histdata[s];
        
        //arma::mat rt = rt_data[s];
        
        arma::vec cumu = cumsum(n_intervals);
        int J = n_intervals[s];
        
        arma::vec lmd;
        if (s==0){
          lmd = lambda_h.subvec(0,J-1);
        }else{
          lmd = lambda_h.subvec(cumu[s-1],cumu[s-1]+J-1);
        }
        for(int j = 0; j < J; j++){
          
          arma::mat tb_j = tb.rows(find(tb.col(P)==(j+1)));
          arma::mat X = tb_j.cols(0,P-1);
          arma::mat stats = tb_j.cols(arma::span(tb_j.n_cols-2, tb_j.n_cols-1));
          
          ll +=   a0*(sum(stats.col(1)) * log(lmd[j])
            - lmd[j] * sum(exp(X * beta) % stats.col(0))
            +  sum((X * beta) % stats.col(1)));
        }
     
      }
    }
  }
  
  //Rcout << "ll: " << ll << "\n";
  //Rcout << "p="<< p << "\n";
  if(p < P){
    //Rcout << "p="<< p << "\n";
    if(prior_beta == "Normal"){
      ll += R::dnorm(parm0[p], prior_beta_mu[p], prior_beta_sigma[p], TRUE);
    }
    //Rcout << "prior_beta: " << R::dnorm(parm0[p], prior_beta_mu[p], prior_beta_sigma[p], TRUE) << "\n";
  }else{
    
    if(p < (P+n_int)){
      
      if(prior_lambda == "Gamma"){ll +=  R::dgamma(parm0[p], prior_lambda_hp1[p-P], 1/(prior_lambda_hp2[p-P]), TRUE);}
      //Rcout << "prior_lambda: " << p << " - " << R::dgamma(parm0[p], prior_lambda_hp1[p-P], 1/(prior_lambda_hp2[p-P]), TRUE) << "\n";
      if(prior_lambda == "Log-normal"){ll += R::dlnorm(parm0[p], prior_lambda_hp1[p-P], prior_lambda_hp2[p-P], TRUE);}
      if(prior_lambda == "Improper"){ll += -log(parm0[p]);}
    }
    if(shared_blh==FALSE and (p >= (P+n_int))){
        if(prior_lambda == "Gamma"){ll +=  R::dgamma(parm0[p], prior_lambdah_hp1[p-P-n_int], 1/(prior_lambdah_hp2[p-P-n_int]), TRUE);}
        //Rcout << "prior_lambda0: " << p << " - " << R::dgamma(parm0[p], prior_lambdah_hp1[p-P-n_int], 1/(prior_lambdah_hp2[p-P-n_int]), TRUE) << "\n";
        if(prior_lambda == "Log-normal"){ll += R::dlnorm(parm0[p], prior_lambdah_hp1[p-P-n_int], prior_lambdah_hp2[p-P-n_int], TRUE);}
        if(prior_lambda == "Improper"){ll += -log(parm0[p]);}
      
    }
    
  }

  //Rcout << "post: " << ll << "\n";
  return  ll;
}


// slice sampler
void slice( arma::vec & parms, pch & b)
{

  double b0, f0, f0_L, f0_R, f0_x1, h0, L, R, V, J, K,w,lower,upper;
  arma::vec parm0;

  
  for (int p = 0; p < (b.D); p++)
  {
    //Rcout << "p: " << p << "\n";
    //Rcout << "parm " << parm0 << "\n";
    // create vector of parameters to modify for slice sampling;
    parm0 = parms;

    // extract slice width and parameter bounds;
    w     = b.slice_widths[p];
    lower = b.lower_limits[p];
    upper = b.upper_limits[p];

    // skip over fixed parameter values;
    if (lower==upper){parms(p) = lower;}
    else
    {
      // current value of the parameter in question;
      b0 = parm0(p);

      // calculate current full conditional value;
      f0 = b.logFC(parm0,p);

      // calculate height of the horizontal slice;
      h0 = f0 - R::rexp(1.0);

      // Calculate initial horizontal interval;
      L = parm0(p) - R::runif(0.0,1.0)*w;
      R = L+w;

      // Truncate bounds to support of the parameter space;
      L = std::max(L,lower);
      R = std::min(R,upper);

      // Step out;
      V = R::runif(0.0,1.0);
      J = floor(b.m*V);
      K = (b.m-1)-J;

      // compute log of full conditional at current boundaries;
      parm0(p) = L; f0_L = b.logFC(parm0,p);
      parm0(p) = R; f0_R = b.logFC(parm0,p);

      while(J>0 and h0<f0_L and L>=lower)
      {
        L        = L-w; if (L<=lower) {L=lower;}
        J        = J-1;
        parm0(p) = L;
        f0_L     = b.logFC(parm0,p);
      }
      while(K>0 and h0<f0_R and R<=upper)
      {
        R        = R+w; if (R>=upper) {R=upper;}
        K        = K-1;
        parm0(p) = R;
        f0_R     = b.logFC(parm0,p);
      }

      // perform rejection sampling;
      int stop  = 0;
      while(stop == 0)
      {
        parm0(p)     = L + R::runif(0.0,1.0)*(R-L);
        f0_x1        = b.logFC(parm0,p);

        if      ( h0       <  f0_x1 ) { parms(p) = parm0(p); stop = 1;  }
        else if ( parm0(p) <  b0    ) { L = parm0(p);                     }
        else if ( parm0(p) >= b0    ) { R = parm0(p);                     }

        if (-0.0000000001 <= L-R and L-R <= 0.0000000001)
        {
          parms(p)= 0.5*(L+R);
          stop      = 1;
        }
      }
    }
  }
}


// [[Rcpp::export]]
Rcpp::List phm_fixed_a0(Rcpp::List & curr_tables0, Rcpp::List & hist_tables0, arma::vec & a0_vec0,
                       arma::vec & n_intervals0, bool & shared_blh0, int P0, 
                       std::string prior_beta0, arma::vec prior_beta_mu0, arma::vec prior_beta_sigma0,
                       std::string prior_lambda0, arma::vec prior_lambda_hp10, arma::vec prior_lambda_hp20, 
                       arma::vec prior_lambdah_hp10, arma::vec prior_lambdah_hp20,
                       arma::vec & lower_limits0, arma::vec & upper_limits0, arma::vec & slice_widths0, 
                       int nMC, int nBI, bool & dCurrent0){
  
  // declare object and set values;
    pch b(curr_tables0,hist_tables0,a0_vec0,n_intervals0,shared_blh0,P0,prior_beta0,prior_beta_mu0,prior_beta_sigma0,
          prior_lambda0,prior_lambda_hp10,prior_lambda_hp20,
          prior_lambdah_hp10, prior_lambdah_hp20,
          lower_limits0,upper_limits0,slice_widths0,dCurrent0);

    
    // Construct container for mcmc samples;
    arma::mat samples(nMC,b.D);

    // create parameter vector container and initial values;
    arma::vec parms(b.D);
    for (int p=0; p<b.D; p++)
    {

      parms[p]= R::runif(0,1);

    }

    for (int s=-nBI;s<nMC;s++)
    {
      slice(parms,b);

      if (s>=0){	samples.row(s) = parms.t();	}
    }
    
    arma::mat beta_samps = samples.cols(0,P0-1);
    arma::mat lambda_samps;
    arma::mat lambda_h_samps;
    if(shared_blh0 == TRUE){
      lambda_samps = samples.cols(P0,b.D-1);
    }else{
      lambda_samps = samples.cols(P0,P0+sum(n_intervals0)-1);
      lambda_h_samps = samples.cols(P0+sum(n_intervals0), b.D-1);
    }
    
    
    Rcpp::List lams(n_intervals0.size());
    arma::vec cumu = cumsum(n_intervals0);
    
    for(int s = 0; s < n_intervals0.size(); s++){
      int J = n_intervals0[s];
      if (s==0){
        lams[s] = lambda_samps.cols(0,J-1);
      }else{
        lams[s] = lambda_samps.cols(cumu[s-1],cumu[s-1]+J-1);
      }
    }
    
    Rcpp::List lams_h(n_intervals0.size());
    if(shared_blh0 == FALSE){
      arma::vec cumu = cumsum(n_intervals0);
      
      for(int s = 0; s < n_intervals0.size(); s++){                            
        int J = n_intervals0[s];
        if (s==0){
          lams_h[s] = lambda_h_samps.cols(0,J-1);
        }else{
          lams_h[s] = lambda_h_samps.cols(cumu[s-1],cumu[s-1]+J-1);
        }
      }
    }

    if(shared_blh0 == TRUE){
      return(Rcpp::List::create(
             Rcpp::Named("beta_samples") = beta_samps,
             Rcpp::Named("lambda_samples") = lams));
    }else{
      return(Rcpp::List::create(
          Rcpp::Named("beta_samples") = beta_samps,
          Rcpp::Named("lambda_samples") = lams,
          Rcpp::Named("lambda0_samples") = lams_h));
    }

  
}


