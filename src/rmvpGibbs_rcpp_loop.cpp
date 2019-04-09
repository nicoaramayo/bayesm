#include "bayesm.h"
 
//---------EDITED FOR MULTIVARIATE ORDERED PROBIT GIBBS SAMPLER & BAYESIAN SIMULTANEOUS DEMAND AND SUPPLY ESTIMATION -------

//EXTRA FUNCTIONS SPECIFIC TO THE MAIN FUNCTION--------------------------------------------

double rtrunSc(double mu, double sigma, double a, double b){
  
// N. Aramayo  

//function to draw from univariate truncated normal distribution
//a the lower bound for truncation
//b the upper bound for truncation

  double FA;
  double FB;
  double unif_sample;
  double out;
	
  FA = R::pnorm((a-mu)/sigma,0,1,1,0);
  FB = R::pnorm((b-mu)/sigma,0,1,1,0);
  //fix condition to avoid underflow when FA = FB = 1 (sampling between two numbers in a normal
  //distribution with prob = 0)
  if(FA == 1 && FB == 1){Rcout << "entered underflow condition: FA = FB = 1"  << endl;
    //FA = .9999999;
    //FB = .9999998;
    unif_sample = R::runif(0,1);
    if(a<b){out = unif_sample*b + (1-unif_sample)*a;}
    else{Rcout << "ERROR: LOWER BOUND GREATER THAN UPPER BOUND" << endl;
      out = mu+sigma*R::qnorm(R::runif(0,1)*(FB-FA)+FA,0,1,1,0);}
    Rcout << "lower (a) = " << a << " sample (out) = " << out << " upper (b) = " << b <<
             " mu = " << mu << " sigma = " << sigma <<  endl;
  }
  //if(FA == 0 && FB == 0){Rcout << "entered underflow condition: FA = FB = 0"  << endl;}
  else{out = mu+sigma*R::qnorm(R::runif(0,1)*(FB-FA)+FA,0,1,1,0);}
  
  return(out);
}

// ---------------- QUICK SORT IMPLEMENTATION -------------------
// This function sorts a vector in ascending order and saves the original position before the sorting in the first vector
void swap(ivec &v, vec &vy, int x, int y);

void quicksort(ivec &vecx,vec &vecy, int L, int R) {
    int i, j, mid, piv;
    i = L;
    j = R;
    mid = L + (R - L) / 2;
    piv = vecx[mid];

    while (i<R || j>L) {
        while (vecx[i] < piv)
            i++;
        while (vecx[j] > piv)
            j--;

        if (i <= j) {
            swap(vecx, vecy, i, j); //error=swap function doesnt take 3 arguments
            i++;
            j--;
        }
        else {
            if (i < R)
                quicksort(vecx,vecy, i, R);
            if (j > L)
                quicksort(vecx,vecy, L, j);
            return;
        }
    }
}

void swap(ivec &v,vec &vy, int x, int y) {
    int temp = v[x];
    v[x] = v[y];
    v[y] = temp;
    temp = vy[x];
    vy[x] = vy[y];
    vy[y] = temp;
}

//------version where not answered options are specified to have negative utility----------
vec drawwi_mvop(vec const& w, vec const& mu, mat const& sigmai, int p, ivec y, vec y_index){
  //function to draw w_i as in an ordered multivariate probit fashion

  int ny = y.size();
  
  vec outwi = w;
  
  for(int i = 0; i < ny; i++){
	  
	  if(i == 0 && y[i] != 100 && i+1 <ny && y[i+1] != 100){
	  //if(i == 0 && y[i] == 1 && i+1 <ny && y[i+1] != 100){
		// if it's the first observed response, sample from a truncated normal from below by the utility of the next response
		// (and it's not the last one, and the following is a ranked response)
	  	vec Cmout = condmom(outwi, mu, sigmai, p, y_index[i]+1);
		  outwi[y_index[i]] = trunNorm(Cmout[0], Cmout[1], outwi[y_index[i+1]], 0);
		
		 
    }else if(i == 0 &&  y[i] != 100){
	  //}else if(i == 0 && y[i] == 1 && i+1 <ny && y[i+1] == 100){
		// if it's the first observed response, sample from a positive truncated normal
		// (and there isn't a next ranked response)
	  	vec Cmout = condmom(outwi, mu, sigmai, p, y_index[i]+1);
		  outwi[y_index[i]] = trunNorm(Cmout[0], Cmout[1], 0.0, 0);
		
	  }else if(y[i] != 100 && i+1 < ny && y[i+1] != 100){
		// if it's an observed response that it's not the first one, and the following response it's a ranked response,
		// sample from a double-sided truncated normal
  		vec Cmout = condmom(outwi, mu, sigmai, p, y_index[i]+1);
  		outwi[y_index[i]] = rtrunSc(Cmout[0], Cmout[1], outwi[y_index[i+1]], outwi[y_index[i-1]]);
	
    }else if(y[i] != 100){
		// if it's an observed response that it's not the first one and there isn't a next ranked response,
		// sample from a double-sided truncated normal truncated below by 0
  		vec Cmout = condmom(outwi, mu, sigmai, p, y_index[i]+1);
  		outwi[y_index[i]] = rtrunSc(Cmout[0], Cmout[1], 0.0, outwi[y_index[i-1]]);
		  	
	 }else{
	  	// if it's not a ranked response sample from a truncated normal, truncated above by 0
      vec Cmout = condmom(outwi, mu, sigmai, p, y_index[i]+1);
	  	outwi[y_index[i]] = trunNorm(Cmout[0], Cmout[1], 0.0, 1);
	  }
  }
	return (outwi);
}


vec draww_mvop(vec const& w, vec const& mu, mat const& sigmai, ivec const& y){
  //function to draw all w vector for all n obs
  
  int p = sigmai.n_cols;
  int n = w.size()/p;
  int ind; 
  vec outw = zeros<vec>(w.size());
	
  ivec y_ordered;
  vec y_subindex;
  
  for(int i = 0; i < n; i++){
    y_subindex = zeros<vec>(p);
    for(int j=0; j < p; j++){
       	 y_subindex[j] = j;
    }
 
    ind = p*i;
	  
    y_ordered = y.subvec(ind,ind+p-1);
    quicksort(y_ordered, y_subindex, 0, p-1);
    
    outw.subvec(ind,ind+p-1) = drawwi_mvop(w.subvec(ind,ind+p-1),mu.subvec(ind,ind+p-1),sigmai,p,y_ordered,y_subindex);
  }
  
  return (outw);
}


//vec price_density(int p, vec const& sigma_s, vec const& price_s, vec const& fo_demand_s, vec const& demand_s,
//		              vec const& gamma, vec const& z_s, vec const& fo_cost_s){
  //density of price for bayesian simultaneous demand and supply estimation
  //NOT BEING USED AT THE MOMENT

//  vec price_density = zeros<vec>(p);
//  double pi = 3.1415926;

//  for(int s=0; s<p; s++){
//	if(price_s[s] > 0){
//	price_density[s] = 1/(sqrt(2*pi*sigma_s[0]))*exp(-1/(2*sigma_s[0])*(log(price_s[s] + pow(fo_demand_s[s], -1)*demand_s[s])
//					                  - gamma[0]*z_s[s]))*eps(fo_cost_s[s]);
//	}
// } return (price_density);
//}

vec expected_demand(vec const& beta, mat const& X, mat const& sigmai){
  //expected demand for the multivariate ordered probit
 
  int p = sigmai.n_cols;
  vec demand = zeros<vec>(p);
  double pi = 3.1415926;
  int n_students = int(X.n_rows/p);
	
  for(int s = 0; s < p; s++){
  	for(int i = 0; i < n_students; i++){
  		demand[s] = demand[s] + exp(2*sqrt(2/pi)*dot(beta,X.row(i*p + s))/sigmai(s,s)) /
  		            (1 + exp(2*sqrt(2/pi)*dot(beta,X.row(i*p + s))/sigmai(s,s)));
  	} //Rcout << "s = " << s << " beta = ";
  	  //for(int b = 0; b < beta.n_rows; b++){Rcout << beta[b] << " ";}
  	  //Rcout << " X.row(0*p + s) = " << X.row(0*p + s) << " dot(beta,X) = " << dot(beta,X.row(0*p + s))
      //      << " sigmai(s,s) = " << sigmai(s,s) << " demand[s] = " << demand[s] << endl;
  } return (demand);
}

double expected_demand_s(int s, vec const& beta, mat const& X, mat const& sigmai){
  //expected demand for the multivariate ordered probit
 
  int p = sigmai.n_cols;
  double demand_s = 0;
  double pi = 3.1415926;
  int n_students = int(X.n_rows/p);
	
  for(int i = 0; i < n_students; i++){
  	demand_s = demand_s + exp(2*sqrt(2/pi)*dot(beta,X.row(i*p + s))/sigmai(s,s)) /
  			       (1 + exp(2*sqrt(2/pi)*dot(beta,X.row(i*p + s))/sigmai(s,s)));
  } //Rcout << "s = " << s << " beta = ";
    //for(int b = 0; b < beta.n_rows; b++){Rcout << beta[b] << " ";}
    //Rcout << " X.row(0*p + s) = " << X.row(0*p + s) << " dot(beta,X) = " << dot(beta,X.row(0*p + s))
    //      << " sigmai(s,s) = " << sigmai(s,s) << " demand_s = " << demand_s << endl;
  return (demand_s);
}

vec first_order_demand(vec const& beta, mat const& X, mat const& sigmai,
                       vec const& pri, vec const& pie, vec const& lvl,
                       int n_price_interactions){
  //expected demand for the multivariate ordered probit
  
  int p = sigmai.n_cols;
  int k = beta.n_rows;
  vec fo_demand = zeros<vec>(p);
  double pi = 3.1415926;
  int n_students = int(X.n_rows/p);
  vec interactions = zeros<vec>(n_price_interactions);

  for(int s = 0; s < p; s++){
  	for(int i = 0; i < n_students; i++){
  	  interactions[0] = pri[i*p + s];
  	  interactions[1] = pie[i*p + s];
  	  interactions[2] = lvl[i*p + s];
  	  fo_demand[s] = fo_demand[s] + (exp(2*sqrt(2/pi)*dot(beta,X.row(i*p + s))/sigmai(s,s)) * 
  	                 2 * sqrt(2/pi) * (1/sigmai(s,s)) * (beta(k-n_price_interactions-1) +
					           dot(beta.subvec(k-n_price_interactions,k-1), interactions))) /
			               pow(1 + exp(2*sqrt(2/pi)*dot(beta,X.row(i*p + s))/sigmai(s,s)), 2);
  	}
  } return (fo_demand);
}

double first_order_demand_s(int s, vec const& beta, mat const& X, mat const& sigmai,
                            vec const& pri, vec const& pie, vec const& lvl,
                            int n_price_interactions){
  //expected demand for the multivariate ordered probit
  
  int p = sigmai.n_cols;
  int k = beta.n_rows;
  double fo_demand_s = 0;
  double pi = 3.1415926;
  int n_students = int(X.n_rows/p);
  vec interactions = zeros<vec>(n_price_interactions);
	
  for(int i = 0; i < n_students; i++){
    interactions[0] = pri[i*p + s];
    interactions[1] = pie[i*p + s];
    interactions[2] = lvl[i*p + s];
  	fo_demand_s = fo_demand_s + (exp(2*sqrt(2/pi)*dot(beta,X.row(i*p + s))/sigmai(s,s)) * 
  	              2 * sqrt(2/pi) * (1/sigmai(s,s)) * (beta(k-n_price_interactions-1) +
  					      dot(beta.subvec(k-n_price_interactions,k-1), interactions))) /
  			          pow(1 + exp(2*sqrt(2/pi)*dot(beta,X.row(i*p + s))/sigmai(s,s)), 2);
  } return (fo_demand_s);
}

vec second_order_demand(vec const& beta, mat const& X, mat const& sigmai,
                        vec const& pri, vec const& pie, vec const& lvl,
                        int n_price_interactions){
  //expected demand for the multivariate ordered probit

  int p = sigmai.n_cols;
  int k = beta.n_rows;
  vec so_demand = zeros<vec>(p);
  double pi = 3.1415926;
  int n_students = int(X.n_rows/p);
  vec interactions = zeros<vec>(n_price_interactions);
	
  for(int s = 0; s < p; s++){
  	for(int i = 0; i < n_students; i++){
  	  interactions[0] = pri[i*p + s];
  	  interactions[1] = pie[i*p + s];
  	  interactions[2] = lvl[i*p + s];
  	  so_demand[s] = so_demand[s] - (pow(2, 2) * pow(sqrt(2/pi), 2) * pow(sigmai(s,s), -2) *
                	   pow(beta(k-n_price_interactions-1) +
                	   dot(beta.subvec(k-n_price_interactions,k-1), interactions), 2) *
                     exp(2*sqrt(2/pi)*dot(beta,X.row(i*p + s))/sigmai(s,s)) *
                     (exp(2*sqrt(2/pi)*dot(beta,X.row(i*p + s))/sigmai(s,s)) - 1)) *
                     pow(1 + exp(2*sqrt(2/pi)*dot(beta,X.row(i*p + s))/sigmai(s,s)), -3);
  	}
  }
  return (so_demand);
}

double second_order_demand_s(int s, vec const& beta, mat const& X, mat const& sigmai,
                             vec const& pri, vec const& pie, vec const& lvl,
                             int n_price_interactions){
  //expected demand for the multivariate ordered probit

  int p = sigmai.n_cols;
  int k = beta.n_rows;
  double so_demand_s = 0;
  double pi = 3.1415926;
  int n_students = int(X.n_rows/p);
  vec interactions = zeros<vec>(n_price_interactions);
	
  for(int i = 0; i < n_students; i++){
    interactions[0] = pri[i*p + s];
    interactions[1] = pie[i*p + s];
    interactions[2] = lvl[i*p + s];
    so_demand_s = so_demand_s - (pow(2, 2) * pow(sqrt(2/pi), 2) * pow(sigmai(s,s), -2) *
                  pow(beta(k-n_price_interactions-1) +
                  dot(beta.subvec(k-n_price_interactions,k-1), interactions), 2) *
                  exp(2*sqrt(2/pi)*dot(beta,X.row(i*p + s))/sigmai(s,s)) *
                  (exp(2*sqrt(2/pi)*dot(beta,X.row(i*p + s))/sigmai(s,s)) - 1)) *
                  pow(1 + exp(2*sqrt(2/pi)*dot(beta,X.row(i*p + s))/sigmai(s,s)), -3);
  }
  return (so_demand_s);
}

vec first_order_costshifter(vec const& beta, mat const& X, mat const& sigmai,
                            vec const& pri, vec const& pie, vec const& lvl,
                            int n_price_interactions){
  //expected demand for the multivariate ordered probit

  int p = sigmai.n_cols;
  int k = beta.n_rows;
  vec fo_costshifters = zeros<vec>(p);
  vec demand = zeros<vec>(p);
  vec fo_demand = zeros<vec>(p);
  vec so_demand = zeros<vec>(p);
	
  demand = expected_demand(beta, X, sigmai);
  fo_demand = first_order_demand(beta, X, sigmai, pri, pie, lvl, n_price_interactions);
  so_demand = second_order_demand(beta, X, sigmai, pri, pie, lvl, n_price_interactions);

  for(int s = 0; s < p; s++){
	fo_costshifters[s] = (2 - pow(fo_demand[s], -2) * so_demand[s] * demand[s]) /
		                   (X(s, k-1) + pow(fo_demand[s], -1) * demand[s]);
  }
  return (fo_costshifters);
}

double first_order_costshifter_s(int s, vec const& beta, mat const& X, mat const& sigmai,
                                 vec const& pri, vec const& pie, vec const& lvl,
                                 int n_price_interactions){
  //expected demand for the multivariate ordered probit

  int k = beta.n_rows;
  double fo_costshifters_s = 0;
  double demand_s = 0;
  double fo_demand_s = 0;
  double so_demand_s = 0;
	
  demand_s = expected_demand_s(s, beta, X, sigmai);
  fo_demand_s = first_order_demand_s(s, beta, X, sigmai, pri, pie, lvl, n_price_interactions);
  so_demand_s = second_order_demand_s(s, beta, X, sigmai, pri, pie, lvl, n_price_interactions);

  fo_costshifters_s = (2 - pow(fo_demand_s, -2) * so_demand_s * demand_s) /
		                  (X(s, k-1) + pow(fo_demand_s, -1) * demand_s);
  
  return (fo_costshifters_s);
}

double uniform_density(double x, double up_lim, double low_lim){
	double punif = 0;
	if(x <= up_lim && x >= low_lim){
		punif = double(1)/(up_lim - low_lim);
	}
	else{punif = 0;}
  return (punif);
}
  
double exponential_density(double x, double lambda){
  return lambda * std::exp(-lambda*x);
}

double normal_density(double x, double m, double s){
  static const double inv_sqrt_2pi = 0.3989422804014327;
  double a = (x - m) / s;
  
  return inv_sqrt_2pi / s * std::exp(-0.5f * a * a);
}

double one_normal_sample(double m, double s){
  const int n_samples = 1;
  vec sample_n; sample_n.randn(n_samples);
  arma_rng::set_seed(sample_n[0]);
  return sample_n[0]*s + m;
}

vec exponential_sample(double lambda, int n_samples){
  vec sample_u; sample_u.randu(n_samples);
  return -log(sample_u)/lambda;
}

double one_exponential_sample(double lambda){
  const int n_samples = 1;
  vec sample_u; sample_u.randu(n_samples);
  arma_rng::set_seed(sample_u[0]);
  return -log(1-sample_u[0])/lambda;
}

double price_density_s(int s, vec const& beta, mat const& X, mat const& sigmai, vec const& sigma_s,
                       double price_s, vec const& gamma, vec const& z_s,
                       vec const& pri, vec const& pie, vec const& lvl,
                       int n_price_interactions){
  //density of price for bayesian simultaneous demand and supply estimation

  double pprice_s = 0;
  double pi = 3.1415926;
  double demand_s = 0;
  double fo_demand_s = 0;
  double fo_costshifters_s = 0;
	
  demand_s = expected_demand_s(s, beta, X, sigmai);
  fo_demand_s = first_order_demand_s(s, beta, X, sigmai, pri, pie, lvl, n_price_interactions);
  fo_costshifters_s = first_order_costshifter_s(s, beta, X, sigmai, pri, pie, lvl, n_price_interactions);
  
	
  pprice_s = double(1)/(sqrt(2*pi*sigma_s[0]))*exp(-1/(2*sigma_s[0])*
                       (log(price_s + demand_s/fo_demand_s)
					             - gamma[0]*z_s[s]))*std::abs(fo_costshifters_s);
  
  return (pprice_s);
}

mat rejection_price_sampler(int p, vec const& sigma_s, vec const& price_s,
		  		vec const& gamma, vec const& z_s, vec const& beta, 
		  		mat X, mat const& sigmai, int n_price_interactions,
		  		vec const& pri, vec const& pie, vec const& lvl){
	
	const int n_samples = 1000;
  //vec sample_x; sample_x.randn(n_samples); sample_x = sample_x*5000 + 1000;  // price range
  double lambda = 0.000001;
  vec sample_x = zeros<vec>(n_samples);
  vec sample_u; sample_u.randu(n_samples);
  double M = 2;
  double pprice_s = 0;
  vec pexp = zeros<vec>(n_samples);
  vec price_samples = zeros<vec>(n_samples);
  vec pprices = zeros<vec>(n_samples);
  int error_counter = 0;
  
  mat accept_mask = zeros<mat>(n_samples, 2*p); //contains on the left-side matrix the sampled prices and on the right-side the acceptance mask
  double condition = 0;
  int k = X.n_cols;
  int n_students = int(X.n_rows/p);
  
  mat domain = zeros<mat>(n_samples, 5*p);
  double precio = 0;
  double density = 0;
  
  for(int s = 0; s < p; s++){
    //THIS FOR LOOP IS FOR OBTAINING THE IMAGE SET OF THE PRICE DENSITY FUNCTION FOR EACH SCHOOL
    if(price_s[s] > 0){
      Rcout << "....calculando para el colegio " << s << endl;
      precio = 1000;
      for(int a = 0; a < n_samples; a++){
        for(int j = 0; j < n_students; j++){
          X(j*p + s,k-n_price_interactions-1) = precio;
          X(j*p + s,k-n_price_interactions) = precio * pri[j*p + s];
          X(j*p + s,k-n_price_interactions+1) = precio * pie[j*p + s];
          X(j*p + s,k-n_price_interactions+2) = precio * lvl[j*p + s];
        }
        density = price_density_s(s, beta, X, sigmai, sigma_s, precio,
                                  gamma, z_s, pri, pie, lvl, n_price_interactions);
        domain(a,s) = precio;
        domain(a,s+p) = density;
        domain(a,s+2*p) = expected_demand_s(s, beta, X, sigmai);
        domain(a,s+3*p) = first_order_demand_s(s, beta, X, sigmai, pri, pie, lvl, n_price_interactions);
        domain(a,s+4*p) = second_order_demand_s(s, beta, X, sigmai, pri, pie, lvl, n_price_interactions);
        precio = precio + 1000;
      }
    } 
  } return (domain);
	  //if(price_s[s] > 0){
	    //sample_x = exponential_sample(lambda, n_samples);
	    //error_counter = 0;
	    
		  //for(int i = 0; i < n_samples; i++){
		    
			  //for(int j = 0; j < n_students; j++){
				  //X(s*p + j,k-1) = sample_x[i];   //replace in X the sampled price for school s
			  //}
			  
			  //pexp[i] = exponential_density(sample_x[i], lambda);
			  //pprice_s = price_density_s(s, beta, X, sigmai, sigma_s, sample_x[i], gamma, z_s, n_price_interactions);
			  //if(M*pexp[i] < pprice_s){
			    //error_counter = error_counter + 1;
			    //}
			  //condition = pprice_s/(M*pexp[i]);
			  
			  //if(sample_u[i] <= condition){
				  //accept_mask(i,s) = sample_x[i];
				  //accept_mask(i,s + p) = 1;
			  //} else{accept_mask(i,s) = sample_x[i];}
		  //} Rcout << "School: " << s << " Wrong samples: " << double(error_counter)/n_samples*100 << "%" << endl;
	  //} //else{for(int i = 0; i < n_samples; i++){accept_mask(i,s) = sample_x[i];}}
  //} //return (accept_mask);
}

vec sample_prices(int p, vec const& sigma_s, vec const& price_s,
                  vec const& gamma, vec const& z_s, vec const& beta, mat X, mat const& sigmai, 
                  int n_price_interactions, vec const& pri, vec const& pie, vec const& lvl){
  
  double lambda = 0.00001;                  //exponential dist rate
  double sample_x = 0;                     //price samples from rejection sampling
  const int n_samples = 1;
  vec sample_u; sample_u.randu(n_samples);  //uniform (0,1) samples
  double M = 2;
  double pprice_s = 0;                     //price density
  double pexp = 0;                         //instrumental density
  vec price_samples = zeros<vec>(p);       // accepted sampled prices
  
  int error_counter = 0;
  
  vec condition = zeros<vec>(1);
  int k = X.n_cols; 
  int n_students = int(X.n_rows/p);
  
  vec warning_counter = zeros<vec>(p);
  
  for(int s = 0; s < p; s++){
    if(price_s[s] > 0){
      
      do{
        sample_x = one_exponential_sample(lambda);
        sample_u.randu(n_samples);
        
        for(int j = 0; j < n_students; j++){
          X(j*p + s,k-n_price_interactions-1) = sample_x;
          X(j*p + s,k-n_price_interactions) = sample_x * pri[j*p + s];
          X(j*p + s,k-n_price_interactions+1) = sample_x * pie[j*p + s];
          X(j*p + s,k-n_price_interactions+2) = sample_x * lvl[j*p + s];
        }
        
        pexp = exponential_density(sample_x, lambda);
        pprice_s = price_density_s(s, beta, X, sigmai, sigma_s, sample_x, gamma, z_s,
                                   pri, pie, lvl, n_price_interactions);
        condition[0] = pprice_s/(M*pexp);
        arma_rng::set_seed(sample_u[0]);
        error_counter = error_counter + 1;
        if(error_counter%100 == 0){
          Rcout << "llevas " << error_counter << " samples rechazados para el colegio " << s << endl;
          Rcout << "price is= " << sample_x << " and the target density= " << pprice_s << " and instrumental density = " << M*pexp << endl;
        }
      }while (sample_u[0] > condition[0]);

      if(pprice_s > M*pexp){
        Rcout << "warning: instrumental density is below the target density for school: " << s << endl;
        Rcout << "price is= " << sample_x << " and the target density= " << pprice_s << " and instrumental density = " << M*pexp << endl;
      }
      price_samples[s] = sample_x;
    }
  } return (price_samples);
}

//MAIN FUNCTION---------------------------------------------------------------------------------------
//[[Rcpp::export]]
List rmvpGibbs_rcpp_loop(int R, int keep, int nprint, int p, 
                         ivec const& y, mat const& X, vec const& beta0, mat const& sigma0, 
                         mat const& V, double nu, vec const& betabar, mat const& A) {

  int n = y.size()/p;   //number of students
  int k = X.n_cols;
  int n_cost_shifters = 1;
  int n_price_interactions = 3;
	
  vec demand = zeros<vec>(p);
  vec fo_demand = zeros<vec>(p);
  vec so_demand = zeros<vec>(p);
  vec fo_cost = zeros<vec>(p);
  vec gamma = zeros<vec>(n_cost_shifters); gamma[0] = gamma[0] - 0.5293939;    //EDITAR ESTO CUANDO HAYAN MÁS COST SHIFTERS
  vec sigma_s = zeros<vec>(1); sigma_s[0] = sigma_s[0] + 1;
  vec price_density = zeros<vec>(p);
  mat sampled_prices_mask;
  vec sampled_prices = zeros<vec>(p);
	
  mat A_mod;  A_mod.eye(k-n_cost_shifters,k-n_cost_shifters)*0.05;  //edited for BSSD
  
  //allocate space for draws
  mat sigmadraw = zeros<mat>(R/keep, p*p);
  mat betadraw = zeros<mat>(R/keep,k-n_cost_shifters);  //betas without considering the cost shifters
  mat pricedraw = zeros<mat>(R/keep,p);  
	
  vec wnew = zeros<vec>(X.n_rows);
  int suma;
	
  vec price = zeros<vec>(p);
  vec pri = zeros<vec>(X.n_rows);
  vec pie = zeros<vec>(X.n_rows);
  vec lvl = zeros<vec>(X.n_rows);
  vec cost_shifter = zeros<vec>(p); //CAMBIAR A MATRIZ LUEGO CUANDO HAYAN VARIOS COST SHIFTERS
  
  for(int i = 0; i < p; i++){
  	price[i] = X(i,k-n_price_interactions-n_cost_shifters-1);
	  cost_shifter[i] = X(i,k-1);  //CAMBIAR A MATRIZ LUEGO CUANDO HAYAN VARIOS COST SHIFTERS
  }
	
  int price_column = X.n_cols-1-n_price_interactions-n_cost_shifters; 
  // 1 column of price, the next 3 are price interactions
  int cost_shifter_column = X.n_cols-n_cost_shifters;

  mat X_copy = zeros<mat>(X.n_rows, X.n_cols-n_cost_shifters); //remove from X the cost shifter columns
  int xrows = X_copy.n_rows;
	
  for(int i = 0; i < xrows; i++){
	  for(int j = 0; j < cost_shifter_column; j++){     //do not go through the last column, as it contains the cost shifters
      X_copy(i,j) = X(i,j);
	    if(j == price_column+1){
	      pri[i] = X(i,j);
	      if(pri[i] == 1 && X(i,price_column) > 0){
	        X_copy(i,j) = X(i,price_column);
	      } else{X_copy(i,j) = 0;}
	    } if(j == price_column+2){
	      pie[i] = X(i,j);
	      if(pie[i] == 1 && X(i,price_column) > 0){
	        X_copy(i,j) = X(i,price_column);
	      } else{X_copy(i,j) = 0;}
	    } if(j == price_column+3){
	      lvl[i] = X(i,j);
	      if(lvl[i] == 1 && X(i,price_column) > 0){
	        X_copy(i,j) = X(i,price_column);
	      } else{X_copy(i,j) = 0;}
	    }
	  }
  }
 //X.submat( first_row, first_col, last_row, last_col ) for subsetting matrix
		  

  // create initial vector of utilities w
  for(int i=0; i<n; i++){
  	suma = 0;
	// count how many responses are not ranked (i.e., y[i]=100) in suma
  	for(int k=0; k<p; k++){
    		if(y[i*p + k] != 100){
      			suma = suma + 1;}}
  	for(int j=0; j<p; j++){
		// for every observed response, sample in a ordered and uniform fashion between 0 and 5 the utility w
    		if(y[i*p + j] != 100){
      			wnew[i*p + j] = 5.00001 - 5*(y[i*p + j]-1)/(double)suma;
		// for every not answered option, sample from negative uniform distribution between 0 and -10
    		}else{
      			wnew[i*p + j] = runif(1, -10, 0)[0];}}
}

  //set initial values of w, beta, sigma (or root of inv)
  vec wold = wnew;
  vec betaold = beta0.subvec(0,X_copy.n_cols-n_cost_shifters);   //update size of beta vector according to modified X matrix
  //vec betaold = beta0;
  mat C = chol(solve(trimatu(sigma0),eye(sigma0.n_cols,sigma0.n_cols))); //C is upper triangular root of sigma^-1 (G) = C'C
                                                                         //trimatu interprets the matrix as upper triangular and makes solve more efficient
  
  mat sigmai, zmat, epsilon, S, IW, ucholinv, VSinv; 
  vec betanew;
  List W;
  k = X_copy.n_cols;  //reassign k
  vec betabar_mod = betabar.subvec(0,k-1);   //update size of beta vector according to modified X matrix
  
  // start main iteration loop
  int mkeep = 0;
  
  if(nprint>0) startMcmcTimer();
  
    for(int rep = 0; rep<R; rep++) {
    
      //draw w given beta(rep-1),sigma(rep-1)
      sigmai = trans(C)*C;
	    
	    
      //draw latent vector
      
      //w is n x (p-1) vector
      //   X ix n(p-1) x k  matrix
      //   y is n x (p-1) vector of ranked 1 to p responses, not ranked responses are 100 by default
      //   beta is k x 1 vector
      //   sigmai is (p-1) x (p-1) 
	    
      // create a copy of the vector of responses as to not modify the original order of the vector
      ivec y_copy = ivec(y);

      //Rcout << "wold: " << wold << endl;
      wnew = draww_mvop(wold,X_copy*betaold,sigmai,y_copy);

      //draw beta given w(rep) and sigma(rep-1)
      //  note:  if Sigma^-1 (G) = C'C then Var(Ce)=CSigmaC' = I
      //  first, transform w_i = X_ibeta + e_i by premultiply by C
      
      zmat = join_rows(wnew,X_copy); //similar to cbind(wnew,X)
      zmat.reshape(p,n*(k+1));
      zmat = C*zmat;
      zmat.reshape(X_copy.n_rows,k+1);
      betanew = breg(zmat(span::all,0),zmat(span::all,span(1,k)),betabar_mod,A_mod);  //A prior modified to A_mod, betabar also updated
      //draw sigmai given w and beta
      epsilon = wnew-X_copy*betanew;
      epsilon.reshape(p,n);  
      S = epsilon*trans(epsilon);
      
      //same as chol2inv(chol(V+S))
      ucholinv = solve(trimatu(chol(V+S)), eye(p,p));
      VSinv = ucholinv*trans(ucholinv);
      
      W = rwishart(nu+n,VSinv);
      C = as<mat>(W["C"]); //conversion from Rcpp to Armadillo requires explict declaration of variable type using as<>
	    
      //price_density = price_density(p, sigma_s, price, fo_demand, demand, gamma, cost_shifter, fo_cost);
      sampled_prices = sample_prices(p, sigma_s, price, gamma, cost_shifter, betanew, X_copy, sigmai,
                                     n_price_interactions, pri, pie, lvl);
      
      for(int s = 0; s < p; s++){
        if(sampled_prices[s] > 0){
          for(int i = 0; i < n; i++){
            X_copy(i*p + s, price_column) = sampled_prices[s];
            X_copy(i*p + s, price_column + 1) = sampled_prices[s] * pri[i*p + s];
            X_copy(i*p + s, price_column + 2) = sampled_prices[s] * pie[i*p + s];
            X_copy(i*p + s, price_column + 3) = sampled_prices[s] * lvl[i*p + s];
          }
        }
      }
      
      demand = expected_demand(betanew, X_copy, sigmai);
      fo_demand = first_order_demand(betanew, X_copy, sigmai, pri, pie, lvl, n_price_interactions);
      so_demand = second_order_demand(betanew, X_copy, sigmai, pri, pie, lvl, n_price_interactions);
      fo_cost = first_order_costshifter(betanew, X_copy, sigmai, pri, pie, lvl, n_price_interactions);
      
      //print time to completion
      if (nprint>0) if ((rep+1)%nprint==0) infoMcmcTimer(rep, R);
      
      //save every keepth draw
        if((rep+1)%keep==0){
          mkeep = (rep+1)/keep;
          betadraw(mkeep-1,span::all) = trans(betanew);
          IW  = as<mat>(W["IW"]);
          sigmadraw(mkeep-1,span::all) = trans(vectorise(IW));
          pricedraw(mkeep-1,span::all) = trans(sampled_prices);
         }
		
      wold = wnew;
      betaold = betanew;
    }
    
  sampled_prices_mask = rejection_price_sampler(p, sigma_s, price, gamma, cost_shifter, betanew,
                                                X_copy, sigmai, n_price_interactions, pri, pie, lvl);
  
  if(nprint>0) endMcmcTimer();
      
  return List::create(
    Named("betadraw") = betadraw, 
    Named("sigmadraw") = sigmadraw,
    //Named("wdraw") = wdraw);
    //use to save only the last w draw:
    Named("wdraw") = wnew,
    Named("modified_X") = X_copy,
    Named("expected_demand") = demand,
    Named("fo_demand") = fo_demand,
    Named("so_demand") = so_demand,
    Named("fo_cost") = fo_cost,
    Named("pricedraw") = pricedraw,
    Named("rejection_sampling") = sampled_prices_mask);
}
