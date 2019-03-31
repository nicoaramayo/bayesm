#include "bayesm.h"
 
//---------EDITED FOR MULTIVARIATE ORDERED PROBIT GIBBS SAMPLER & BAYESIAN SIMULTANEOUS DEMAND AND SUPPLY ESTIMATION -------

//EXTRA FUNCTIONS SPECIFIC TO THE MAIN FUNCTION--------------------------------------------

void print_double_in_C(double beta) {  
    Rcout <<  beta << ";";
}

void print_vec_in_C(vec beta) {  
    Rcout <<  beta << ";";
}

void print_int_in_C(int beta) {  
    Rcout <<  beta << ";";
}

void print_line_in_C() {  
    Rcout <<  endl;
}

double rtrunSc(double mu, double sigma, double a, double b){
  
// N. Aramayo  

//function to draw from univariate truncated normal distribution
//a the lower bound for truncation
//b the upper bound for truncation

  double FA;
  double FB;
  double out;
	
  FA = R::pnorm((a-mu)/sigma,0,1,1,0);
  FB = R::pnorm((b-mu)/sigma,0,1,1,0);
  out = mu+sigma*R::qnorm(R::runif(0,1)*(FB-FA)+FA,0,1,1,0);
  
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


vec price_density(int p, vec const& sigma_s, vec const& price_s, vec const& fo_demand_s, vec const& demand_s,
		  vec const& gamma, vec const& z_s, vec const& fo_cost_s){
  //density of price for bayesian simultaneous demand and supply estimation

  vec price_density = zeros<vec>(p);
  double pi = 3.1415926;

  for(int s=0; s<p; s++){
	if(price_s[s] > 0){
	price_density[s] = 1/(sqrt(2*pi*sigma_s[0]))*exp(-1/(2*sigma_s[0])*(log(price_s[s] + pow(fo_demand_s[s], -1)*demand_s[s])
					                  - gamma[0]*z_s[s]))*eps(fo_cost_s[s]);
	}
 } return (price_density);
}

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
  	}
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
  } return (demand_s);
}

vec first_order_demand(vec const& beta, mat const& X, mat const& sigmai){
  //expected demand for the multivariate ordered probit
  // version without price interactions  (with price interaction in column k-1 is beta(k-1)*X(i*p + s, k-1))

  int p = sigmai.n_cols;
  int k = beta.n_cols;
  vec fo_demand = zeros<vec>(p);
  double pi = 3.1415926;
  int n_students = int(X.n_rows/p);
	
  for(int s = 0; s < p; s++){
  	for(int i = 0; i < n_students; i++){
		fo_demand[s] = fo_demand[s] + (exp(2*sqrt(2/pi)*dot(beta,X.row(i*p + s))/sigmai(s,s)) * 
					       sigmai(s,s) * beta(k-1)) /
			pow(1 + exp(2*sqrt(2/pi)*dot(beta,X.row(i*p + s))/sigmai(s,s)), 2);
  	}
  } return (fo_demand);
}

double first_order_demand_s(int s, vec const& beta, mat const& X, mat const& sigmai){
  //expected demand for the multivariate ordered probit
  // version without price interactions  (with price interaction in column k-1 is beta(k-1)*X(i*p + s, k-1))

  int p = sigmai.n_cols;
  int k = beta.n_cols;
  double fo_demand_s = 0;
  double pi = 3.1415926;
  int n_students = int(X.n_rows/p);
	
  for(int i = 0; i < n_students; i++){
	fo_demand_s = fo_demand_s + (exp(2*sqrt(2/pi)*dot(beta,X.row(i*p + s))/sigmai(s,s)) * 
					       sigmai(s,s) * beta(k-1)) /
			pow(1 + exp(2*sqrt(2/pi)*dot(beta,X.row(i*p + s))/sigmai(s,s)), 2);
  } return (fo_demand_s);
}

vec second_order_demand(vec const& beta, mat const& X, mat const& sigmai){
  //expected demand for the multivariate ordered probit
  // version without price interactions

  int p = sigmai.n_cols;
  int k = beta.n_cols;
  vec so_demand = zeros<vec>(p);
  double pi = 3.1415926;
  int n_students = int(X.n_rows/p);
	
  for(int s = 0; s < p; s++){
  	for(int i = 0; i < n_students; i++){
		so_demand[s] = so_demand[s] - (pow(sigmai(s,s), 2) * 2 * sqrt(2/pi)*pow(beta(k-1), 2) * 
			       exp(2*sqrt(2/pi)*dot(beta,X.row(i*p + s))/sigmai(s,s))) /
			       pow(1 + exp(2*sqrt(2/pi)*dot(beta,X.row(i*p + s))/sigmai(s,s)), 3);
  	}
  }
  return (so_demand);
}

double second_order_demand_s(int s, vec const& beta, mat const& X, mat const& sigmai){
  //expected demand for the multivariate ordered probit
  // version without price interactions

  int p = sigmai.n_cols;
  int k = beta.n_cols;
  double so_demand_s = 0;
  double pi = 3.1415926;
  int n_students = int(X.n_rows/p);
	
  for(int i = 0; i < n_students; i++){
	so_demand_s = so_demand_s - (pow(sigmai(s,s), 2) * 2 * sqrt(2/pi)*pow(beta(k-1), 2) * 
			       exp(2*sqrt(2/pi)*dot(beta,X.row(i*p + s))/sigmai(s,s))) /
			       pow(1 + exp(2*sqrt(2/pi)*dot(beta,X.row(i*p + s))/sigmai(s,s)), 3);
  }
  return (so_demand_s);
}

vec first_order_costshifter(vec const& beta, mat const& X, mat const& sigmai){
  //expected demand for the multivariate ordered probit
  // version without price interactions

  int p = sigmai.n_cols;
  int k = beta.n_cols;
  vec fo_costshifters = zeros<vec>(p);
  vec demand = zeros<vec>(p);
  vec fo_demand = zeros<vec>(p);
  vec so_demand = zeros<vec>(p);
	
  demand = expected_demand(beta, X, sigmai);
  fo_demand = first_order_demand(beta, X, sigmai);
  so_demand = second_order_demand(beta, X, sigmai);

  for(int s = 0; s < p; s++){
	fo_costshifters[s] = (2 - pow(fo_demand[s], -2) * so_demand[s] * demand[s]) /
		             (X(s, k-1) + pow(fo_demand[s], -1) * demand[s]);
  	}
  return (fo_costshifters);
}

double first_order_costshifter_s(int s, vec const& beta, mat const& X, mat const& sigmai){
  //expected demand for the multivariate ordered probit
  // version without price interactions

  int k = beta.n_cols;
  double fo_costshifters_s = 0;
  double demand_s = 0;
  double fo_demand_s = 0;
  double so_demand_s = 0;
	
  demand_s = expected_demand_s(s, beta, X, sigmai);
  fo_demand_s = first_order_demand_s(s, beta, X, sigmai);
  so_demand_s = second_order_demand_s(s, beta, X, sigmai);

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
                       double price_s, vec const& gamma, vec const& z_s){
  //density of price for bayesian simultaneous demand and supply estimation

  double pprice_s = 0;
  double pi = 3.1415926;
  double demand_s = 0;
  double fo_demand_s = 0;
  double fo_costshifters_s = 0;
	
  demand_s = expected_demand_s(s, beta, X, sigmai);
  fo_demand_s = first_order_demand_s(s, beta, X, sigmai);
  fo_costshifters_s = first_order_costshifter_s(s, beta, X, sigmai);
  
	
  pprice_s = double(1)/(sqrt(2*pi*sigma_s[0]))*exp(-1/(2*sigma_s[0])*(log(price_s + demand_s/fo_demand_s)
					                  - gamma[0]*z_s[s]))*std::abs(fo_costshifters_s);
  
  return (pprice_s);
}

mat rejection_price_sampler(int p, vec const& sigma_s, vec const& price_s,
		  		vec const& gamma, vec const& z_s, vec const& beta, mat X, mat const& sigmai){
	
	const int n_samples = 30000;
  //vec sample_x; sample_x.randn(n_samples); sample_x = sample_x*5000 + 1000;  // price range
  double lambda = 0.00001;
  vec sample_x = zeros<vec>(n_samples);
  vec sample_u; sample_u.randu(n_samples);
  int M = 2000;
  double pprice_s = 0;
  vec pexp = zeros<vec>(n_samples);
  vec price_samples = zeros<vec>(n_samples);
  vec pprices = zeros<vec>(n_samples);
  int error_counter = 0;
  
  mat accept_mask = zeros<mat>(n_samples, 2*p); //contains on the left-side matrix the sampled prices and on the right-side the acceptance mask
  double condition = 0;
  int k = X.n_cols; 
  int n_students = int(X.n_rows/p);
  
  //vec beta_conv = zeros<vec>(k);
  //mat domain = zeros<mat>(n_samples, 2*p);
  //double precio = 0;
  //double density = 0;
  
  //beta_conv[0] = 0.55134;
  //beta_conv[1] = 0.302344;
  //beta_conv[2] = 0.007699741;
  //beta_conv[3] = 0.7725113;
  //beta_conv[4] = -1.199107;
  //beta_conv[5] = -0.1773656;
  //beta_conv[6] = -1.519115;
  //beta_conv[7] = 0.1854648;
  //beta_conv[8] = -0.02959265;
  //beta_conv[9] = -0.3439252;
  //beta_conv[10] = 0.1120288;
  //beta_conv[11] = -2.148379e-06;
  
  for(int s = 0; s < p; s++){
    
    //if(price_s[s] > 0){
      //precio = 0;
      //for(int a = 0; a < 30000; a++){
      //  for(int j = 0; j < n_students; j++){
      //    X(s*p + j,k-1) = precio;  
      //  }
      //  density = price_density_s(s, beta_conv, X, sigmai, sigma_s, precio, gamma, z_s);
      //  domain(a,s) = precio;
      //  domain(a,s+p) = density;
      //  precio = precio + 10;
      //}
      
    //}
    
	  if(price_s[s] > 0){
	    sample_x = exponential_sample(lambda, n_samples);
	    error_counter = 0;
	    
		  for(int i = 0; i < n_samples; i++){
		    
			  for(int j = 0; j < n_students; j++){
				  X(s*p + j,k-1) = sample_x[i];   //replace in X the sampled price for school s
			  }
			  
			  pexp[i] = exponential_density(sample_x[i], lambda);
			  pprice_s = price_density_s(s, beta, X, sigmai, sigma_s, sample_x[i], gamma, z_s);
			  if(M*pexp[i] < pprice_s){
			    error_counter = error_counter + 1;
			    }
			  condition = pprice_s/(M*pexp[i]);
			  
			  if(sample_u[i] <= condition){
				  accept_mask(i,s) = sample_x[i];
				  accept_mask(i,s + p) = 1;
			  } else{accept_mask(i,s) = sample_x[i];}
		  } Rcout << "School: " << s << " Wrong samples: " << double(error_counter)/n_samples*100 << "%" << endl;
	} //else{for(int i = 0; i < n_samples; i++){accept_mask(i,s) = sample_x[i];}}
  } return (accept_mask);
}

vec sample_prices(int p, vec const& sigma_s, vec const& price_s,
                  vec const& gamma, vec const& z_s, vec const& beta, mat X, mat const& sigmai){
  
  double lambda = 0.00001;                  //exponential dist rate
  double sample_x = 0;                     //price samples from rejection sampling
  const int n_samples = 1;
  vec sample_u; sample_u.randu(n_samples);  //uniform (0,1) samples
  int M = 2000;
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
      sample_x = one_exponential_sample(lambda);
      error_counter = 0;
      sample_u.randu(n_samples);
      
        for(int j = 0; j < n_students; j++){
          X(j*p + s,k-1) = sample_x;   //replace in X the sampled price for school s
        }
        
        pexp = exponential_density(sample_x, lambda);
        pprice_s = price_density_s(s, beta, X, sigmai, sigma_s, sample_x, gamma, z_s);
        condition[0] = pprice_s/(M*pexp);
        arma_rng::set_seed(sample_u[0]);

        while(sample_u[0] > condition[0]){
          sample_x = one_exponential_sample(lambda);
          pexp = exponential_density(sample_x, lambda);
          pprice_s = price_density_s(s, beta, X, sigmai, sigma_s, sample_x, gamma, z_s);
          condition[0] = pprice_s/(M*pexp);
          sample_u.randu(n_samples);
          arma_rng::set_seed(sample_u[0]);
          
          error_counter = error_counter + 1;
          
          if(error_counter%100 == 0){
            Rcout << "llevas " << error_counter << " samples rechazados para el colegio " << s << endl;
            Rcout << "price is= " << sample_x << " and the target density= " << pprice_s << " and instrumental density = " << M*pexp << endl;
          }
          
          if(pprice_s > M*pexp){
            Rcout << "warning: instrumental density is below the target density for school: " << s << endl;
            Rcout << "price is= " << sample_x << " and the target density= " << pprice_s << " and instrumental density = " << M*pexp << endl;
          }
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
  //int z = 1; //set number of cost shifters manually
	
  vec demand = zeros<vec>(p);
  vec fo_demand = zeros<vec>(p);
  vec so_demand = zeros<vec>(p);
  vec fo_cost = zeros<vec>(p);
  vec gamma = zeros<vec>(1); gamma[0] = gamma[0] - 1;
  vec sigma_s = zeros<vec>(1); sigma_s[0] = sigma_s[0] + 1;
  vec price_density = zeros<vec>(p);
  mat sampled_prices_mask;
  vec sampled_prices = zeros<vec>(p);
	
  mat A_mod;  A_mod.eye(k-1,k-1)*0.01;  //edited for BSSD
  
  //allocate space for draws
  mat sigmadraw = zeros<mat>(R/keep, p*p);
  mat betadraw = zeros<mat>(R/keep,k-1);  //betas without considering the cost shifters
  mat pricedraw = zeros<mat>(R/keep,p);  
	
	
  vec wnew = zeros<vec>(X.n_rows);
  int suma;
	
  vec price = zeros<vec>(p);
  vec cost_shifter = zeros<vec>(p);
  for(int i = 0; i < p; i++){
  	price[i] = X(i,k-2);
	cost_shifter[i] = X(i,k-1);
  }
	
  int price_column = X.n_cols-2;
  int cost_shifter_column = X.n_cols-1;

  mat X_copy = zeros<mat>(X.n_rows, X.n_cols-1);
  int xrows = X_copy.n_rows;
	
  for(int i = 0; i < xrows; i++){
	  for(int j = 0; j < cost_shifter_column; j++){     //do not go through the last column, as it contains the cost shifters
		  if(j == price_column){             // if I'm at the price column
			  //X_copy(i,j) = 0;         //initialize the price variable at 0
			  X_copy(i,j) = X(i,j);
		  } else{X_copy(i,j) = X(i,j);}
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
  vec betaold = beta0.subvec(0,X_copy.n_cols-1);   //update size of beta vector according to modified X matrix
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
	    
      demand = expected_demand(betanew, X_copy, sigmai);
      fo_demand = first_order_demand(betanew, X_copy, sigmai);
      so_demand = second_order_demand(betanew, X_copy, sigmai);
      fo_cost = first_order_costshifter(betanew, X_copy, sigmai);
      //price_density = price_density(p, sigma_s, price, fo_demand, demand, gamma, cost_shifter, fo_cost);
      sampled_prices = sample_prices(p, sigma_s, price, gamma, cost_shifter, betanew, X_copy, sigmai);
      
      for(int s = 0; s < p; s++){
        if(sampled_prices[s] > 0){
          for(int i = 0; i < n; i++){
            X_copy(i*p + s,price_column) = sampled_prices[s];
          }
        }
      }
      
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
    
  sampled_prices_mask = rejection_price_sampler(p, sigma_s, price, gamma, cost_shifter, betanew, X_copy, sigmai);
  
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
