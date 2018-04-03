#include "bayesm.h"
 
//---------------------------------EDITED FOR MULTIVARIATE ORDERED PROBIT GIBBS SAMPLER ------------------------------------

//EXTRA FUNCTIONS SPECIFIC TO THE MAIN FUNCTION--------------------------------------------

void print_double_in_C(double beta) {  
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
vec drawwi_mvop_n(vec const& w, vec const& mu, mat const& sigmai, int p, ivec y, vec y_index){
  //function to draw w_i as in an ordered multivariate probit fashion

  int ny = y.size();
  
  vec outwi = w;
  
  for(int i = 0; i < ny; i++){
	  //print_in_C(outwi[y_index[i]]);
	  //print_in_C(condmom(outwi, mu, sigmai, p, y_index[i]+1)[0]);
	  //print_in_C(condmom(outwi, mu, sigmai, p, y_index[i]+1)[1]);
	  
	  if(i == 0 && i+1 <ny && y[i+1] != 100){
	  //if(i == 0 && y[i] == 1 && i+1 <ny && y[i+1] != 100){
		// if it's the first observed response, sample from a truncated normal from above by the utility of the next response
		// (and it's not the last one, and the following is a ranked response)
	  	vec Cmout = condmom(outwi, mu, sigmai, p, y_index[i]+1);
		outwi[y_index[i]] = trunNorm(Cmout[0], Cmout[1], outwi[y_index[i+1]], 0);
		 
          }else if(i == 0 && i+1 <ny && y[i+1] == 100){
	  //}else if(i == 0 && y[i] == 1 && i+1 <ny && y[i+1] == 100){
		// if it's the first observed response, sample from a positive truncated normal
		// (and it's not the last one, and the following is a not-ranked response)
	  	vec Cmout = condmom(outwi, mu, sigmai, p, y_index[i]+1);
		outwi[y_index[i]] = trunNorm(Cmout[0], Cmout[1], 0.0, 0);
		
	  }else if(y[i] != 100 && i+1 < ny && y[i+1] != 100){
		// if it's another observed response, and the following response it's a ranked response, sample from a double-sided
		// truncated normal
		vec Cmout = condmom(outwi, mu, sigmai, p, y_index[i]+1);
		outwi[y_index[i]] = rtrunSc(Cmout[0], Cmout[1], outwi[y_index[i+1]], outwi[y_index[i-1]]);
	
          }else if(y[i] != 100 && i+1 < ny && y[i+1] == 100){
		// if it's another observed response, and the following response it's not a ranked response, sample from a
		//double-sided truncated normal, and truncated below by 0
		vec Cmout = condmom(outwi, mu, sigmai, p, y_index[i]+1);
		outwi[y_index[i]] = rtrunSc(Cmout[0], Cmout[1], 0.0, outwi[y_index[i-1]]);
		  
	 }else if(y[i] != 100 && i == ny-1){
	  	// if it's another observed response, and it's the last response, sample from a
		// truncated normal from above
		vec Cmout = condmom(outwi, mu, sigmai, p, y_index[i]+1);
		outwi[y_index[i]] = trunNorm(Cmout[0], Cmout[1], outwi[y_index[i-1]], 1);
		  	
	 }else{
	  	// if it's not a ranked response sample from a truncated normal, truncated above by 0
          	vec Cmout = condmom(outwi, mu, sigmai, p, y_index[i]+1);
	  	outwi[y_index[i]] = trunNorm(Cmout[0], Cmout[1], 0.0, 1);
	  }
	  
  //print_in_C(outwi[y_index[i]]);
  //print_line_in_C();
  }
	return (outwi);
}


vec drawwi_mvop(vec const& w, vec const& mu, mat const& sigmai, int p, ivec y, vec y_index){
  //function to draw w_i as in an ordered multivariate probit fashion

  int ny = y.size();
  
  vec outwi = w;
	
  int last_response_index = 1000;
  
  for(int i = 0; i < ny; i++){
	  
	  if(i == 0 && y[i] == 1 && i+1 <ny && y[i+1] != 100){
		// if it's the first observed response, sample from a truncated normal from above by the utility of the next response
		// (and it's not the last one, and the following is a ranked response)
	  	vec Cmout = condmom(outwi, mu, sigmai, p, y_index[i]+1);
		outwi[y_index[i]] = trunNorm(Cmout[0], Cmout[1], outwi[y_index[i+1]], 0);
		  
	  }else if(i == 0 && y[i] == 1 && i+1 <ny && y[i+1] == 100){
		// if it's the first observed response, sample from a normal distribution
		// (and it's not the last one, but it's the last ranked response)
	  	vec Cmout = condmom(outwi, mu, sigmai, p, y_index[i]+1);
		outwi[y_index[i]] = rtrunSc(Cmout[0], Cmout[1], -1000000000000000.0, 1000000000000000.0);
		last_response_index = y_index[i];
		
	  }else if(y[i] != 100 && i+1 < ny && y[i+1] != 100){
		// if it's another observed response, and the following response it's a ranked response, sample from a double-sided
		// truncated normal
		vec Cmout = condmom(outwi, mu, sigmai, p, y_index[i]+1);
		outwi[y_index[i]] = rtrunSc(Cmout[0], Cmout[1], outwi[y_index[i+1]], outwi[y_index[i-1]]);
		  
	  }else if(y[i] != 100 && i+1 < ny && y[i+1] == 100){
		// if it's another observed response, and the following response it's not a ranked response, sample from a
		// truncated normal (this is the last observed response, so we save the index)
		vec Cmout = condmom(outwi, mu, sigmai, p, y_index[i]+1);
		outwi[y_index[i]] = trunNorm(Cmout[0], Cmout[1], outwi[y_index[i-1]], 1);
		last_response_index = y_index[i];
		  
	 }else if(y[i] != 100 && i == ny-1){
	  	// if it's another observed response, and it's the last response, sample from a
		// truncated normal from above
		vec Cmout = condmom(outwi, mu, sigmai, p, y_index[i]+1);
		outwi[y_index[i]] = trunNorm(Cmout[0], Cmout[1], outwi[y_index[i-1]], 1);
	  
	  }else{
	  	// if it's not a ranked response, sample from a normal truncated above by the utility of the last ranked response
          	vec Cmout = condmom(outwi, mu, sigmai, p, y_index[i]+1);
	  	outwi[y_index[i]] = trunNorm(Cmout[0], Cmout[1], outwi[last_response_index], 1);
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
    
    outw.subvec(ind,ind+p-1) = drawwi_mvop_n(w.subvec(ind,ind+p-1),mu.subvec(ind,ind+p-1),sigmai,p,y_ordered,y_subindex);
  }
  
  return (outw);
}

//MAIN FUNCTION---------------------------------------------------------------------------------------
//[[Rcpp::export]]
List rmvpGibbs_rcpp_loop(int R, int keep, int nprint, int p, 
                         ivec const& y, mat const& X, vec const& beta0, mat const& sigma0, 
                         mat const& V, double nu, vec const& betabar, mat const& A) {
                       

  int n = y.size()/p;
  int k = X.n_cols;
  
  //allocate space for draws
  mat sigmadraw = zeros<mat>(R/keep, p*p);
  mat betadraw = zeros<mat>(R/keep,k);
  mat wdraw = zeros<mat>(R/10,y.size());
	
  vec wnew = zeros<vec>(X.n_rows);
  int suma;

  // create initial vector of utilities w
  for(int i=0; i<n; i++){
  	suma = 0;
	// count how many responses are not ranked (i.e., y[i]=100) in suma
  	for(int k=0; k<p; k++){
    		if(y[i*p + k] != 100){
      			suma = suma + 1;}}
  	for(int j=0; j<p; j++){
		// for every observed response, sample in a ordered and uniform fashion between 0 and 2 the utility w
    		if(y[i*p + j] != 100){
      			wnew[i*p + j] = 2.00001 - 2*(y[i*p + j]-1)/(double)suma;
		// for every not answered option, sample from negative uniform distribution between 0 and -2
    		}else{
      			wnew[i*p + j] = runif(1, -2, 0)[0];}}
}

  //set initial values of w, beta, sigma (or root of inv)
  vec wold = wnew;
  vec betaold = beta0;
  mat C = chol(solve(trimatu(sigma0),eye(sigma0.n_cols,sigma0.n_cols))); //C is upper triangular root of sigma^-1 (G) = C'C
                                                                         //trimatu interprets the matrix as upper triangular and makes solve more efficient
  
  mat sigmai, zmat, epsilon, S, IW, ucholinv, VSinv; 
  vec betanew;
  List W;
  
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

      wnew = draww_mvop(wold,X*betaold,sigmai,y_copy);

      //draw beta given w(rep) and sigma(rep-1)
      //  note:  if Sigma^-1 (G) = C'C then Var(Ce)=CSigmaC' = I
      //  first, transform w_i = X_ibeta + e_i by premultiply by C
      
      zmat = join_rows(wnew,X); //similar to cbind(wnew,X)
      zmat.reshape(p,n*(k+1));
      zmat = C*zmat;
      zmat.reshape(X.n_rows,k+1);
      
      betanew = breg(zmat(span::all,0),zmat(span::all,span(1,k)),betabar,A);
	   
      //draw sigmai given w and beta
      epsilon = wnew-X*betanew;
      epsilon.reshape(p,n);  
      S = epsilon*trans(epsilon);
      
      //same as chol2inv(chol(V+S))
      ucholinv = solve(trimatu(chol(V+S)), eye(p,p));
      VSinv = ucholinv*trans(ucholinv);
      
      W = rwishart(nu+n,VSinv);
      C = as<mat>(W["C"]); //conversion from Rcpp to Armadillo requires explict declaration of variable type using as<>
      
      //print time to completion
      if (nprint>0) if ((rep+1)%nprint==0) infoMcmcTimer(rep, R);
      
      //save every keepth draw
        if((rep+1)%keep==0){
          mkeep = (rep+1)/keep;
          betadraw(mkeep-1,span::all) = trans(betanew);
          IW  = as<mat>(W["IW"]);
          sigmadraw(mkeep-1,span::all) = trans(vectorise(IW));
         }
	    
      //save w draws every 10th draw
        if((rep+1)%10==0){
	  mkeep = (rep+1)/10;
	  wdraw(mkeep-1,span::all) = trans(wnew);
	 }
		
      wold = wnew;
      betaold = betanew;
    }
  
  if(nprint>0) endMcmcTimer();
      
  return List::create(
    Named("betadraw") = betadraw, 
    Named("sigmadraw") = sigmadraw,
    Named("wdraw") = wdraw);
    //use to save only the last w draw:
    //Named("wdraw") = wnew);
}
