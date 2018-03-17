#include "bayesm.h"
 
//---------------------------------EDITED FOR MULTIVARIATE ORDERED PROBIT GIBBS SAMPLER ------------------------------------

//EXTRA FUNCTIONS SPECIFIC TO THE MAIN FUNCTION--------------------------------------------

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

vec drawwi_mvp(vec const& w, vec const& mu, mat const& sigmai, int p, ivec y,
                  mat const& X, vec const& betahat, vec y_index){
  //function to draw w_i as in the ordered probit model

  int ny = y.size();
  
  vec beta = betahat;
  
  vec outwi = w;
  
  for(int i = 0; i < ny; i++){
	  
	  if(i == 0 && i+1 <ny && y[i] == 1 && y[i+1] != 100){
		// if it's the first observed response, sample from a truncated normal from above the previous iteration draw of w_i
	  	vec Cmout = condmom(outwi, mu, sigmai, p, i+1);
		outwi[y_index[i]] = trunNorm(Cmout[0], Cmout[1], outwi[y_index[i+1]], 0);
		  
	  }else if(i == 0 && i+1 <ny && y[i] == 1){
		// if it's the first observed response, sample from a truncated normal from above the previous iteration draw of w_i
	  	vec Cmout = condmom(outwi, mu, sigmai, p, i+1);
		outwi[y_index[i]] = trunNorm(Cmout[0], Cmout[1], 0.0, 0);
		
	  }else if(y[i] != 100 && i+1 < ny && y[i+1] != 100){
		// if it's another observed response, and the following response it's a ranked response, sample from a double-sided
		// truncated normal
		vec Cmout = condmom(outwi, mu, sigmai, p, i+1);
		outwi[y_index[i]] = rtrunSc(Cmout[0], Cmout[1], outwi[y_index[i-1]], outwi[y_index[i+1]]);
		  
	  }else if(y[i] != 100 && i+1 < ny && y[i+1] == 100){
		// if it's another observed response, and the following response it's not a ranked response, sample from a
		// double-sided truncated normal, and truncated below by 0
		vec Cmout = condmom(outwi, mu, sigmai, p, i+1);
		outwi[y_index[i]] = rtrunSc(Cmout[0], Cmout[1], outwi[y_index[i-1]], 0.0);
		  
	 }else if(y[i] != 100 && i == ny-1){
	  	// if it's another observed response, and it's the last response, sample from a
		// double-sided truncated normal, and truncated below by 0
		vec Cmout = condmom(outwi, mu, sigmai, p, i+1);
		outwi[y_index[i]] = rtrunSc(Cmout[0], Cmout[1], outwi[y_index[i-1]], 0.0);
	
	 }else{
	  // if it's not a ranked response sample from a truncated normal, truncated above by 0
          vec Cmout = condmom(outwi, mu, sigmai, p, i+1);
	  outwi[y_index[i]] = trunNorm(Cmout[0], Cmout[1], 0.0, 1);
	  }
	  
  }
	return (outwi);
}


vec draww_mvp(vec const& w, vec const& mu, mat const& sigmai, ivec const& y,
                 mat const& X, vec const& betahat){

  //function to draw all w vector for all n obs
  
  int p = sigmai.n_cols;
  int n = w.size()/p;
  int ind; 
  vec outw = zeros<vec>(w.size());
  ivec y_ordered;
  vec y_subindex;
  
  for(int i = 0; i<n; i++){
    vec y_index = zeros<vec>(y.size());
    for(int j=0; j < p; j++){
       	 y_index[j] = j;
    }
 
    ind = p*i;
	  
    y_ordered = y.subvec(ind,ind+p-1);
    quicksort(y_ordered, y_index, 0, p-1);
    y_subindex = y_index.subvec(ind,ind+p-1);
    
    outw.subvec(ind,ind+p-1) = drawwi_mvp(w.subvec(ind,ind+p-1),mu.subvec(ind,ind+p-1),sigmai,p,y_ordered,
                X, betahat, y_subindex);
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
		// for every observed response, sample in a ordered and uniform fashion between 0 and 1 the utility w
    		if(y[i*p + j] != 100){
      			wnew[i*p + j] = 1 - (y[i*p + j]-1)/(double)suma;
		// for every not answered option, sample from negative uniform distribution between 0 and -1
    		}else{
      			wnew[i*p + j] = runif(1, -1, 0)[0];}}
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
	    
      ivec y_copy = ivec(y);
          
      wnew = draww_mvp(wold,X*betaold,sigmai,y_copy,X,betabar);

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
        
      wold = wnew;
      betaold = betanew;
    }
  
  if(nprint>0) endMcmcTimer();
	
  //vec y_index;
  //ivec y_copy = ivec(y);  
   //   for(int i = 0; i<n; i++){
  // 		 vec y_index = zeros<vec>(y.size());
  //               for(int j=0; j < p; j++){
  //     	 		y_index[j] = j;
  //               }
  //    }
  //quicksort(y_copy, y_index, 0, p-1);
      
  return List::create(
    Named("betadraw") = betadraw, 
    Named("sigmadraw") = sigmadraw,
    Named("wdraw") = wnew);
   // Named("ydraw") = y,
   // Named("yorddraw") = y_copy,
   // Named("yinddraw") = y_index);
}
