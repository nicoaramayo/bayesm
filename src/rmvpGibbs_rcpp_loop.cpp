#include "bayesm.h"
 
//EDITED FOR MULTIVARIATE ORDERED PROBIT GIBBS SAMPLER -------------------------------------

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

vec drawwi_mvp(vec const& w, vec const& mu, mat const& sigmai, int p, ivec y,
                  mat const& X, vec const& betahat){
  //function to draw w_i as in the ordered probit model

  int ny = y.size();
  
  vec beta = betahat;
  
  vec outwi = w;
  
  for(int i = 0; i<ny; i++){	
	  
	  if (i<1 && y[i]<100){
		// if it's i's first observed response, sample from a truncated normal from above the previous iteration draw of w_i
		// y[i] = 100 is a default value for not choosing alternative i
	  	vec Cmout = condmom(outwi, mu, sigmai, p, i+1);
	  	outwi[i] = trunNorm(Cmout[0], Cmout[1], outwi[i], 0);
		  
	  } else if (y[i]<100 && y[i+1]<100){
	  	// if it's one of i's observed responses (and it's not the last observed response),
		// sample from a truncated normal by the last w_i draw above and below from the previous iteration 
		// draw of the following choice
                  vec Cmout = condmom(outwi, mu, sigmai, p, i+1);
		  outwi[i] = rtrunSc(Cmout[0], Cmout[1], outwi[i-1], outwi[i+1]);
		  //outwi[i] = rtrunVec(Cmout[0], Cmout[1], outwi[i-1], outwi[i+1]);
		  
	  } else if (y[i]<100){
	  	// if it's one of i's observed responses (and it's the last observed response),
		// sample from a truncated normal by the last w_i draw above and below from 0
                  vec Cmout = condmom(outwi, mu, sigmai, p, i+1);
		  outwi[i] = rtrunSc(Cmout[0], Cmout[1], outwi[i-1], 0.0);
		  //outwi[i] = rtrunVec(Cmout[0], Cmout[1], outwi[i-1], 0.0);
		  
	  } else {
		// if it's a non-selected choice, sample from a negative truncated normal
	 	vec Cmout = condmom(outwi, mu, sigmai, p, i+1);
	  	outwi[i] = trunNorm(Cmout[0], Cmout[1], 0.0, 1);
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
  
  for(int i = 0; i<n; i++){
    ind = p*i;
    outw.subvec(ind,ind+p-1) = drawwi_mvp(w.subvec(ind,ind+p-1),mu.subvec(ind,ind+p-1),sigmai,p,y.subvec(ind,ind+p-1),
                X, betahat);
  }
  
  return (outw);
}

//MAIN FUNCTION---------------------------------------------------------------------------------------
[[Rcpp::export]]
List my_rmvpGibbs_rcpp_loop(int R, int keep, int nprint, int p, 
                         ivec const& y, mat const& X, vec const& beta0, mat const& sigma0, 
                         mat const& V, double nu, vec const& betabar, mat const& A,
			vec const& w_0) {
                       

  int n = y.size()/p;
  int k = X.n_cols;
  
  //allocate space for draws
  mat sigmadraw = zeros<mat>(R/keep, p*p);
  mat betadraw = zeros<mat>(R/keep,k);
  vec wnew = w_0; // initial utility vector
  
  //set initial values of w,beta, sigma (or root of inv)
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
      //   y is n x (p-1) vector of binary (0,1) outcomes 
      //   beta is k x 1 vector
      //   sigmai is (p-1) x (p-1) 
          
      wnew = draww_mvp(wold,X*betaold,sigmai,y,X,betabar);
  
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
      
  return List::create(
    Named("betadraw") = betadraw, 
    Named("sigmadraw") = sigmadraw,
    Named("wdraw") = wnew);
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
  
  //set initial values of w,beta, sigma (or root of inv)
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
      //   y is n x (p-1) vector of binary (0,1) outcomes 
      //   beta is k x 1 vector
      //   sigmai is (p-1) x (p-1) 
          
      wnew = draww_mvp(wold,X*betaold,sigmai,y,X,betabar);
  
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
      
  return List::create(
    Named("betadraw") = betadraw, 
    Named("sigmadraw") = sigmadraw,
    Named("wdraw") = wnew);
}
