// Survey-based statistical catch-at-length model (SURBAL) without process error (R follows AR1)
// Known catchability
// create a discrete dnorm function for pla.
// constant Linf and vbk
// add error in estimation of initial number
// dev_logF follows two-way AR1 process
// add plus group
// use a more effcient discrete_pnorm function (taylor series expansion)
// model tailored to yellowtail flounder

#include <TMB.hpp>
#include <iostream>

#include "model_math.hpp"

template<class Type>
Type objective_function<Type>::operator() ()
{

  // input data;
  DATA_MATRIX(logN_at_len); // observed survey catch at length (L by Y)
  DATA_MATRIX(na_matrix); // indicating elements where N_at_len is NA.
  DATA_VECTOR(log_q); // catchability at length
  DATA_VECTOR(len_border); // breaking border between length bins
  DATA_VECTOR(age); // ages
  DATA_INTEGER(Y); // number of years
  DATA_INTEGER(A); // number of ages
  DATA_INTEGER(L); // number of length bins
  DATA_MATRIX(weight); // weight at length (L by Y)
  DATA_MATRIX(mat); // maturity at length (L by Y)
  DATA_SCALAR(M); // assumed constant natural mortality

  Type zero = 0.0;
  Type one = 1.0;

  //define parameters;

  // fixed effects
  PARAMETER(log_init_Z); // intital total mortality used to set initial age structure
  PARAMETER(log_std_log_N0); // variability of initial number at age

  PARAMETER(mean_log_R); // mean recruitment
  PARAMETER(log_std_log_R); // recruitment standard error
  PARAMETER(logit_log_R); // define the AR1 coefficient of recruitment variation

  PARAMETER(mean_log_F);
  PARAMETER(log_std_log_F);
  PARAMETER(logit_log_F_y);
  PARAMETER(logit_log_F_a);

  PARAMETER(log_vbk); // VB parameter
  PARAMETER(log_Linf); // VB parameter
  PARAMETER(t0); // VB parameter
  PARAMETER(log_cv_len); // cv of length-at-age

  PARAMETER(log_std_index); // measurement error of survey catch at legnth

  // random effects
  PARAMETER_VECTOR(dev_log_R); // recruitment deviation
  PARAMETER_ARRAY(dev_log_F); // deviation of fishing moratlity across years and ages
  PARAMETER_VECTOR(dev_log_N0); // deviation of initial number at age (except for the first age), so size = A-1

  // derived parameters
  Type init_Z = exp(log_init_Z);
  Type std_log_N0 = exp(log_std_log_N0);

  Type std_log_R = exp(log_std_log_R);
  Type phi_log_R = invlogit(logit_log_R);

  Type std_log_F = exp(log_std_log_F);
  Type phi_log_F_y = invlogit(logit_log_F_y);
  Type phi_log_F_a = invlogit(logit_log_F_a);

  Type vbk = exp(log_vbk);
  Type Linf = exp(log_Linf);
  Type cv_len = exp(log_cv_len);

  Type std_index = exp(log_std_index);

  vector<Type> N(Y);
  matrix<Type> NA(A,Y);
  matrix<Type> log_NA(A,Y);
  matrix<Type> NL(L,Y);
  matrix<Type> log_NL(L,Y);
  matrix<Type> Z(A,Y);
  matrix<Type> F(A,Y);
  matrix<Type> pla(L,A);

  using namespace density; // call functions in TMB density namespace

  pla = length_age_key(len_border, age, Linf, vbk, t0, cv_len);

  //compute Z,F and M
  for(int i = 0;i < A;++i){
    for(int j = 0;j < Y;++j){
        F(i,j) = exp(mean_log_F + dev_log_F(i,j));
        Z(i,j)=F(i,j)+M;
     }
  }

  // get annual recruitment
  vector<Type> log_Rec = dev_log_R + mean_log_R; // may need to plug in SR model here
  vector<Type> Rec = exp(log_Rec); // - half * std_log_R * std_log_R); // bias correction

  // The cohort model;
  //initializing first year
  log_NA(0,0) = log_Rec(0);
  NA(0,0) = Rec(0);
  //Type cz=zero;
  for(int i = 1;i < A;++i){
    log_NA(i,0) = log_NA(i-1,0) -init_Z  + dev_log_N0(i-1);
	NA(i,0) = exp(log_NA(i,0));
  }
  //compute numbers at age
  for(int j = 1;j < Y;++j){
    log_NA(0,j) = log_Rec(j);
	NA(0,j) = Rec(j);
    for(int i = 1;i < (A-1);++i){
		log_NA(i,j) = log_NA(i-1,j-1) - Z(i-1,j-1);
		NA(i,j) = exp(log_NA(i,j));
	}
	NA(A-1,j) = NA(A-2,j-1)*exp(-Z(A-2,j-1)) + NA(A-1,j-1)*exp(-Z(A-1,j-1)); // plus group
	log_NA(A-1,j) = log(NA(A-1,j));
  }

  // compute number at length
  for(int i=0; i<Y; ++i){
	NL.col(i) = pla * NA.col(i);
	N(i)=NA.col(i).sum();
  }
  log_NL= log(NL.array() + Type(1e-300));

  // compute biomass at age and length
  vector<Type> B(Y);
  vector<Type> SSB(Y);
  matrix<Type> BL(L,Y);
  matrix<Type> SBL(L,Y);
  vector<Type> CN(Y); // annual catch in numbers
  matrix<Type> CNA(A,Y);

  for(int i=0; i<L; ++i){
	  for(int j=0; j<Y; ++j){
		BL(i,j)=NL(i,j)*weight(i,j);
		SBL(i,j)=BL(i,j)*mat(i,j);
	  }
  }

  for(int i=0; i<A; ++i){
	  for(int j=0; j<Y; ++j){
		CNA(i,j)=NA(i,j)*(one - exp(-Z(i,j)))*(F(i,j)/Z(i,j));
	  }
  }

  for(int i=0; i<Y; ++i){
	  CN(i)=CNA.col(i).sum();
	  B(i)=BL.col(i).sum();
	  SSB(i)=SBL.col(i).sum();
  }

  //calculate index
  matrix<Type> Elog_index(L,Y);
  matrix<Type> resid_index(L,Y);
  for(int i = 0;i < L;++i){
	for(int j=0; j<Y; ++j){
  	Elog_index(i,j) = log_q(i) + log_NL(i,j);
		resid_index(i,j) = zero;
		if (na_matrix(i,j)>0) resid_index(i,j) = (logN_at_len(i,j) - Elog_index(i,j)) * na_matrix(i,j);
	}
  }

  //negative log-likelihoods
  // parallel_accumulator enables TMB::openmp() multi-thread gradient computation
  // When openmp(1), behaves identically to Type nll = 0
  parallel_accumulator<Type> nll(this);

  //index, the measurement error
  for(int i=0;i<L;++i){
	for(int j=0;j<Y;++j){
		if(na_matrix(i,j)>0){
			nll -= dnorm(resid_index(i,j),zero,std_index,true);
		}
	}
  }

  // logF
  nll += SCALE(SEPARABLE(AR1(phi_log_F_y),AR1(phi_log_F_a)),std_log_F)(dev_log_F);

  //recruitment estimates
  nll += SCALE(AR1(phi_log_R),std_log_R)(dev_log_R); // AR1 evaluation to generate nll for logR

  // initial population size
  nll -= dnorm(dev_log_N0,zero,std_log_N0,true).sum(); // iid evaluation of the initial population size

 // report results
  REPORT(pla);
  REPORT(N)
  REPORT(NL);
  REPORT(NA);
  REPORT(dev_log_N0);
  REPORT(std_log_N0);
  REPORT(B);
  REPORT(BL);
  REPORT(SSB);
  REPORT(SBL);
  REPORT(Z);
  REPORT(F);
  REPORT(CN);
  REPORT(CNA);
  REPORT(mean_log_F);
  REPORT(dev_log_F);
  REPORT(phi_log_F_y);
  REPORT(phi_log_F_a);
  REPORT(Rec);
  REPORT(mean_log_R);
  REPORT(dev_log_R);
  REPORT(std_log_R);
  REPORT(phi_log_R);
  REPORT(Elog_index);
  REPORT(resid_index);
  REPORT(std_index);

  REPORT(vbk);
  REPORT(Linf);
  REPORT(t0);
  REPORT(cv_len);

  ADREPORT(Rec);
  ADREPORT(N);
  ADREPORT(NA);
  ADREPORT(NL);
  ADREPORT(B);
  ADREPORT(SSB);
  ADREPORT(F);
  ADREPORT(CN);

  return nll;
}
