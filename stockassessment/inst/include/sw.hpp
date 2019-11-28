template <class Type>
Type nllSW(dataSet<Type> &dat, confSet &conf, paraSet<Type> &par, array<Type> &logSW, vector<Type> &logCSW, data_indicator<vector<Type>,Type> &keep, objective_function<Type> *of){
  Type nll=0; 
  int stateDimSW=logSW.dim[0];
  int timeSteps=logSW.dim[1];
  Type sdLogSW=exp(par.logSdLogSW(0));
  matrix<Type> swvar(stateDimSW,stateDimSW);
  matrix<Type> swcor(stateDimSW,stateDimSW);
  vector<Type> swsd(stateDimSW);  

  
  swcor.setZero();

  for(int i=0; i<stateDimSW; ++i){
    swcor(i,i)=1.0;
  }


  int i,j;
  for(i=0; i<stateDimSW; ++i){
    swsd(i)=sdLogSW;
  }
 
  for(i=0; i<stateDimSW; ++i){
    for(j=0; j<stateDimSW; ++j){
      swvar(i,j)=swsd(i)*swsd(j)*swcor(i,j);
    }
  }
  
  //density::MVNORM_t<Type> neg_log_densityF(fvar);
  MVMIX_t<Type> neg_log_densitySW(swvar,Type(0.0));
  for(int i=1;i<timeSteps;i++){
    nll+=neg_log_densitySW(logSW.col(i)-logSW.col(i-1)); // F-Process likelihood
  }
  
  if(conf.stockWeightProcess==2){
    for(int i= 0; i<logCSW.size(); ++i){
      nll += -dnorm(logCSW(i), Type(0), exp(par.logSdLogCSW(0)), true);
    }
  }

  return nll;
}

template <class Type>
Type nllSWobs(dataSet<Type> &dat, confSet &conf, paraSet<Type> &par, array<Type> &logSW, vector<Type> &logCSW, data_indicator<vector<Type>,Type> &keep, objective_function<Type> *of){
  Type nll=0; 
  Type pred;
  array<Type> SW=dat.stockMeanWeight;
  for(int a=0; a<SW.dim[1]; ++a){
    for(int y=0; y<SW.dim[0]; ++y){
      if(y<a){
        pred=logSW(a,y);
      }else{
        pred=logSW(a,y)+logCSW(y-a);
      }
      nll += -dnorm(log(SW(y,a)), pred, exp(par.logSdLogSWObs(0)), true);
    }
  }
  return nll;
}
