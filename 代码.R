
# total iteration
pdd=200

#  value 1 
rhob_0g=rep(0,pdd)
rhob_maxg=rep(0,pdd)
rhob_errorg=rep(0,pdd)
sigmafg=rep(0,pdd)
betafg=rep(0,pdd)

#  value 2 
rhob_0gg=rep(0,pdd)
rhob_maxgg=rep(0,pdd)
rhob_errorgg=rep(0,pdd)
sigmafgg=rep(0,pdd)
betafgg=rep(0,pdd)

for(r_8 in 1:pdd){ 
  
  #basic parameters
  
  #Matrix dimension 
  
  n_0=200
  
  
  #spatial weighted matrix
  
  # standard matrix
  I_1=diag(n_0)
  
  #previously defined matrix2
  
  
  W_3=matrix(0, nrow=n_0, ncol=n_0)
  d_3=matrix(0, nrow=n_0, ncol=n_0)
  sum3_nrowd=rep(0,n_0)
  #31
  for(i_0 in 1:n_0){
    for(j_0 in 1:n_0){
      if(j_0==i_0){d_3[i_0,j_0]=d_3[i_0,j_0]
      }
      else{
        d_3[i_0,j_0]=1/abs(i_0-j_0)
      }
    }
    sum3_nrowd[i_0]=sum(d_3[i_0,])
  }
  #32
  for(i_1 in 1:n_0){
    for(j_1 in 1:n_0){
      W_3[i_1,j_1]=d_3[i_1,j_1]/(sum3_nrowd[i_1])
    }
  }
  
  ym=eigen(W_3)
  
  
  
  ymm=-1/(abs(min(ym$val))) 
  
  
  #real values of parameters
  
  beta_0=1
  sigma_0=1
  rho_0=0.5
  
  
  
  #noise
  
  epsilon_sn=matrix(0, nrow=n_0, ncol=1)
  
  for(i_2 in 1:n_0){
    
    epsilon_sn[i_2,1]=rnorm(1,0,1)
    
  }
  
  #No contamination, generated data of model1 
  
  S_1n=matrix(0, nrow=n_0, ncol=n_0)
  S_1n=solve(I_1-rho_0*W_3)
  
  X_3n=matrix(beta_0, nrow=n_0, ncol=1)
  
  Y_1=matrix(0, nrow=n_0, ncol=1)
  
  
  Y_1= S_1n %*% (X_3n+   epsilon_sn)
  
  
  # Huber function and turning parameter  
  c_1=1.65   #1.65   2
  
  Huberf=function(A_1){
    for(r_6 in 1:n_0){ 
      if(abs(A_1[r_6,1])<c_1){  A_1[r_6,1]= A_1[r_6,1]}
      else{A_1[r_6,1]=c_1*sign(A_1[r_6,1])}
    }
    A_1
  }
  
  
  
  # complemented function 
  
  htilde=function(c){2*(c^2)*(1-pnorm(c,mean=0,sd=1))-2*c*dnorm(c,mean=0,sd=1)-1+2*pnorm(c,mean=0,sd=1)}
  
  # ps: test: htilde(0) pnorm(0,mean=0,sd=1)
  
  # alpha matrix function
  D_alpha=diag(n_0)
  alphaf=function(z){
    for(r_9 in 1:n_0 ){ 
      if(z[r_9,1]!=0){z[r_9,1]=Huberf(z)[r_9,1]/(z[r_9,1])} 
      else{z[r_9,1]=1}
      D_alpha[r_9,r_9]=z[r_9,1]
    }
    D_alpha
  }
  
  # G(rho) function
  Gf= function(rhob){
    W_3  %*% solve(I_1- rhob*W_3)
  }
  
  
  
  # error variable
  errorrenew_0=rep(0,infty+1)
  
  # Robust MLE  with repeat method
  

  
  eror=0.00001
  
  
  #  initial value  1
  

 
  #0.3677446
  #0.7
  #0.3
  rhob_0=rep(0,infty+1)
  rhob_max=rep(0,infty+1)
  
  #  initial value 2

  betaf_0=rep(0,infty+1)
  
  #  initial value 3
  sigmafsqueref_0=rep(0,infty+1)

  sigmaf_0= rep(0,infty+1)
  
  #  likelihood function value
  
  
  # terminal variable error
  
  
  pse=1000
  segma=rep(0, pse+1)
  
  
  # commulator 
  n_2=rep(0,1)
  
  # residual summ and mean value
  summ_n1=rep(0,1)
  mean_n1=rep(0,1)
  
  
  rhoerror=rep(0,1)
  likeli_zmax=rep(0,1)
  
  
  ##########################################################################################################################
  
  # use Binary segmentation method  
  n4p=4    # 4 5
  n3=2^(n4p)
  
  # determined lower bound  
  cl=-0.18
  loo=cl
  up=1
  
  # stopping errror (<0.0012 when n3=10)
  er13=1
  er3=0.1
  
  # number
  infty=200
  
  #step1
  rhoop=rep(0,n3)
  betoop=rep(0,n3)
  sigmoop=rep(0,n3)
  rhomoa=rep(0,n3)
  rhoerroroop=rep(0,n3)
  
  #scale
  uppbound=0.95
  scalev=(uppbound-cl)/(n3)
  #############################################################################################################################
  for(r_11 in 1:n3){
    lo=(r_11-1)*scalev+loo
    up= (r_11)*scalev+loo
    
    #Test for fixed initial value
    rhob_0[1]=lo
     
    # S operator
    S_1o= I_1- rhob_0[1]*W_3
    betaf_0[1]=solve(t(X_3n) %*%  X_3n) %*% t(X_3n) %*% S_1o %*% Y_1
    
    #residuals variable
    zf_1=S_1o  %*%  Y_1 - betaf_0[1]* X_3n
    sigmaf_0[1]=(1/n_0)* t(zf_1) %*% zf_1
    
    zf_z3=  (1/ sigmaf_0[1])*(( I_1- rhob_0[1]*W_3)  %*%  Y_1 - betaf_0[1]*X_3n ) 
    zf_zh4= Huberf( zf_z3)
    likeli_zmax[1]=(-n_0)* log(sigmaf_0[1])+log(abs(det(I_1- rhob_0[1]*W_3))) -(1/(2*htilde(c_1)))*t(zf_zh4) %*% zf_zh4
    rhoerror[1]= abs((1/sigmaf_0[1])*t(betaf_0[1]* ( Gf(rhob_0[1]) %*% X_3n) ) %*% zf_zh4 + t(zf_zh4) %*% t(Gf(rhob_0[1])) %*% zf_zh4  -sum(diag(Gf(rhob_0[1])))*htilde(c_1))
    
    ##########################################################################################################################
    
      
      #STEP1 MLE estimates 
      
      for(r_10 in 1:infty){
        # original  
        betaf=betaf_0[r_10]
        
        #beta
        betatest_01[1]=betaf
        for(r_13 in 1:pse){
          zf_0= (1/sigmaf_0[r_10])*(( I_1- rhob_0[r_10]*W_3)  %*%  Y_1 - betatest_01[r_13]* X_3n )
          betatest_01[r_13+1]=solve(t(X_3n) %*%  X_3n) %*% t(X_3n) %*% alphaf(zf_0) %*% ( I_1- rhob_0[r_10]*W_3)  %*% Y_1
          betaf_0[r_10+1]= betatest_01[r_13+1]
          erenew_b0=  abs(betatest_01[r_13+1]-betatest_01[r_13])
          if(erenew_b0 < eror){break}
        }
        
        betaf=betaf_0[r_10+1]
        #residuals variable
        zf_1= (1/(sigmaf_0[r_10]))*(( I_1- rhob_0[r_10]*W_3)  %*%  Y_1 - betaf * X_3n )
        zf_h2=Huberf(zf_1)
        
        #sigma
        sigmafsqueref_0[r_10+1]=((sigmaf_0[r_10])^2/(n_0*htilde(c_1)))*t(zf_h2) %*% zf_h2 
        sigmasqf =  sigmafsqueref_0[r_10+1]
        sigmaf_0[r_10+1]= (sigmasqf)^(1/2)
        sigmaf=sigmaf_0[r_10+1]
        
        #stop
        errorrenew_0[r_10]=abs(sigmaf_0[r_10+1]-sigmaf_0[r_10]) + abs(rhob_0[r_10+1]-rhob_0[r_10])+ abs(betaf_0[r_10+1]-betaf_0[r_10])
        errorrenew=errorrenew_0[r_10]
        if(errorrenew < eror){break}
        
        #residuals variable2 renew
        zf_3= function(rhob){(1/( sigmaf_0[r_10+1] ))*(( I_1- rhob*W_3)  %*%  Y_1 - betaf * X_3n )}
        zf_h4= function(rhob){Huberf( zf_3(rhob))}
        sigmafsqueref_h5=function(rhob){((sigmaf_0[r_10+1])^2/(n_0*htilde(c_1)))*t(zf_h4(rhob)) %*% zf_h4(rhob) }
        
        
        
        MLE_1=function(rhobb){
          #modified logistic likelihood function
          (-n_0/2)* log(sigmafsqueref_h5(rhobb))+log(abs(det(I_1- rhobb*W_3)))-(1/(2*htilde(c_1)))*t(zf_h4(rhobb)) %*% zf_h4(rhobb)
        }  
        
        rhop_10b=optimize(MLE_1, interval = c(lo, up), maximum = TRUE)
        
        rhob_0[r_10+1]=rhop_10b[[1]]
        rhob_ter=rhop_10b[[1]]
        rhob_max[r_10+1]=rhop_10b[[2]]
      }
      
      #STEP3 index of test
      # Residual error
      zf_z3=  (1/ sigmaf )*(( I_1- rhob_ter*W_3)  %*%  Y_1 - betaf * X_3n ) 
      zf_zh4= Huberf( zf_z3)
      zf_zh5=abs(zf_zh4)
      
      # rho minmum 1
      rhoerror[1]= abs((1/sigmaf)*t(betaf* ( Gf(rhob_ter) %*% X_3n) ) %*% zf_zh4 + t(zf_zh4) %*% t(Gf(rhob_ter)) %*% zf_zh4  -sum(diag(Gf(rhob_ter)))*htilde(c_1))
      likeli_zmax[1]=  (-n_0)* log(sigmaf)+log(abs(det(I_1- rhob_ter*W_3))) -(1/(2*htilde(c_1)))*t(zf_zh4) %*% zf_zh4
      # taking value
      betoop[r_11]=betaf
      sigmoop[r_11]=sigmaf
      rhoop[r_11]=rhob_ter
      rhomoa[r_11]=likeli_zmax[1]
      rhoerroroop[r_11]=rhoerror[1]
  
  print(r_11)
  }
  
  ##############################################################################################################################
  
  #global maximum
  
  ngg=which.max(rhomoa)
  rhob_maxg[r_8]=rhomoa[ngg]
  rhob_0g[r_8]=rhoop[ngg]
  betafg[r_8]=betoop[ngg]
  sigmafg[r_8]=sigmoop[ngg]
  rhob_errorg[r_8]=rhoerroroop[ngg]
  
  #give value
  rhob_0gg[r_8]=rhob_0g[r_8]
  rhob_maxgg[r_8]=rhob_maxg[r_8]
  rhob_errorgg[r_8]=rhob_errorg[r_8]
  sigmafgg[r_8]=sigmafg[r_8]
  betafgg[r_8]=betafg[r_8]
  
  ###############################################################################################################################
  
  # strategy for optimum
  
  
  
  # inital position
  nllwo=1
  
  for(r_11 in 1:n4p){ 
    nscale=(n3/(2^(r_11)))
    nlmi1o=nllwo+nscale-1
    nlmi2o=nllwo+nscale
    nlupo=nllwo+2*nscale-1
    bjiao=rep(0,2)
    bjiao1=rep(0,2)
    #1
    n11c=which.max(rhomoa[nllwo:nlmi1o])+nllwo-1
    bjiao1[1]=rhomoa[n11c]
    print(bjiao1[1])
    bjiao[1]=rhoerroroop[n11c]
    print(bjiao[1])
    #2
    n21c=which.max(rhomoa[nlmi2o:nlupo])+nlmi2o-1
    bjiao1[2]=rhomoa[n21c]
    print(bjiao1[2])
    bjiao[2]=rhoerroroop[n21c]
    print(bjiao[2])
    # taking value
    
    
    #previous review
    if(r_11==1){errorrenew1=abs(bjiao1[2]- bjiao1[1])}
    if(errorrenew1 < er13){break}
    errorrenew2=abs(bjiao[2]- bjiao[1])
    if(errorrenew2 < er3){break}
    #choose real value
    if(bjiao[1]<bjiao[2]){
      rhob_0g[r_8]=rhoop[n11c]
      betafg[r_8]=betoop[n11c]
      sigmafg[r_8]=sigmoop[n11c]
      rhob_maxg[r_8]=rhomoa[n11c]
      rhob_errorg[r_8]=rhoerroroop[n11c]
      nllwo=nllwo
    }
    else{
      rhob_0g[r_8]=rhoop[n21c]
      betafg[r_8]=betoop[n21c]
      sigmafg[r_8]=sigmoop[n21c]
      rhob_maxg[r_8]=rhomoa[n21c]
      rhob_errorg[r_8]=rhoerroroop[n21c]
      nllwo= nlmi2o
    }
    print(rhob_0g[r_8])
    # print(betafg[r_8])
    # print(sigmafg[r_8])
    print(rhob_maxg[r_8])
    print(rhob_errorg[r_8])
    
    print('-------------------------------------------------')
    if(errorrenew1 < er13){break}
  }
  
  
  ###############################################################################################################################
  
  
  # tu 1x2
  par(mfrow=c(1,2))
  
  
  plot(rhomoa[1:n3],col = "green", type = "b")
  plot(rhoerroroop[1:n3], col = "blue",type = "b")
  
  print('--------------------------------------------------------------------------------')
  
  
}

# terminal value
#1
rhob_0g
betafg
sigmafg

rhob_maxg
rhob_errorg

#2
rhob_0gg
betafgg
sigmafgg

rhob_maxgg
rhob_errorgg

rhoop
betoop
sigmoop
rhomoa
rhoerroroop

###############################################################################################################################

mean(rhob_0g)-rho_0
mean(betafg)-beta_0
mean(sigmafg)-sigma_0



