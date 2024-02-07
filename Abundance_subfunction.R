

Dq_in<-function(x1,x2,mm,q){
  out = apply(mm,1,function(mm) {
    mm1 = mm[1]
    mm2 = mm[2]
    D_share(x1, x2, mm1, mm2, q)
  })
}



Dq0_2_ext_cpp<-function(p1_hat,p2_hat,n1,n2,mmext){
  out = apply(mmext,1,function(mm) {
    mm1 = mm[1]
    mm2s = mm[2]
    h0hat<-h0_hat_cpp(p1_hat, p2_hat, mm1, mm2s,n1,n2)
  })
}


Dq0_mix_ext_cpp<-function(p1_hat,p2_hat,n1,n2,mmext){
  out = apply(mmext,1,function(mm) {
    mm1 = mm[1]
    mm2s = mm[2]
    h0hat<-h0_hat_cpp(p1_hat, p2_hat, mm1, mm2s,n1,n2)
  })
}

h0hat<-function( pi1, pi2, m1,  m2s, n2 ){
  output = 0
  output = (1-pi1)^m1*(1-pi2)^n2*(1-(1-pi2)^m2s)
  return (output)
}


Dq0_2_ext_old<-function(x1,x2,p1_hat,p2_hat,n1,n2,mmext){
  out = apply(mmext,1,function(mm) {
    mm1 = mm[1]
    mm2s = mm[2]
    
    p1a = p1_hat[(x1>=1)&(x2>=1)]
    p2a = p2_hat[(x1>=1)&(x2>=1)]
    p1b = p1_hat[(x1==0)&(x2>=1)]
    p2b = p2_hat[(x1==0)&(x2>=1)]
    p2a = p2_hat[(x2>=1)]
    sub1 = sum((1-p1a)^mm1*(1-p2a)^n2*(1-(1-p2a)^(mm2s))/((1-(1-p1a)^n1)*(1-(1-p2a)^n2)))
    # sub1 = sum(h0hat(p1a,p2a,mm1,mm2s,n2)/((1-(1-p1a)^n1)*(1-(1-p2a)^n2)))
  })
}


Dq0_2_ext<-function(x1,x2,p1_hat,p2_hat,n1,n2,mmext){
  out = apply(mmext,1,function(mm) {
    mm1 = mm[1]
    mm2s = mm[2]
    
    p1a = p1_hat[(x1>=1)&(x2>=1)]
    p2a = p2_hat[(x1>=1)&(x2>=1)]
    p1b = p1_hat[(x1==0)&(x2>=1)]
    p2b = p2_hat[(x1==0)&(x2>=1)]
    #p2a = p2_hat[(x2>=1)]
    p1anew<<-p1a
    p2anew<<-p2a
    sub1=sum((1-p1a)^mm1*(1-p2a)^n2*(1-(1-p2a)^(mm2s))/((1-(1-p1a)^n1)*(1-(1-p2a)^n2)))
    
    # sub1 = sum((1-p1a)^mm1*(1-p2a)^n2*(1-(1-p2a)^(mm2s))/((1-(1-p1a)^n1)*(1-(1-p2a)^n2)))
    # sub1 = sum(h0hat(p1a,p2a,mm1,mm2s,n2)/((1-(1-p1a)^n1)*(1-(1-p2a)^n2)))
  })
}


Dq0_2_ext_new<-function(x1,x2,p1_hat,p2_hat,n1,n2,mmext){
  out = apply(mmext,1,function(mm) {
    mm1 = mm[1]
    mm2s = mm[2]
    
    p1a = p1_hat[(x1>=1)&(x2>=1)]
    p2a = p2_hat[(x1>=1)&(x2>=1)]
    p1b = p1_hat[(x1==0)&(x2>=1)]
    p2b = p2_hat[(x1==0)&(x2>=1)]
    p2a = p2_hat[(x2>=1)]
    
    p1anew<<-p1a
    p2anew<<-p2a
    sub1=sum((1-p2a)^n2*(1-(1-p2a)^(mm2s))/((1-(1-p2a)^n2)))
    
    # sub1 = sum((1-p1a)^mm1*(1-p2a)^n2*(1-(1-p2a)^(mm2s))/((1-(1-p1a)^n1)*(1-(1-p2a)^n2)))
    # sub1 = sum(h0hat(p1a,p2a,mm1,mm2s,n2)/((1-(1-p1a)^n1)*(1-(1-p2a)^n2)))
  })
}



Dq2<-function(x1,x2,mm,n1,n2){
  
  out = apply(mm,1,function(mm) {
    m1 = mm[1]
    m2 = mm[2]
    
    if(m1 !=0 || m2!=0){
      
      D2 = 1/(m1+m2)+sum(x1*(x1-1)/(n1*(n1-1)))*m1*(m1-1)/(m1+m2)^2+sum(x2*(x2-1)/(n2*(n2-1)))*m2*(m2-1)/(m1+m2)^2+2*sum(x1*x2/(n1*n2))*m1*m2/(m1+m2)^2
      D2<-1/D2
      
    }
    else if(m1==0 && m2==0)
      D2<-0
  })
}




Dq1_ext<-function(x1,x2,p1_hat,p2_hat,n1,n2,mm){
  
  out = apply(mm,1,function(mm) {
    m1 = mm[1]
    m2 = mm[2]
    # print(paste("m1=",m1,"m2=",m2,Sys.time()))
    h1hat<-h1_hat(p1_hat, p2_hat, x1,x2,m1, m2,n1,n2)
    # print(paste("h1hat=",h1hat,Sys.time()))
    dq1_ext<-exp(h1hat)
  })
}


Dq1_ext_choose<-function(x1,x2,p1_hat,p2_hat,n1,n2,mm){
  maxx1<-max(x1>0)
  maxx2<-max(x2>0)
  out = apply(mm,1,function(mm) {
    m1 = mm[1]
    m2 = mm[2]
    loopm1<-min(m1,maxx1)
    loopm2<-min(m2,maxx2)
    h1hat<-h1_hat_choose(p1_hat, p2_hat, x1,x2,m1, m2,n1,n2)
    print(paste("h1_hat_choose=",h1hat,Sys.time()))
    dq1_ext<-exp(h1hat)
  })
}



Chat1_f0Fun <-function(f1, f2, n) {
  if (f2 > 0) {
    f0 <- (n - 1) / n * f1^2 / (2 * f2)
    #C<-1-f0*f1/(n*f0+f1)
    #C <- 1 - f1 / n * ((n - 1) * f1 / ((n - 1) * f1 + 2 * f2))
    C <- 1 - f1 / n * (n-1)*f1/((n-1)*f1+2*f2)
    
  } else if (f2 == 0 & f1 != 0) {
    f0 <- (n - 1) / n * f1 * (f1 - 1) / 2
    #C<-1-f0*f1/(n*f0+f1)
    C <- 1 - f1 / n * ((n - 1) * (f1 - 1) / ((n - 1) * (f1 - 1) + 2))
    
  } else {
    f0 <- (n - 1) / n * f1 * (f1 - 1) / 2
    #f0 <- 0
    C <- 1
  }
  f0 <- ceiling(f0)
  return(c(C, f0))
}


##for create bootstrap sample#####
Abun_CreatBootstrapSample <- function(data, nboots = 0){
  data <- data[rowSums(data)>0,]
  if(sum(data[,1]>0) < sum(data[,2]>0)){
    data = data[,c(2,1)]
  }
  
  if(nboots>1){
    x1 = data[, 1]
    x2 = data[, 2]
    n1 = sum(x1)
    n2 = sum(x2)
    
   
    ##bootstrap p from JADE: two parameters
    # p1 = DetAbu(x1, zero = T)
    # p2 = DetAbu(x2, zero = T)
    # data_p = data
    # data_p[,1]<-p1
    # data_p[,2]<-p2
    # undetec1 <- UndAbu(x1)
    # undetec2 <- UndAbu(x2)
    
    ## bootstrap p from 2013: one parameter
    p1_est = boot_p_abu(x1)
    p2_est = boot_p_abu(x2)
    data_p = data
    data_p[,1]<-p1_est[1:nrow(data_p)]
    data_p[,2]<-p2_est[1:nrow(data_p)]
    undetec1 <- p1_est[-c(1:nrow(data_p))]
    undetec2 <- p2_est[-c(1:nrow(data_p))]
    
    
    
    f1 = sum(rowSums(data)==1)
    f2 = sum(rowSums(data)==2)
    n=n1+n2
    f0hat<-Chat1_f0Fun(f1,f2,n)[2]
    
    data_p = matrix(data = 0,nrow = f0hat, ncol = 2,
                    dimnames = list(NULL, names(data_p))) %>% rbind(data_p,.)
    zero1 <-  which(data_p[,1]==0)
    zero2 <-  which(data_p[,2]==0)
    
    data_boot <- list() 
    for(k in 1:nboots){
      # fill1 <- sample(x = zero1,size = length(undetec1), replace = F)
      # fill2 <- sample(x = zero2,size = length(undetec2), replace = F)
      if (length(zero1) >= 1) fill1 <- sample(x = zero1,size = length(undetec1), replace = F) else fill1 <- NULL   ## 2024/02/07 modify
      if (length(zero2) >= 1) fill2 <- sample(x = zero2,size = length(undetec2), replace = F) else fill2 <- NULL   ## 2024/02/07 modify
      
      data_pboot <- data_p 
      data_pboot[fill1,1] <- undetec1
      data_pboot[fill2,2] <- undetec2
      nrow(data_pboot)
      # bootx1 <- sapply(data_pboot[,1], function(i) rbinom(n = 1, size = n1, prob = i)) 
      # bootx2 <- sapply(data_pboot[,2], function(i) rbinom(n = 1, size = n2, prob = i)) 
      bootx1<- rmultinom (n = 1, size = n1, prob = data_pboot[,1])
      bootx2<- rmultinom (n = 1, size = n2, prob = data_pboot[,2])
      
      tmp <- sum(bootx1>0) > sum(bootx2>0)
      tmp2 <- ((bootx1==0) & (bootx1>0)) %>% sum 
      while( (tmp==F) | tmp2==0 ){
        # fill1 <- sample(x = zero1,size = length(undetec1), replace = F)
        # fill2 <- sample(x = zero2,size = length(undetec2), replace = F)
        if (length(zero1) >= 1) fill1 <- sample(x = zero1,size = length(undetec1), replace = F) else fill1 <- NULL   ## 2024/02/07 modify
        if (length(zero2) >= 1) fill2 <- sample(x = zero2,size = length(undetec2), replace = F) else fill2 <- NULL   ## 2024/02/07 modify
        
        data_pboot <- data_p 
        data_pboot[fill1,1] <- undetec1
        data_pboot[fill2,2] <- undetec2
        bootx1<- rmultinom (n = 1, size = n1, prob = data_pboot[,1])
        bootx2<- rmultinom (n = 1, size = n2, prob = data_pboot[,2]) 
        tmp <- sum(bootx1>0) > sum(bootx2>0)
        tmp2 <- ((bootx1==0) & (bootx2>0)) %>% sum 
      }
      data_b<-cbind(bootx1,bootx2)
      data_b <- data_b[rowSums(data_b)>0,] 
      colnames(data_b)<-c(colnames(data))
      data_boot[[k]]<-data_b
      
    }
    
  }
  return(data_boot)
}


###for calculate abundance bootstrap confidence interval####

cal_estboot_CI<-function(estboot=NULL,est=NULL){
  out_boot012 <- sapply(1:length(estboot), function(j){
    tmp<-estboot[[j]]
    out012tmp<-matrix(c(tmp$q0[,5], tmp$q1[,5], tmp$q2[,5]), ncol = 3,dimnames = list(NULL,c("q0","q1","q2")))
    
  }, simplify = "array")
  
  
  out_bootana <- sapply(1:length(estboot), function(j){
    tmp<-estboot[[j]]
    outanatmp<-matrix(c(tmp$q0_ana[,3], tmp$q0_ana[,4], tmp$q0_ana[,5]),ncol = 3,dimnames=list(NULL,colnames(tmp$q0_ana)[3:5]))
    
  }, simplify = "array")
  sd_boot_012 <- apply(out_boot012,MARGIN = c(1,2), sd) 
  sd_boot_ana <- apply(out_bootana,MARGIN = c(1,2), sd) 
  
  # mean_boot_012 <- apply(out_boot012,MARGIN = c(1,2), mean) 
  #  mean_boot_ana <- apply(out_bootana,MARGIN = c(1,2), mean) 
  estfinal<-est
  estfinal$q0 <- estfinal$q0 %>% cbind(., LCL = (.[,5] - 1.96*sd_boot_012[,1]) , UCL = (.[,5] + 1.96*sd_boot_012[,1]),s.e.=sd_boot_012[,1])
  estfinal$q1 <- estfinal$q1 %>% cbind(., LCL = (.[,5] - 1.96*sd_boot_012[,2]) , UCL = (.[,5] + 1.96*sd_boot_012[,2]),s.e.=sd_boot_012[,2])
  estfinal$q2 <- estfinal$q2 %>% cbind(., LCL = (.[,5] - 1.96*sd_boot_012[,3]) , UCL = (.[,5] + 1.96*sd_boot_012[,3]),s.e.=sd_boot_012[,3])
  estfinal$q0_ana <- estfinal$q0_ana %>% cbind(., LCL = (.[,3:5] - 1.96*sd_boot_ana[,1:3]) , UCL = (.[,3:5] + 1.96*sd_boot_ana[,1:3]),s.e=sd_boot_ana[,1:3])
  estfinal$q0_ana[ estfinal$q0_ana < 0 ] = 0
  
  return(estfinal)
  
}



