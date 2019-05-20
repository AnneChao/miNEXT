#*******************************************************************************
#********************************************************************************
#
## R scipts "Incidence" for Chao et al. (2019) paper on proportional mixture of rarefaction/extrapolation. 
## This R code for incidence data to compute and plot proportional mixture of two rarefaction/extrapolation curves, and
# the number of shared and unique species in any two mixture of rarefaction curves.
# This code includes three parts:
# (1) Comupute composite diversity of any sample and species composition (shared and unique species);
# (2) Plot individual rarefaction/extrapolation curves and mixture curves;
# (3) Example.
#
# NOTE: The packages "ggplot2", "dplyr", "reshape2", "ggpubr", "Rcpp", must be 
# installed and loaded before running the scripts. 
# 
#
#*******************************************************************************
#*******************************************************************************


####################################################################################
#
# (1). Comupute composite diversity of any sample and species composition (shared and unique species)
#
####################################################################################

library(Rcpp)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(reshape2)
sourceCpp('function_inci.cpp')
source('JADE.r')
source('bootstrap_p.r')
source('plot_function.r')

#' Incidence(data, allpts = FALSE, size = NULL, knots = 20 ,nboots = 0) for incidence data, comupute composite diversity of any sample and species composition (shared and unique species)
#' along with their confidence interval.
#' @param data a Sx2 dataframe, the intact assemblage (main) assemblage should be the first column.
#' @param allpts specifying whether to compute all combinations of sampling units of two assemblages. Default is FALSE.
#' @param size a vector specifying the smapling units of intact (main) assemblage. Default is NULL.
#' @param knots the number of points that the mixture diveristy will be computed. Default is 20.
#' @param nboots the number of replication bootstrap times. Use 0 to skip bootstrap which might take more time. Default is 0.                 
#' @return a list containing 4 tables. Th first 3 are diversities of the two assemblages and the mixed one. The 4th table is the species composition of the mixed assemblage. All estimators
#' are presented along with their confidence interval.
Incidence <- function(data, allpts = FALSE, size = NULL, knots = 20, nboots = 0){
  data <- data[rowSums(data)>0,]
  if(sum(data[,1]>0) < sum(data[,2]>0)){
    data = data[,c(2,1)]
  }
  esti <- Inci(data, allpts, size, knots )
  if(nboots>1){
    nT1 <- data[1,1]
    nT2 <- data[1,2]
    datap = data[-1,]

    p1_est = boot_p_inc(data[,1])
    p2_est = boot_p_inc(data[,2])
    # p1_est = c(DetInc(data[,1],zero = T),UndInc(data[,1]))
    # p2_est = c(DetInc(data[,2],zero = T),UndInc(data[,2]))
    datap[,1]<-p1_est[1:nrow(datap)]
    datap[,2]<-p2_est[1:nrow(datap)]
    undetec1 <- p1_est[-c(1:nrow(datap))]
    undetec2 <- p2_est[-c(1:nrow(datap))]
    
    Q1 = sum(rowSums(data[-1,])==1)
    Q2 = sum(rowSums(data[-1,])==2)
    nT = nT1+nT2
    Q0.hat <- ceiling(ifelse(Q2 == 0, (nT - 1) / nT * Q1 * (Q1 - 1) / 2, (nT - 1) / nT * Q1 ^ 2/ 2 / Q2))  #Chao2 estimator
    datap = matrix(data = 0,nrow = Q0.hat, ncol = 2,
                   dimnames = list(NULL, names(datap))) %>% rbind(datap,.)
    zero1 <-  which(datap[,1]==0)
    zero2 <-  which(datap[,2]==0)
    data_boot <- sapply(1:nboots,function(k){
      fill1 <- sample(x = zero1,size = length(undetec1), replace = F)
      fill2 <- sample(x = zero2,size = length(undetec2), replace = F)
      datapp <- datap 
      datapp[fill1,1] <- undetec1
      datapp[fill2,2] <- undetec2
      data1 <- sapply(datapp[,1], function(i) rbinom(n = 1, size = nT1, prob = i)) 
      data2 <- sapply(datapp[,2], function(i) rbinom(n = 1, size = nT2, prob = i)) 
      tmp <- sum(data1>0) > sum(data2>0)
      tmp2 <- ((data1==0) & (data2>0)) %>% sum 
      while( (tmp==F) | tmp2==0 ){
        fill1 <- sample(x = zero1,size = length(undetec1), replace = F)
        fill2 <- sample(x = zero2,size = length(undetec2), replace = F)
        datapp <- datap 
        datapp[fill1,1] <- undetec1
        datapp[fill2,2] <- undetec2
        data1 <- sapply(datapp[,1], function(i) rbinom(n = 1, size = nT1, prob = i)) 
        data2 <- sapply(datapp[,2], function(i) rbinom(n = 1, size = nT2, prob = i)) 
        tmp <- sum(data1>0) > sum(data2>0)
        tmp2 <- ((data1==0) & (data2>0)) %>% sum 
      }
      c(data1,data2)
    })
    sepe <- nrow(data_boot)
    out_boot <- sapply(1:nboots, function(j){
      #print(j)
      d1 <- c(nT1,data_boot[1:(sepe/2),j])
      d2 <- c(nT2,data_boot[(sepe/2+1):sepe,j])
      data_b <- matrix( c(d1,d2), ncol = 2 ,dimnames = list(NULL, names(datap) ))
      data_b <- data_b[rowSums(data_b)>0,] 
      out = Inci(data_b, allpts, size, knots )
      out012 <- matrix(c(out$q0[,5], out$q1[,5], out$q2[,5]),
                       ncol = 3, dimnames = list(NULL,c("q0","q1","q2")))
      out_p <- matrix(c(out$q0_ana[,3], out$q0_ana[,4], out$q0_ana[,5]),
                      ncol = 3, dimnames = list(NULL,colnames(out$q0_ana)[3:5])) %>% rbind(., .[1:(nrow(out$q0)-nrow(out$q0_ana)),])
      cbind(out012, out_p)
    }, simplify = "array")
    sd_boot <- apply(out_boot,MARGIN = c(1,2), sd) 
    esti$q0 <- esti$q0 %>% cbind(., LCL = (.[,5] - 1.96*sd_boot[,1]) , UCL = (.[,5] + 1.96*sd_boot[,1]), s.e. = sd_boot[,1])
    esti$q1 <- esti$q1 %>% cbind(., LCL = (.[,5] - 1.96*sd_boot[,2]) , UCL = (.[,5] + 1.96*sd_boot[,2]), s.e. = sd_boot[,2])
    esti$q2 <- esti$q2 %>% cbind(., LCL = (.[,5] - 1.96*sd_boot[,3]) , UCL = (.[,5] + 1.96*sd_boot[,3]), s.e. = sd_boot[,3])
    tmp <- nrow(esti$q0_ana)
    esti$q0_ana <- esti$q0_ana %>% cbind(., LCL = (.[1:tmp,3:5] - 1.96*sd_boot[1:tmp,4:6]) , UCL = (.[1:tmp,3:5] + 1.96*sd_boot[1:tmp,4:6]),
                                         s.e. = sd_boot[1:tmp,4:6])
    esti$q0_ana[ esti$q0_ana < 0 ] = 0
    esti
  }else{
    esti
  }
}

#' Inci(data, allpts = FALSE, size = NULL, knots = 20 ) for incidence data, comupute composite diversity of any sample and species composition (shared and unique species)
#' @param data a Sx2 dataframe, the intact assemblage (main) assemblage should be the first column.
#' @param allpts specifying whether to compute all combinations of sampling units of two assemblages. Default is FALSE.
#' @param size a vector specifying the smapling units of intact (main) assemblage. Default is NULL.
#' @param knots the number of points that the mixture diveristy will be computed. Default is 20.
#' @return a list containing 4 tables. Th first 3 are diversities of the two assemblages and the mixed one. The 4th table is the species composition of the mixed assemblage. 
Inci <- function(data, allpts = FALSE, size = NULL, knots = 20 ){
  d2_hat  <- function(x1,x2,t1,t2){
    if(t1!=0 || t2!=0){
      u1 = sum(x1)
      u2 = sum(x2)
      u1_hat = t1*u1/T1;
      u2_hat = t2*u2/T2;
      D2 = 1/(u1_hat+u2_hat)+sum(x1*(x1-1)/(T1*(T1-1)))*t1*(t1-1)/(u1_hat+u2_hat)^2+sum(x2*(x2-1)/(T2*(T2-1)))*t2*(t2-1)/(u1_hat+u2_hat)^2+2*sum(x1*x2)*t1*t2/(T1*T2)/(u1_hat+u2_hat)^2
      1/D2
    }
    else if (t1==0 && t2==0){
      D2 = 0
      D2
    }
    
  }
  
  ts = as.numeric(data[1,])
  T1 = ts[1]
  T2 = ts[2]
  
  d0ob = c(sum(data[,1]>0),sum(data[,2]>0))
  if(T1 == T2){
    extr = 0
    if(d0ob[1]<d0ob[2]){
      data = data[,c(2,1)]
    }
  }else if(((T1-T2)*(d0ob[1]-d0ob[2])) > 0){
    extr = 1
    if(T1<T2){
      data = data[,c(2,1)]
    }
  }else if (((T1-T2) * (d0ob[1]-d0ob[2])) < 0 ){
    extr = 0
    if(T1>T2){
      data = data[,c(2,1)]
    }
  }
  ts = as.numeric(data[1,])
  T1 = ts[1]
  T2 = ts[2]
  if(allpts == TRUE){
    t1 = seq(T1,0)
    t2 = T1-t1
  }else if (allpts == FALSE & !is.null(size)){
    if(max(size)>T1)
      stop(paste0("Because the community ",colnames(data)[1]," has higher divrsity, its sampling units should be viewed main."))
    
    t1 = sort(c(size,(T1-T2),T1),decreasing = T)
    t1 = t1[!duplicated(t1[t1>=0])]
    t2 = T1-t1
  }else if (allpts == FALSE & is.null(size)){
    t1 = round(seq(T1,0,length.out = knots),0)
    t1 = sort(c(t1,(T1-T2)),decreasing = T)
    t1 = t1[t1>=0] %>% unique
    t2 = T1-t1
  }
  
  datap = data[-1,]
  datap[,1] = DetInc(data[,1],zero = T)
  datap[,2] = DetInc(data[,2],zero = T)
  
  data = data[-1,]
  D = c(sum(data[,1]>0),sum(data[,2]>0))
  t1_in = sort(t1[t1>=(T1-T2)],decreasing = T) 
  t2_in = t1_in[1] - t1_in
  d01_int = t(sapply(1:length(t1_in),function(i) D0_rare(xi = c(T1,data[,1]),yi = c(T2,data[,2]),t1 = t1_in[i],t2 = t2_in[i])))
  
  d01 = d01_int
  sites = colnames(data)
  datash = data[(data[,1]>0 & data[,2]>0), , drop = F]
  dataun1 = data[(data[,1]>0 & data[,2]==0), , drop = F]
  dataun2 = data[(data[,1]==0 & data[,2]>0), , drop = F]
  
  un1 = sapply(t1_in,function(i){
    sum(un_inci(yi = dataun1[,1], T = T1, t = i))
  })
  un2 = sapply(t2_in,function(i){
    sum(un_inci(yi = dataun2[,2], T = T2, t = i))
  })
  sh12 = sapply(1:length(t1_in),function(i){
    sum(sh_inci(yi1 = datash[,1],yi2 = datash[,2], T1 = T1, t1 =  t1_in[i],
                T2 = T2, t2 = t2_in[i]))
  })
  
  q0_ana = data.frame(m1 = t1_in, m2 = t2_in,
                      q0_un1 = un1, q0_un2 = un2,
                      q0_sh = sh12)
  
  
  colnames(q0_ana)[3:5] = c(paste0("Uni.",sites[1]),
                            paste0("Uni.",sites[2]),
                            "Share")
  
  
  #assemblage I
  a1_q01 =  t(sapply(1:length(t1[t1<=T1]),function(i) D0_rare(xi = c(T1,data[,1]),yi = c(T2,data[,2]),t1 = t1[t1<=T1][i],t2 = 0))) 
  a1_q2 =  sapply(1:length(t1),function(i) d2_hat(x1 = data[,1],x2 = data[,2],sort(t1,decreasing = T)[i],0))
  
  
  if(extr == 1){ 
    t1_ex = sort(t1[t1<(T1-T2)],decreasing = T) 
    t2_ex = T1-t1_ex
    
    d01_ext = t(sapply(1:length(t1_ex),function(i) D0_rare(xi = c(T1,data[,1]),yi = c(T2,data[,2]),t1 = t1_ex[i],t2 = T2)))
    d0_ext_2 = sapply(1:length(t1_ex),function(i) h0_hat(datap[,1],datap[,2],t1_ex[i],t2_ex[i],T1,T2))
    d1_ext_2 = sapply(1:length(t1_ex),function(i) h1_hat(datap[,1],datap[,2],data[,1],data[,2],t1_ex[i],t2_ex[i],T1,T2))
    d1_ext_2 = exp(d1_ext_2)
    d01_ext[,1] = d01_ext[,1] + d0_ext_2
    d01_ext[,2] = d01_ext[,2] * d1_ext_2
    
    
    d01 = rbind(d01_int,d01_ext)
    
    #assemblage II extrapolation
    a2_q01_ext = t(sapply(1:length(t1_ex),function(i) D0_rare(xi = c(T1,data[,1]),yi = c(T2,data[,2]),t1 =0,t2 = T2)))
    a2_q0_ext_2 = t(sapply(1:length(t2_ex),function(i) h0_hat(datap[,1],datap[,2],0,t2_ex[i],T1,T2)))[,1]
    a2_q1_ext_2 = t(sapply(1:length(t1_ex),function(i) h1_hat(datap[,1],datap[,2],data[,1],data[,2],0,t2_ex[i],T1,T2)))[,1]
    a2_q1_ext_2 = exp(a2_q1_ext_2)
    a2_q01_ext[,1] = a2_q01_ext[,1] + a2_q0_ext_2
    a2_q01_ext[,2] = a2_q01_ext[,2] * a2_q1_ext_2
  }
  d2_in_ex = sapply(1:length(t1),function(i) d2_hat(x1 = data[,1],x2 = data[,2],sort(t1,decreasing = T)[i],sort(t2)[i]))
  d012 = cbind(t1 = t1, d01,d2_in_ex)
  
  #assemblage II interpolation
  if(T2 > T1) {
    t2_exceed = seq(T1,T2,length.out = 5) %>% round 
    t2_in = c(t2_in, t2_exceed) %>% unique
    t2 = c(t2, t2_exceed) %>% unique
  }
  a2_q01 = t(sapply(1:length(t2_in),function(i) D0_rare(xi = c(T1,data[,1]),yi = c(T2,data[,2]),t1 = 0,t2 = t2_in[i])))
  a2_q2 =  sapply(1:length(t2),
                  function(i) d2_hat(x1 = data[,1],x2 = data[,2],0,t2[i])) 
  if(extr == 1){ a2_q01 = rbind(a2_q01,a2_q01_ext) }
  
  
  d0 = cbind(t1 = t1, t2 = t2, a1_q01[,1], a2_q01[,1], d012[,2]);colnames(d0) = c("t1","t2",colnames(data)[1],colnames(data)[2],"Mixture")
  d1 = cbind(t1 = t1, t2 = t2, a1_q01[,2], a2_q01[,2], d012[,3]);colnames(d1) = c("t1","t2",colnames(data)[1],colnames(data)[2],"Mixture")
  d2 = cbind(t1 = t1, t2 = t2, a1_q2, a2_q2, d012[,4]);colnames(d2) = c("t1","t2",colnames(data)[1],colnames(data)[2],"Mixture")
  
  if(extr == 1){ 
    cross = which(d0[,2] == T2)
    d0 = rbind(d0,d0[cross,]);d0 = d0[order(d0[,2]),];rownames(d0) = NULL
    d1 = rbind(d1,d1[cross,]);d1 = d1[order(d1[,2]),];rownames(d1) = NULL
    d2 = rbind(d2,d2[cross,]);d2 = d2[order(d2[,2]),];rownames(d2) = NULL
    # q0_ana = rbind(q0_ana,q0_ana[cross,]); q0_ana = q0_ana[order(q0_ana[,2]),]
    # rownames(q0_ana) = NULL
    
  }
  
  
  
  colnames(d0)[c(1,2)] = c("m1","m2")
  
  colnames(d1)[c(1,2)] = c("m1","m2")
  colnames(d2)[c(1,2)] = c("m1","m2")
  
  if(T2 > T1){
    d0[d0[,2]>T1,c(3,5)]=NA
    d1[d1[,2]>T1,c(3,5)]=NA
    d2[d2[,2]>T1,c(3,5)]=NA
  }
  return(list(q0 = d0, q1 = d1, q2 = d2, q0_ana = q0_ana))
}
####################################################################################
#
# (2). Example
#
####################################################################################
bird = read.table("Bird_Incidence_Data.txt")
result_inci_bird = Incidence(bird, knots = 10, nboots = 50)
multi.plot(data = bird, ans1 = result_inci_bird, type = "incidence")
