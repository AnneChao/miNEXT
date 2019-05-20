#*******************************************************************************
#********************************************************************************
#
## R scipts "Abundance" for Chao et al. (2019) paper on proportional mixture of rarefaction/extrapolation. 
## This R code for abundance data to compute and plot proportional mixture of two rarefaction/extrapolation curves, and
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
sourceCpp('function_abun.cpp')
source("Abundance_subfunction.R")
source("JADE.R")
source("bootstrap_p.R")
source("plot_function.R")

#' Abundance(data, knots = 20, size=NULL, nboots = 0) for abundance data, comupute composite diversity of any sample and species composition (shared and unique species)
#' @param data1 a Sx2 dataframe, the intact assemblage (main) assemblage should be placed in the first column.
#' @param knots the number of points that the mixture diveristy will be computed in rarefaction and extrapolation, respectively. Default is 10.
#' @param size a vector specifying the smapling units of intact (main) assemblage. Default is NULL.
#' @param nbbots the number of replication times for bootstrap. Use 0 to skip bootstrap which might take more time. Default is 0. 
#' @return a list containing 4 tables. Th first 3 are diversities of the two assemblages and the mixed one. The 4th table is the species composition of the mixed assemblage. 
Abundance<-function(data=NULL,knots=10,size=NULL,nboots=0){
  result_abun<-NULL
  result_abun_CI<-NULL
  result_abun<-Abun(data1=data,knots=knots,size=size)
  if(nboots>0){
    print("bootstrap start")
   # print("create bootstrap sample")
    boot_sample<-Abun_CreatBootstrapSample(data=data,nboots=nboots)
    
    boots<-0
    boot_est<-NULL
    boots<-length(boot_sample)
    boot_est <- lapply(1:boots,function(k){
    #  print(paste("boots",k))
      bdata<-boot_sample[[k]]
      out =  Abun(bdata, knots = knots,size=size)
    })
    
   # boot_est_out<<-boot_est
  #  result_abun_out<<-result_abun
    result_abun_CI<-cal_estboot_CI(estboot=boot_est,est=result_abun)
    result_abun_CI
  }
  else{
    result_abun
  }
}

#' Abun(data, knots = 20) for abundance data, comupute composite diversity of any sample and species composition (shared and unique species)
#' @param data1 a Sx2 dataframe, the intact assemblage (main) assemblage should be placed in the first column.
#' @param size a vector specifying the smapling units of intact (main) assemblage. Default is NULL.
#' @param knots the number of points that the mixture diveristy will be computed in rarefaction and extrapolation, respectively. Default is 10.
#' @return a list containing 4 tables. Th first 3 are diversities of the two assemblages and the mixed one. The 4th table is the species composition of the mixed assemblage. 

Abun <- function(data1, knots = 10, size = NULL){
  D1 = sum(data1[, 1]>0)
  D2 = sum(data1[, 2]>0)
  if(D1<D2){
    data1 <- data1[,c(2,1)]
  }else{
    data1 <- data1
  }
  x1 = data1[, 1]
  x2 = data1[, 2]
  n1 = sum(x1)
  n2 = sum(x2)
  D1 = sum(data1[, 1]>0)
  D2 = sum(data1[, 2]>0)
  n = c(n1, n2)
  D = c(D1, D2)
  tm_n = which.max(n)
  tm_D = which.max(D)
  if(tm_D==1){
    if(n1<=n2){
      m2 = round(seq(0, (n1), length.out = knots))
      m1 = n1 - m2
      m = cbind(m1, m2)
      m2 = round(seq(0, (n2), length.out = knots))
    }else{
      m2 = round(c(seq(0, n2, length.out = knots),
                   seq(n2, (n1), length.out = knots)))
      m1 = n1 - m2
      m = cbind(m1, m2)
    }
  }else{
    if(n1<n2){
      m1 = round(c(seq(0, n1, length.out = knots),
                   seq(n1, (n2), length.out = knots)))
      m2 = n2- m1 
      m = cbind(m1, m2)
    }else{
      m1 = round(seq(0, (n2), length.out = knots))
      m2 = m2 - m1 
      m = cbind(m1, m2)
      m1 = round(seq(0, (n2), length.out = knots))
    }
  }
  if(!is.null(size)){
    m1 <- size
    m2 <- n1- m1 
    m = cbind(m1, m2)
  }
  
  maxx1<-max(x1)
  maxx2<-max(x2)
  
  
  ###Assam I,q0,q1,q2############
  
 # print(paste("start  q0_1,q1_1,q2_1",Sys.time())) 
  mm<-cbind(m1,0)
  q0_1<-Dq_in(x1,x2,mm,0)
  q1_1<-Dq_in(x1,x2,mm,1)
  q2_1<-Dq_in(x1,x2,mm,2)
  q2_1_new<-Dq2(x1,x2,mm,n1,n2)
  
  
  ###Assam II,q0 in,ext ok,q1 in ok,ext ok, q2 in ext ok
#  print(paste("start  q0_2_in,q1_2_in",Sys.time())) 
  m2tmp<-m2[m2<=n2]
  mm<-cbind(0,m2tmp)
  q0_2_in<-Dq_in(x1,x2,mm,0)
  q1_2_in<-Dq_in(x1,x2,mm,1)
  
  
  
  p1 = DetAbu(data1[,1], zero = T)
  p2 = DetAbu(data1[,2], zero = T)
  p = cbind(p1, p2)
  p1_hat = p[, 1]
  p2_hat = p[, 2]
  
  m2tmp <-m2[m2>n2]
  m2tmp_s<-m2tmp-n2
  mmext<-cbind(0,m2tmp_s)
  mm<-cbind(0,n2)
  
  q0_2_ext<-Dq_in(x1,x2,mm,0)
  #  h0_2_ext<-Dq0_2_ext(x1,x2,p1_hat,p2_hat,n1,n2,mmext)
#  print(paste("start  q0_2_ext_cpp",Sys.time())) 
  h0_2_ext_cpp<-Dq0_2_ext_cpp(p1_hat,p2_hat,n1,n2,mmext)
  # h0_2_ext_new<-Dq0_2_ext_new(x1,x2,0,p2_hat,n1,n2,mmext)
  q0_2_ext<-q0_2_ext+h0_2_ext_cpp
  
  q0_2<-c(q0_2_in,q0_2_ext)
  
  #q1_2
  mm<-cbind(0,n2)
  q1_2_ext<-Dq_in(x1,x2,mm,1)
  
  m2tmp <-m2[m2>n2]
  mmext<-cbind(0,m2tmp)
#  print(paste("start  q1_2_ext",Sys.time())) 
  h1_2_ext = Dq1_ext(x1,x2,p1_hat,p2_hat,n1,n2,mmext)
  q1_2_ext<-q1_2_ext*h1_2_ext
  q1_2<-c(q1_2_in,q1_2_ext)
  
  #q=2 in and ext using same formular
  mm<-cbind(0,m2)
#  print(paste("start  q2_2",Sys.time())) 
  q2_2<-Dq2(x1,x2,mm,n1,n2)
  
  
  
  
  ####mix q=0,1,2##################################################  
  mm<-cbind(m1,m2)
  
#  print(paste("start  q2_mix",Sys.time())) 
  q2_mix<-Dq2(x1,x2,mm,n1,n2)
  
  m1tmp<-m1[m2<=n2]
  m2tmp<-m2[m2<=n2]
  mm<-cbind(m1tmp,m2tmp)
  
#  print(paste("start  q0_mix_in q1_mix_in ",Sys.time())) 
  q0_mix_in<-Dq_in(x1,x2,mm,0)
  q1_mix_in<-Dq_in(x1,x2,mm,1)
  
  
  m1tmp<-m1[m2>n2]
  m2tmp <-m2[m2>n2]
  m2tmp_s<-m2tmp-n2
  mmexts<-cbind(m1tmp,m2tmp_s)
  mm<-cbind(m1tmp,n2)
  
  q0_mix_ext<-Dq_in(x1,x2,mm,0)
  
#  print(paste("start  h0_mix_ext_cpp  ",Sys.time())) 
  h0_mix_ext_cpp<-Dq0_2_ext_cpp(p1_hat,p2_hat,n1,n2,mmexts)
  # h0_2_ext_new<-Dq0_2_ext_new(x1,x2,0,p2_hat,n1,n2,mmext)
  q0_mix_ext<-q0_mix_ext+h0_mix_ext_cpp
  #NumericVector pi1, NumericVector pi2, int m1, int m2s, int n1, int n2
  q0_mix<-c(q0_mix_in,q0_mix_ext)
  
  m1tmp<-m1[m2>n2]
  m2tmp <-m2[m2>n2]
  mmext<-cbind(m1tmp,m2tmp)
  
  mm<-cbind(m1tmp,n2)
  q1_mix_ext<-Dq_in(x1,x2,mm,1)
 # print(paste("start  h1_mix_ext  ",Sys.time())) 
  h1_mix_ext<-Dq1_ext(x1,x2,p1_hat,p2_hat,n1,n2,mmext)
  q1_mix_ext1<-q1_mix_ext*h1_mix_ext
  q1_mix<-c(q1_mix_in,q1_mix_ext1)
  
  
  #q=2 in and ext using same formular
  mm<-cbind(m1,m2)
  q2_mix<-Dq2(x1,x2,mm,n1,n2)
  
  output0 = cbind(m1, m2, q0_1, q0_2, q0_mix)
  output1 = cbind(m1, m2, q1_1, q1_2, q1_mix)
  output2 = cbind(m1, m2, q2_1, q2_2, q2_mix)
  colnames(output0) = c("m1", "m2", paste(colnames(data1)), "Mixture")
  colnames(output1) = c("m1", "m2", paste(colnames(data1)), "Mixture")
  colnames(output2) = c("m1", "m2", paste(colnames(data1)), "Mixture")
  
  ###########q0_ana: # q=0 unique and share#############yhc code#######
  # D0 = round(D0, 3)
  print(paste("q0_ana  ",Sys.time())) 
  m1_i = unique(m1[m1 >= n1-n2])
  m2_i = unique(n1 - m1_i)
  datash = data1[(data1[,1]>0 & data1[,2]>0), , drop=F]
  dataun1 = data1[(data1[,1]>0 & data1[,2]==0), , drop=F]
  dataun2 = data1[(data1[,1]==0 & data1[,2]>0), , drop=F]
  un1 = sapply(m1_i,function(i){
    sum(un_abun(xi = dataun1[,1], n = n1, m = i))
  })
  un2 = sapply(m2_i,function(i){
    sum(un_abun(xi = dataun2[,2], n = n2, m = i))
  })
  sh12 = sapply(1:length(m1_i),function(i){
    sum(sh_abun(xi1 = datash[,1],xi2 = datash[,2], n1 = n1,m1 =  m1_i[i],
                n2 = n2, m2 = m2_i[i]))
  })
  # sh = sh1 + sh2 -sh12
  # total = un1 + un2 + sh
  # n1 n2 need adjustment
  
  ##########prof. Chao's new paper doesn't mention it################################################
  # if(max(m2)>n2){
  #   m1_e <- unique(m1[m1 < n1-n2])
  #   m2_e <- unique(n1 - m1_e)
  #   psh1 <-  p[(data1[,1]>0 & data1[,2]>0) , 1]
  #   psh2 <-  p[(data1[,1]>0 & data1[,2]>0) , 2]
  #   pun1 <-  p[(data1[,1]>0 & data1[,2]==0) , 1]
  #   pun2 <-  p[(data1[,1]==0 & data1[,2]>0) , 2]
  #   ex <- sapply(1:length(m1_e), function(i){
  #     #sh12_ex_1 <- (1- exp(lchoose(n1 - datash[,1], m1_e[i]) - lchoose(n1, m1_e[i])) ) %>% sum
  #     sh12_ex_2 <- ( (1-psh1)^m1_e[i] * (1-psh2)^n2 * (1-(1-psh2)^(m2_e[i]-n2))  / (1-(1-psh1)^n1) / (1-(1-psh2)^n2) ) %>% sum()
  #     sh12_ex <- length(psh1) + sh12_ex_2
  #     
  #     un1_ex <- (1- exp(lchoose(n1 - dataun1[,1], m1_e[i])-lchoose(n1, m1_e[i])) ) %>% sum
  #     
  #     un2_ex <- ( ((1-pun2)^n2 - (1-pun2)^m2_e[i])  / (1-(1-pun2)^n2) )%>% sum() + length(pun2)
  #     
  #     rbind(un1_ex, un2_ex, sh12_ex)
  #     
  #   })
  #   un1 <- c(un1, un1[length(sh12)], ex[1,])
  #   un2 <- c(un2, un2[length(sh12)], ex[2,])
  #   sh12 <- c(sh12, sh12[length(sh12)], ex[3,])
  #   m1_i <- m1
  #   m2_i <- m2
  # }
  q0_ana = data.frame(m1 = m1_i, m2 = m2_i,
                      q0_un1 = un1, q0_un2 = un2,
                      q0_sh = sh12)
  colnames(q0_ana)[3:5] = c(paste("Unique to",colnames(data1)[1]),
                            paste("Unique to",colnames(data1)[2]),
                            "Share")  
  
  
  list(q0 = output0, q1 = output1, q2 = output2, q0_ana = q0_ana)
  
  
  
}


####################################################################################
#
# (2). Example
#
####################################################################################
spider = read.table("Spider_Abundance_Data.txt")
result_abun_spider = Abundance(data=spider, knots = 10)

multi.plot(data = spider, ans1 = result_abun_spider, type = "abundance")

