#' @title Calculations for some useful statistics 
#' @name tableCalculation
#' @description Calculate the statistic which is used to test difference between two survival curves
#' @param df the dataframe which consists of three columns:surviving time ,censoring or not(0,1),group
#' @return a dataframe consisists of ten columns,representing:"death occuring time","deathnumber at each time","total number at risk","expected death number in group1","expected death number in group2","number at risk in group1","number at risk in group2","variance of theoretical death number","death number in group1","death number in group2".
#' @examples
#' \dontrun{
#' tableCalculation(aml1)
#' }
#' @import graphics
#' @import stats
#' @import knitr
#' @import rmarkdown
#' @export

tableCalculation<-function(df){
  
  group1=as.vector(unique(df[,3]))                  
  
  n1=length(as.vector(df[which(df[,3]==group1[1]),1]))
  n2=length(as.vector(df[which(df[,3]==group1[2]),1]))    
  
  new.df<-data.frame(matrix(0,n1+n2,17))                 
  indices<-order(df[,1])                                              
  
  new.df[1]<-df[indices,3]
  new.df[2]<-df[indices,1]
  new.df[3]<-df[indices,2]
  new.df[16]<-df[indices,2]
  names(new.df)<-c("group","t","deathnumber","njia","nyi","total","jiatheory","yitheory","njiaji","nyiji","v","dnjia","dnyi","dnjiaji","dnyiji","dnpie","indx")
  new.df$njia[1]=n1
  new.df$nyi[1]=n2
  new.df$total=c((n1+n2):1)
  new.df$njiaji=n1
  new.df$nyiji=n2
  jiafacttotal<-sum(subset(df,df[,3]==group1[1])[,2])
  yifacttotal<-sum(subset(df,df[,3]==group1[2])[,2])
  for(i in 1:(n1+n2)){
    if(new.df$dnpie[i]!=0){
      if(new.df$group[i]==group1[1]){
        new.df$dnjia[i]=new.df$dnjia[i]+1
      }
      if(new.df$group[i]==group1[2]){
        new.df$dnyi[i]=new.df$dnyi[i]+1
      }
    }
  }
  for(i in 1:(n1+n2-1)){
    if(new.df$group[i]==group1[1]){
      new.df$njia[i+1]<-(new.df$njia[i]-1)
      new.df$nyi[i+1]<-new.df$nyi[i]
      
    }
    if(new.df$group[i]==group1[2]){
      new.df$njia[i+1]<-new.df$njia[i]
      new.df$nyi[i+1]<-(new.df$nyi[i]-1)
      
    }
    if(new.df$deathnumber[i]==0){
      new.df$jiatheory[i]=0
      new.df$yitheory[i]=0
    }
    
    if(new.df$t[i]==new.df$t[i+1]){
      
      new.df$jiatheory[i]=0
      new.df$yitheory[i]=0
      new.df$deathnumber[i+1]<-(new.df$deathnumber[i+1]+new.df$deathnumber[i])
      new.df$dnjia[i+1]<-new.df$dnjia[i+1]+new.df$dnjia[i]
      new.df$dnyi[i+1]<-new.df$dnyi[i+1]+new.df$dnyi[i]
      new.df$total[i+1]<-new.df$total[i]
      new.df$njiaji[i+1]=new.df$njiaji[i]
      new.df$nyiji[i+1]=new.df$nyiji[i]
      new.df$dnjiaji[]
    }
    if(new.df$t[i]!=new.df$t[i+1]){
      new.df$njiaji[i+1]=new.df$njia[i+1]
      new.df$nyiji[i+1]=new.df$nyi[i+1]
      new.df$jiatheory[i]=new.df$deathnumber[i]*new.df$njiaji[i]/new.df$total[i]
      new.df$yitheory[i]=new.df$deathnumber[i]*new.df$nyiji[i]/new.df$total[i]
      if(new.df$total[i]!=1){
        new.df$v[i]=new.df$jiatheory[i]*(new.df$total[i]-new.df$deathnumber[i])*new.df$nyiji[i]/(new.df$total[i]-1)/new.df$total[i]
      }
      new.df$dnjiaji[i]=new.df$dnjia[i]
      new.df$dnyiji[i]=new.df$dnyi[i]
      new.df$indx[i]=1
    }
  }
  new.df$jiatheory[n1+n2]=new.df$deathnumber[n1+n2]*new.df$njiaji[n1+n2]/new.df$total[n1+n2]
  new.df$yitheory[n1+n2]=new.df$deathnumber[n1+n2]*new.df$nyiji[n1+n2]/new.df$total[n1+n2]
  new.df$dnjiaji[n1+n2]=new.df$dnjia[n1+n2]
  new.df$dnyiji[n1+n2]=new.df$dnyi[n1+n2]
  new.df$indx[n1+n2]=1
  return(subset(new.df,new.df$indx!=0)[,c(-1,-4,-5,-12,-13,-16,-17)])
}

#' @title Kaplan-Meier estimate for survival function 
#' @name SurvFuncEsti
#' @description Kaplan-Meier estimate for survival function 
#' @param df the dataframe which consists of three columns:surviving time ,censoring or not(0,1),group
#' @return a data frame consisists of eleven columns,adding Kaplan-Meier estimate for survival function at each time to the results in tableCalculation()
#' @examples
#' \dontrun{
#' SurvFuncEsti(aml1)
#' }
#' @export


SurvFuncEsti<-function(df){
  df1<-tableCalculation(df)
  tt<-nrow(df1)
  S<-rep(0,times=tt)
  S[1]=1-df1$deathnumber[1]/df1$total[1]
  for(i in 2:tt){
    S[i]=S[i-1]*(1-df1$deathnumber[i]/df1$total[1])
  }
  df2=cbind(df1,S)
  names(df2$S)<-"S"
  return(df2)
}

#' @title test if there is a difference between two survival curves.
#' @name survtest
#' @description test if there is a difference between two survival curves using Harrington and Fleming or Peto & Peto method,based on big sample theory and permutation test.
#' @param df the dataframe which consists of three columns:surviving time ,censoring or not(0,1),group.
#' @param B times for permutation.
#' @param method "FH" for Fleming and Harrington,"PP" for Peto & Peto method.
#' @param p a nonnegative parameter used in FM method.
#' @param q like p, controls the weight function S(t)^p(1-S(t))^q.
#' @param modified TRUE for modified Peto & Peto method,FALSE for not.
#' @return a list with components:chisq_statistic for test statistic,p_val_chisq for p value based on big samples,p_val_pur for p value based on permutation test,method for the method you used.
#' @examples
#' \dontrun{
#' survtest(aml1,999,"FH",1,2)
#' survtest(aml1,999,"PP",modified=FALSE)
#' }
#' @export



survtest<-function(df,B=999,method="FH",p=1,q=1,modified=TRUE){
  if(ncol(df)!=3){
    return("error:you must put in a dataframe or matrix which have 3 column")
  }                                             
  
  if(length(unique(df[,3]))!=2){
    return("error:the data should be split into 2 groups")
  }                                             
  
  if(!all(unique(df[,2])%in%c(0,1))){
    return("error:the second line must be of 0,1 to indicate if it is censor data")
  }
  if(method=="FH"){
    new.df<-df
    new<-SurvFuncEsti(df)
    n1=new$njiaji[1]
    n2=new$nyiji[1]
    s<-rep(0,times=nrow(new))
    s[1]=1
    for(i in 2:nrow(new)){
      s[i]=new$S[i-1]
    }
    z1<-sum(s^p*(new$jiatheory-new$dnjiaji)*(1-s)^q)                
    sigama<-sum(s^(2*p)*new$v*(1-s)^(2*q))
    chisq_statistic<-z1^2/sigama
    p_val_chisq<-(1-pchisq(chisq_statistic,1))
    F_H_pur<-rep(0,times=B)
    for(j in 1:B){
      new.df$group<-sample(c(rep(1,times=n1),rep(2,times=n2)),replace=FALSE)
      new<-SurvFuncEsti(new.df)
      for(i in 2:nrow(new)){
        s[i]=new$S[i-1]
      }
      z1pie<-sum(s^p*(new$jiatheory-new$dnjiaji)*(1-s)^q)                
      sigamapie<-sum(s^(2*p)*new$v*(1-s)^(2*q))
      F_H_pur[j]<-z1pie^2/sigamapie
    }
    p_val_pur<-sum(F_H_pur>chisq_statistic)/(1+B)
    hist(F_H_pur,breaks=100,main=paste("purmutation statistic in F_L test with p=",p,"q=",q))
    abline(v =chisq_statistic , col="red", lwd =2)
    legend(5,200,paste("p=",p_val_pur),col="green",bty="n")
    return(list(chisq_statistic=chisq_statistic,p_val_chisq=p_val_chisq,p_val_pur=p_val_pur,method=paste("Fleming_Harriton with p=",p,"q=",q)))
  }
  if(method=="PP"){
    new.df<-df
    new<-SurvFuncEsti(df)
    n1=new$njiaji[1]
    n2=new$nyiji[1]
    if(!modified){
      z1<-sum(new$S*(new$jiatheory-new$dnjiaji))                
      sigama<-sum((new$S^2)*new$v)
      chisq_statistic<-z1^2/sigama
      p_val_chisq<-(1-pchisq(chisq_statistic,1))
      B<-999
      P_P_pur<-rep(0,times=B)
      for(j in 1:B){
        new.df$group<-sample(c(rep(1,times=n1),rep(2,times=n2)),replace=FALSE)
        new<-SurvFuncEsti(new.df)
        z1pie<-sum(new$S*(new$jiatheory-new$dnjiaji))                
        sigamapie<-sum((new$S^2)*new$v)
        P_P_pur[j]<-z1pie^2/sigamapie
      }
      p_val_pur<-sum(P_P_pur>chisq_statistic)/(1+B)
      hist(P_P_pur,breaks=100,main="purmutation statistic in Peto_Peto test,not modified")
      abline(v =chisq_statistic , col="red", lwd =2)
      legend(5,200,paste("p=",p_val_pur),col="green",bty="n")
      return(list(chisq_statistic=chisq_statistic,p_val_chisq=p_val_chisq,p_val_pur=p_val_pur,method="Peto_Peto not modified"))
    }
    if(modified){
      z1<-sum(new$S*new$total*(new$jiatheory-new$dnjiaji)/(new$total+1))                
      sigama<-sum((new$S^2)*new$v*new$total^2/(new$total+1)^2)
      chisq_statistic<-z1^2/sigama
      p_val_chisq<-(1-pchisq(chisq_statistic,1))
      B<-999
      P_P_pur<-rep(0,times=B)
      for(j in 1:B){
        new.df$group<-sample(c(rep(1,times=n1),rep(2,times=n2)),replace=FALSE)
        new<-SurvFuncEsti(new.df)
        z1pie<-sum(new$S*new$total*(new$jiatheory-new$dnjiaji)/(new$total+1))
        sigamapie<-sum((new$S^2)*new$v*new$total^2/(new$total+1)^2)
        P_P_pur[j]<-z1pie^2/sigamapie
      }
      p_val_pur<-sum(P_P_pur>chisq_statistic)/(1+B)
      hist(P_P_pur,breaks=100,main="purmutation statistic in Peto_Peto test,modified")
      abline(v =chisq_statistic , col="red", lwd =2)
      legend(5,200,paste("p=",p_val_pur),col="green",bty="n")
      return(list(chisq_statistic=chisq_statistic,p_val_chisq=p_val_chisq,p_val_pur=p_val_pur,method="Peto_Peto modified"))
    }
  }
  
}

#' @title aml1
#' @name aml1
#' @description Survival in patients with Acute Myelogenous Leukemia
#' @examples 
#' \dontrun{
#' data(aml1)
#' }
NULL


