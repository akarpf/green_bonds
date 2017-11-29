library(data.table)
library(stargazer)
library(lme4)
library(boot)
library(MASS)
library(zoo)
library(reshape)
library(oaxaca)

#contact: Andreas Karpf, andreas.karpf@gmail.com

## set directory and read data
setwd('SET PARENT DIRECTORY')
df  = fread('traddisaggcsv1.csv')

## helper functions
tabular.cast_df <- function(xx,...){
  # a bunch of assumptions that must be met for this function to work:
  if(!require(reshape)) stop("The {reshape} package must be installed for this function to work")
  if(!require(tables)) stop("The {tables} package must be installed for this function to work")
  if(! any(class(xx) == "cast_df")) stop("This function only works for cast_df objects")
  # xx is a casted object
  
  m_xx <- melt(xx)
  rdimnames_xx <- attr(xx, "rdimnames")
  if(length(rdimnames_xx)>2) stop("This function only works for 2D tables")
  
  ROWS <- colnames(rdimnames_xx[[1]])
  COLUMNS <- colnames(rdimnames_xx[[2]])
  colnames_m_xx <- colnames(m_xx)
  
  # This is for cases when one of the equations has "(all)" in them due to something like cast(DATA, x ~.)
  if(all(ROWS == "value")) ROWS <- 1
  if(all(COLUMNS == "value")) COLUMNS <- 1
  
  if(any(colnames_m_xx == "value.1")) {	# then we are supposed to have a "(all)" case (e.g: cast(DATA, .~x)  )
    # m_xx <- m_xx[, -c(which(colnames_m_xx == "value")[-1])] # then remove the column with no value but "(all)"	# This would only work for cast(DATA, x~.) and not for cast(DATA, .~x)
    m_xx[,"value"]   <- m_xx[,"value.1"]
    column_where_all_is <- which(colnames_m_xx  == "value.1")
    m_xx <- m_xx[, -column_where_all_is] # then remove the column with no value but "(all)"
    colnames_m_xx <- colnames(m_xx)
  }
  if(sum(colnames_m_xx == "value") > 1 ) {	# then we are supposed to have a "(all)" case (e.g: cast(DATA, x~.)  )
    # m_xx <- m_xx[, -c(which(colnames_m_xx == "value")[-1])] # then remove the column with no value but "(all)"	# This would only work for cast(DATA, x~.) and not for cast(DATA, .~x)
    column_where_all_is <- which(m_xx[1,] == "(all)")
    m_xx <- m_xx[, -column_where_all_is] # then remove the column with no value but "(all)"
    colnames_m_xx <- colnames(m_xx)
  }
  
  LEFT <- paste(ROWS , collapse="*")
  RIGHT <- paste(COLUMNS , collapse="*")
  
  # turn all ROWS/COLUMNS variables into factors - so to make sure that the tabular will work on them as we expect
  column_to_turn_into_factor <- intersect(c(ROWS, COLUMNS), colnames_m_xx)	# this removes the "1"s in case of cast(DATA, x~.)
  for(i in column_to_turn_into_factor) m_xx[,i] <- factor(m_xx[,i])
  
  v <- function(x) x[1L]
  txt <- paste("tabular(value*v*", LEFT , "~" ,RIGHT ,", data = m_xx, suppressLabels  = 2,...)", sep = "")
  # suppressLabels is in order to remove the value and the v labels (which are added so to make sure the information inside the table is presented)
  eval(parse(text = txt ))
}

summary.oa <- function(res,rn){
  npar = length(res$beta$beta.A)
  napar = names(res$beta$beta.A)
  tfov = round(matrix(res$threefold$overall,ncol = 2, byrow = TRUE),rn)
  tfvar = round(res$threefold$variables,rn)
  tfvar.end = tfvar[,1:2]
  tfvar.coef = tfvar[,3:4]
  tfvar.inte = tfvar[,5:6]
  twf = round(res$twofold$overall,rn)
  twfov1 = matrix(twf[1,2:ncol(twf)],ncol = 2,byrow = TRUE)
  twfov2 = matrix(twf[2,2:ncol(twf)],ncol = 2,byrow = TRUE)
  twfov3 = matrix(twf[3,2:ncol(twf)],ncol = 2,byrow = TRUE)
  twfov4 = matrix(twf[4,2:ncol(twf)],ncol = 2,byrow = TRUE)
  twfov5 = matrix(twf[5,2:ncol(twf)],ncol = 2,byrow = TRUE)
  twfov6 = matrix(twf[5,2:ncol(twf)],ncol = 2,byrow = TRUE)
  twfvar = lapply(res$twofold$variables, function(x) round(x,rn))
  restb = c("Threefold Decomposition", "","")
  restb = rbind(restb,cbind("","Coefficient","Standard Error"))
  restb = rbind(restb,cbind(rbind("Endowment","Coefficients","Interaction"),tfov))
  restb = rbind(restb,c("","",""))
  restb = rbind(restb,c("Detail - Endowment:","",""))
  restb = rbind(restb,cbind(napar,tfvar.end))
  restb = rbind(restb,c("","",""))
  restb = rbind(restb,c("Detail - Coefficients:","",""))
  restb = rbind(restb,cbind(napar,tfvar.coef))
  restb = rbind(restb,c("","",""))
  restb = rbind(restb,c("Detail - Interaction:","",""))
  restb = rbind(restb,cbind(napar,tfvar.inte))
  restb = rbind(restb,c("","",""))
  restb = rbind(restb,c("Twofold Decomposition", "",""))
  restb = rbind(restb,c("","",""))
  restb = rbind(restb,c("Reference Group B", "",""))
  restb = rbind(restb,cbind(c("Explained","Unexplained","Unexplained A","Unexplained B"),twfov1))
  restb = rbind(restb,c("","",""))
  restb = rbind(restb,c("Detail - Explained:","",""))
  restb = rbind(restb,cbind(napar,twfvar[[1]][,2:3]))
  restb = rbind(restb,c("","",""))
  restb = rbind(restb,c("Detail - Unexplained:","",""))
  restb = rbind(restb,c("","",""))
  restb = rbind(restb,cbind(napar,twfvar[[1]][,4:5]))
  restb = rbind(restb,c("","",""))
  restb = rbind(restb,c("Reference Group A", "",""))
  restb = rbind(restb,cbind(c("Explained","Unexplained","Unexplained A","Unexplained B"),twfov2))
  restb = rbind(restb,c("","",""))
  restb = rbind(restb,c("Detail - Explained:","",""))
  restb = rbind(restb,c("","",""))
  restb = rbind(restb,cbind(napar,twfvar[[2]][,2:3]))
  restb = rbind(restb,c("","",""))
  restb = rbind(restb,c("Detail - Unexplained:","",""))
  restb = rbind(restb,cbind(napar,twfvar[[2]][,4:5]))
  restb = rbind(restb,c("","",""))
  restb = rbind(restb,c("Equally weighted (Reimers 1983)", "",""))
  restb = rbind(restb,c("","",""))
  restb = rbind(restb,cbind(c("Explained","Unexplained","Unexplained A","Unexplained B"),twfov3))
  restb = rbind(restb,c("","",""))
  restb = rbind(restb,c("Detail - Explained:","",""))
  restb = rbind(restb,c("","",""))
  restb = rbind(restb,cbind(napar,twfvar[[3]][,2:3]))
  restb = rbind(restb,c("","",""))
  restb = rbind(restb,c("Detail - Unexplained:","",""))
  restb = rbind(restb,c("","",""))
  restb = rbind(restb,cbind(napar,twfvar[[3]][,4:5]))
  restb = rbind(restb,c("","",""))
  restb = rbind(restb,c("Weighted Average (Cotton 1988)", "",""))
  restb = rbind(restb,c("","",""))
  restb = rbind(restb,cbind(c("Explained","Unexplained","Unexplained A","Unexplained B"),twfov4))
  restb = rbind(restb,c("","",""))
  restb = rbind(restb,c("Detail - Explained:","",""))
  restb = rbind(restb,c("","",""))
  restb = rbind(restb,cbind(napar,twfvar[[4]][,2:3]))
  restb = rbind(restb,c("","",""))
  restb = rbind(restb,c("Detail - Unexplained:","",""))
  restb = rbind(restb,c("","",""))
  restb = rbind(restb,cbind(napar,twfvar[[4]][,4:5]))
  restb = rbind(restb,c("","",""))
  restb = rbind(restb,c("Pooled Regression w/o group (Neumark 1988)", "",""))
  restb = rbind(restb,c("","",""))
  restb = rbind(restb,cbind(c("Explained","Unexplained","Unexplained A","Unexplained B"),twfov5))
  restb = rbind(restb,c("","",""))
  restb = rbind(restb,c("Detail - Explained:","",""))
  restb = rbind(restb,c("","",""))
  restb = rbind(restb,cbind(napar,twfvar[[5]][,2:3]))
  restb = rbind(restb,c("","",""))
  restb = rbind(restb,c("Detail - Unexplained:","",""))
  restb = rbind(restb,c("","",""))
  restb = rbind(restb,cbind(napar,twfvar[[5]][,4:5]))
  restb = rbind(restb,c("","",""))
  restb = rbind(restb,c("Pooled Regression w/ group (Jann 2008)", "",""))
  restb = rbind(restb,c("","",""))
  restb = rbind(restb,cbind(c("Explained","Unexplained","Unexplained A","Unexplained B"),twfov6))
  restb = rbind(restb,c("","",""))
  restb = rbind(restb,c("Detail - Explained:","",""))
  restb = rbind(restb,c("","",""))
  restb = rbind(restb,cbind(napar,twfvar[[6]][,2:3]))
  restb = rbind(restb,c("","",""))
  restb = rbind(restb,c("Detail - Unexplained:","",""))
  restb = rbind(restb,c("","",""))
  restb = rbind(restb,cbind(napar,twfvar[[6]][,4:5]))
  restb = rbind(restb,c("","",""))
  tval = as.numeric(restb[-c(1,2),2])/as.numeric(restb[-c(1,2),3])
  restb = cbind(restb,c("","t-statistic",round(tval,rn)))
  pval = suppressMessages(2 * pt(abs(as.numeric(restb[,4])), res$reg$reg.pooled.1$df.residual, lower.tail = FALSE))
  restb = cbind(restb,c("","p-value",round(pval[-c(1,2)],rn)))
  restb[is.na(restb)] = ""
  colnames(restb) = restb[2,]
  restb[2,] = c("","","","","")
  rownames(restb) = rep("",nrow(restb))
  restb
}

## oaxaca blinder (Reference Group and Neumark (1988)) with mixed effects - extension of the code in the R library oaxaca
mixedoaxaca = function(data,ind,formula,inde,grpv,rfgrp){
  
  ####
  dfr = data[ind]
  grp = as.vector(dfr[,grpv,with = F] == rfgrp)
  sa  = dfr[grp]
  sb = dfr[!grp]
  cids = intersect(unlist(unique(sa[,inde,with = F])),unlist(unique(sb[,inde,with = F])))
  dfr = dfr[unlist(dfr[,inde,with=F]) %in% cids]
  sa = sa[unlist(sa[,inde,with=F]) %in% cids]
  sb = sb[unlist(sb[,inde,with=F]) %in% cids]
  
  ###
  a <- suppressMessages(lmer(formula,data=sa))
  ra <- model.matrix(a,"random")
  fa <- model.matrix(a,"fixed")
  ca <- as.matrix(coefficients(a)[[1]])
  prog <<- prog + 1
  setTxtProgressBar(pb, prog)
  
  b <- suppressMessages(lmer(formula,data=sb))
  rb <- model.matrix(b,"random")
  fb <- model.matrix(b,"fixed")
  cb <- as.matrix(coefficients(b)[[1]])
  prog <<- prog + 1
  setTxtProgressBar(pb, prog)
  
  f <- suppressMessages(lmer(formula,data=dfr))
  rf <- model.matrix(f,"random")
  ff <- model.matrix(f,"fixed")
  cf <- as.matrix(coefficients(f)[[1]])
  prog <<- prog + 1
  setTxtProgressBar(pb, prog)
  
  #overall decomposition 
  S <- function(b,z,c){
    mean(rowSums(b * (z %*% c)))
  }
  
  #weight function endowment/explained
  wx <- function(b,x1,x2){
    N <- colMeans(t((colMeans(x1) - colMeans(x2)) * t(b)))
    D <- sum(N)
    return(N/D)
  }
  
  #weight function coefficients/unexplained
  wb <- function(x,b1,b2){
    N  <- colMeans(t(colMeans(x) * t(b1-b2)))
    D  <- sum(N)
    return(N/D)
  }
  
  #overall decomposition: reference group
  drg_f  <- (S(fa,ra,ca) - S(fb,rb,ca)) + (S(fb,rb,ca)-S(fb,rb,cb))
  drg_e  <- (S(fa,ra,ca) - S(fb,rb,ca)) #explained
  drg_u  <- (S(fb,rb,ca)-  S(fb,rb,cb)) #unexplained
  
  #overall decomposition: Neumark
  dnm_f  <- (S(fa,ra,ca) - S(fa,ra,cf)) + (S(fb,rb,cf)-S(fb,rb,cb)) + (S(fa,ra,cf)-S(fb,rb,cf))
  dnm_e  <- (S(fa,ra,cf) - S(fb,rb,cf)) #explained                 
  dnm_u  <- (S(fa,ra,ca) - S(fa,ra,cf)) + (S(fb,rb,cf)-S(fb,rb,cb)) #unexplained  
  
  #detailed decomposition
  wtsx_drg <- wx(ca,fa,fb) #endowment: reference group 
  wtsb_drg <- wb(fa,ca,cb) #coefficients: reference group
  drg_e_d <- drg_e*wtsx_drg
  drg_u_d <- drg_u*wtsb_drg
  
  wtsx_dnm <- wx(cf,ff,fb) #endowment: Neumark
  wtsb_dnm <- wb(ff,cf,cb) #coefficients: Neumark
  dnm_e_d <- dnm_e*wtsx_dnm
  dnm_u_d <- dnm_u*wtsb_dnm
  
  res <- c(drg_f,drg_e,drg_u,drg_e_d,drg_u_d,dnm_f,dnm_e,dnm_u,dnm_e_d,dnm_u_d)
  names(res) <- c("rg_diff","rg_exp","rg_unex",paste0("rg_ex_",colnames(cf)),paste0("rg_unex_",colnames(cf)),
                  "nm_diff","nm_exp","nm_unex",paste0("nm_ex_",colnames(cf)),paste0("nm_unex_",colnames(cf)))
  return(res)
}

## above + bootstrapped standard errors
mixedoaxacabs = function(data,formula,inde,grpv,rfgrp,bsiter){
  
  grp = as.vector(data[,grpv,with = F] == rfgrp)
  na  = sum(grp)
  nb = sum(!grp)
  prog <<- 0
  print("Running functions over bootstrap samples")
  pb <<- txtProgressBar(min = 0, max = bsiter * 3, style = 3)
  oab <- boot(data,mixedoaxaca,R = bsiter,formula = formula,inde = inde,grpv = grpv,rfgrp = rfgrp)
  
  res <- cbind(oab$t0,apply(oab$t,2,function(x) sqrt(var(x))))
  
  ttest <- function(x){
    res <- t.test(x, mu = 0, conf.level = 0.95)
    return(rbind(res$statistic,res$conf.int[1],res$conf.int[2],res$p.value))
  }
  
  res <- cbind(res,t(apply(oab$t,2,ttest)))
  res <- as.matrix(res)
  restb <- c("Reference Group","","","","","","","")
  restb <- rbind(restb,c("","Overall",rep("",6)))
  restb <- rbind(restb,cbind(c("","",""),c("Difference","Explained","Unexplained"),res[rownames(res) %in% c("rg_diff","rg_exp","rg_unex"),]))
  restb <- rbind(restb,rep("",8))
  restb <- rbind(restb,c("","Explained - Detail",rep("",6)))
  rgexp <- res[grep("rg_ex_",rownames(res)),]
  restb <- rbind(restb,cbind(rep("",nrow(rgexp)),sub("rg_ex_","",rownames(rgexp)),rgexp))
  restb <- rbind(restb,rep("",8))
  restb <- rbind(restb,c("","Unexplained - Detail",rep("",6))) 
  rgune <- res[grep("rg_unex_",rownames(res)),]
  restb <- rbind(restb,cbind(rep("",nrow(rgune)),sub("rg_unex_","",rownames(rgune)),rgune))
  restb <- rbind(restb,rep("",8))
  restb <- rbind(restb,c("Neumark (1988)","","","","","","",""))
  restb <- rbind(restb,rep("",8))
  restb <- rbind(restb,c("","Overall",rep("",6))) 
  restb <- rbind(restb,cbind(c("","",""),c("Difference","Explained","Unexplained"),res[rownames(res) %in% c("nm_diff","nm_exp","nm_unex"),]))
  restb <- rbind(restb,rep("",8))
  restb <- rbind(restb,c("","Explained - Detail",rep("",6))) 
  nmexp <- res[grep("nm_ex_",rownames(res)),]
  restb <- rbind(restb,cbind(rep("",nrow(nmexp)),sub("nm_ex_","",rownames(nmexp)),nmexp))
  restb <- rbind(restb,rep("",8))
  restb <- rbind(restb,c("","Unexplained - Detail",rep("",6))) 
  nmune <- res[grep("nm_unex_",rownames(res)),]
  restb <- rbind(restb,cbind(rep("",nrow(nmune)),sub("nm_unex_","",rownames(nmune)),nmune))
  restb <- rbind(restb,rep("",8))
  restb <- rbind(restb,c("Sample Size","","","","","","",""))
  restb <- rbind(restb,c("","Group A",na,rep("",5)))
  restb <- rbind(restb,c("","Group B",nb,rep("",5)))
  restb <- rbind(restb,c("Bootstrap Repl.",bsiter,rep("",6)))
  rownames(restb) <- NULL
  restb <- data.frame(restb)
  nacols = (restb[,3] == 0) & !is.na(restb[,3])
  replmat <- restb[nacols,-(1:2)]
  restb <- data.frame(restb[,1:2],sapply(restb[,3:8],function(x) round(as.numeric(as.character(x)),4)))
  restb[nacols,-(1:2)] = matrix('-',ncol = ncol(replmat),nrow = nrow(replmat))
  colnames(restb) <- c(" ", " ","Est.","Std.Err.","t",".25 pct.",".975 pct.","P >|t|")
  out <- list()
  out$restb <- restb
  out$obd <- oab
  out$ra <- suppressMessages(lmer(formula,data=df[grp]))
  out$rb <- suppressMessages(lmer(formula,data=df[!grp]))
  out$rp <- suppressMessages(lmer(formula,data=df))
  return(out)
}

df[,tradedate:=as.IDate(tradedate)]
df[,frequency:={
  SD1 = .SD[,list(tradedate,NA)]
  SD1[,freq:=nrow(SD1[.SD$tradedate >= tradedate & tradedate >= (.SD$tradedate - 30)]),by=1:nrow(SD1)]
  SD1$freq
  },by =idcs]
df$grp = (df$type == "brown bond")
df$iss = as.factor(as.numeric(as.factor(df$issuerid)))
df = df[df$ytm > 0,]
df = df[df$ytm < quantile(df$ytm,0.999),]
df[,ratingA:= as.factor(ratingc == "A")]
df[,ratingB:= as.factor(ratingc == "B")]
df[,ratingN:= as.factor(ratingc == "NR")]
df[state == "AP",state:="MD"]
df[,state:=as.factor(state)]
df[,year := as.numeric(year(as.IDate(tradedate)))]

## add infos about US federal states
statedebt = fread('statedebt.csv')
stateabbr = fread('stateabbrev.csv',header = F,col.names = c("state","abbr"))
statedebt = merge(statedebt,stateabbr,by="state",all.x = TRUE)
setnames(statedebt,c("state","abbr"),c("statename","state"))
df = merge(df,statedebt,by = c("year","state"),all.x = T)
df[,state_local_debt:=100 * statelocaldebt/gsp]
df[,state_debt:=100 * statedebt/gsp]
df[,local_debt:=100 *localdebt/gsp]
df[,volume:=as.numeric(amount/paissuance)]
df = df[volume<1]
df[,tradetype:=as.factor(tradetype)]
df[,broker:=(tradetype == "D")]
df[,otc:=(tradetype == "S" | tradetype == "P")]
#save(df,file="dfdatapre.RData")

sumt = melt.data.table(df[,list(type,ytc,dtm,amount,paissuance,local_debt,state_debt,realgrowth)], id=c("type"), na.rm=TRUE)
sumt1 = copy(sumt)
sumt1[,type:='overall']
sumt = rbind(sumt,sumt1)
setnames(sumt,c("type","variable2","value"))
sumtablep = tabular.cast_df(cast(sumt,type ~ variable2, c(mean,sd)))
sumtablep1 = as.matrix(sumtablep)[-c(1,3),]

## table - descriptive statistics
stargazer(sumtablep1,type = "latex",align = T,column.sep.width="-10pt",summary = F,header = F,title="Descriptive statistics of the most relevant continuous variables used in the Oxaca-Blinder Decomposition for
          the subsets of green and brown bonds respectively as well as the overall data set",
          out = 'summarystat.tex',font.size = 'scriptsize', label = 'tab:summarystat')

df[,population:=scale(population)]
df[,dtm:=scale(dtm)]
df[,paissuance:=scale(paissuance)]
df[,freq:=scale(frequency)]
df[,amount:=scale(amount)]


specifi = "ytc ~ dtm + amount + paissuance + freq + ratingA + ratingB + treasury + broker + local_debt + state_debt + realgrowth + population + (1 | iss)"
#save(specifi, df, file = "datafortimedoaxaca.RData")
#load("datafortimedoaxaca.RData")

oab = mixedoaxacabs(df,specifi,"iss","type","brown bond",100)

## produce oxaca blinder output tables
# regressions
stargazer(oab$restb,type = "latex",summary = F,rownames = F,title="Results of the Oaxaca-Blinder Decomposition - Reference Group and Neumark (1988) (dtm: days to maturity; amount: 
volume of the transaction; A and B Rating: dummy variable indicating whether the bond has a rating of A or B; local\\_debt: the accumulated local debt in a state,
state\\_debt: state debt; realgrowth: real gsp growth of the state, population: population of a state)",
          out = 'oabresult.tex',font.size = 'scriptsize', label = 'tab:oabresult')

# decomposition
stargazer(oab$rp,oab$ra,oab$rb,type = "latex",summary = T,rownames = F,title="Regression Results: Conventional Bonds (1), Green Bonds (2), Pooled Data (3)",
          out = 'regres.tex',font.size = 'scriptsize', label = 'tab:regres')

#stargazer(summary.oa(oab,3),type = "text",summary = F,rownames = F)
