# Trait by environment (example with yield)
ENV=function(trait="yield"){
  data.line=data.check=gen.check=gen.line=NULL
  data(soynam,envir=environment(),package="SoyNAM")
  test=dcast(data.line,strain~environ,value.var=trait,mean)
  Strain=as.character(test[,1])
  test=test[,-1]
  test=data.matrix(test)
  rownames(test)=Strain
  test=test[rownames(gen.line),]
  test[is.nan(test)]=NA
  return(test)}

# Main function
BLUP=function(trait="yield",family="all",env="all",use.check=TRUE,clean.rep=TRUE){
    
  # arguments
  # trait = c("yield","maturity","height","lodging","protein","oil","size","fiber")
  # family = numbers or "all"
  # env = 1:18 or "all"
  # use.check = TRUE or FALSE 
  # clean.rep = TRUE or FALSE (remove repeated lines)
  # model: trait = check(fixed) + env(random) + line(random)
  
  # Load data
  data.line=data.check=gen.check=gen.line=NULL
  data(soynam,envir=environment(),package="SoyNAM")
  
  # FAM
  fam=rownames(gen.line)
  fam=gsub('DS1.-','',fam)
  fam=gsub('...$','',fam,perl = T)
  fam=as.numeric(fam)
  
  # CHR
  chr=rep(NA,20)
  for(i in 1:20) chr[i]=length(grep(
      paste("Gm",sprintf("%02d",i),sep=''),
      colnames(gen.line)));rm(i)
    
  if(is.numeric(family)) data.line = data.line[data.line$family%in%family,]

  
  if(is.numeric(env))  if(length(env)==1) stop("At least two environments where the trait was measured are required")
  
  if(is.numeric(env)){
    E1 = as.numeric(data.line$environ)
    data.line = data.line[E1%in%env,]
    E2 = as.numeric(data.check$environ)
    data.check = data.check[E2%in%env,]
  }
    
  # Check function
  CHECK=function(trait){ test=dcast(data.check,environ+spot~strain,value.var=trait,mean)
                         rownames(test)=test[,2];E=test[,1];test=test[,-c(1,2)];test=data.matrix(test);test[is.nan(test)]=NA;
                         X=function(X) unlist(tapply(X,E,FUN=function(x){m=mean(x,na.rm=T);SD=sd(x,na.rm=T);return((x-m)/SD)}))
                         MEAN=apply(test,2,X);C=rowMeans(MEAN,na.rm=T);names(C)=rownames(test);C[is.nan(C)]=0;return(C)}
    
  # Model terms
  Y = data.line[,trait]
  G = data.line[,"strain"]
  E = data.line[,"environ"]
  
  # BLUP
  if(use.check){
    cat('solving BLUE of checks\n')
    check = CHECK(trait);set = as.character(data.line[,"spot"])
    C = check[set]
    cat('solving BLUP of phenotypes\n')
    blup=lmer(Y~C+(1|E)+(1|G))
  }else{
    cat('solving BLUP of phenotypes\n')
    blup=lmer(Y~(1|E)+(1|G))}
  BV = coef(blup)$G[,1]
  names(BV) = rownames(coef(blup)$G)
  BV = BV[rownames(gen.line)]
  
  if(is.numeric(family)) {
    BV = BV[fam%in%family]
    gen.line = gen.line[fam%in%family,]
    fam = fam[fam%in%family]
  }
    
  if(clean.rep){
    cat('removing repeated genotypes\n')
    d=cleanREP(y = cbind(BV,BV),fam = fam,gen = gen.line)
    LIST = list(Phen=d$y[,1],Gen=d$gen,Chrom=chr,Fam=d$fam)
  }else{
    LIST = list(Phen=BV,Gen=gen.line,Chrom=chr,Fam=fam)
  }
  
  return(LIST)
}

