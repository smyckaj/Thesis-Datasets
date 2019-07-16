formula_builder=function(set){
  #######################
  #a small function for parsing set elements into formula system
  #######################
  if (length(set)>1) {paste(set[2:length(set)], set[1],sep="~")}
}


zero_builder=function(set){
  #######################
  #a small function for parsing set elements into 0
  #######################
  if (length(set)>=1) {paste(set,"0",sep="~")}
}


restrict_classe2geosse=function(table, forbidden=NULL,
                                single.sympatry=F, no.sympatry=F, 
                                single.vicariance=F, no.vicariance=F, 
                                pooled.founder=F, single.founder=F,no.founder=F,
                                single.extinction=F, no.extinction=F,
                                pooled.migration=F, single.migration=F, no.migration=F){

  #######################
  #input - states table
  #output - list of formulas to be used as constrain(lik,formulae=formulas) in diversitree
  #######################  

  #libraries
  library(diversitree)
  library(rje)

  #table into set representation
  tableset=list()
  for (i in 1:length(table[,1])) {tableset[[i]]=names(table)[as.logical(table[i,])]}
  #get parameter names
  phy <- rcoal(100)
  names=names(starting.point.classe(phy, k=length(table[,1])))
  names=setdiff(names,forbidden)
  
  #parameter states from strings
  #############################
  #branching for k>10, to be done for k>100
  if (length(table[,1])>10) {zp=1} else {zp=0}
  
  #lambdas
  lambdas=names[grep("lambda",names)]

  fromlambdas=as.numeric(substr(lambdas, 7, 7+zp*1))
  tolambdas1=as.numeric(substr(lambdas, 8+zp*1, 8+zp*2))
  tolambdas2=as.numeric(substr(lambdas, 9+zp*2, 9+zp*3))

  #mus
  mus=names[grep("mu",names)]
  frommus=as.numeric(substr(mus, 3, 3+zp*1))

  #qs
  qs=names[grep("q",names)]
  fromqs=as.numeric(substr(qs, 2, 2+zp*1))
  toqs=as.numeric(substr(qs, 3+zp*1, 3+zp*2))
  
  #biogeographic processes definition
  ###########################
  #sympatric speciation
  sympatry_logical_list=list()
  sympatry_list=list()
  
  for (j in 1:length(names(table))){
    
    sympatry_logical_list[[j]]=rep(F,length(lambdas))
    
    for (i in 1:length(lambdas)){
      sympatry_logical_list[[j]][i]=(setequal(tableset[[fromlambdas[i]]],tableset[[tolambdas1[i]]]) && is.subset(tableset[[tolambdas2[i]]],tableset[[fromlambdas[i]]]) && setequal(tableset[[tolambdas2[i]]], names(table)[j])) | #2 is the offspring species
                          (setequal(tableset[[fromlambdas[i]]],tableset[[tolambdas2[i]]]) && is.subset(tableset[[tolambdas1[i]]],tableset[[fromlambdas[i]]]) && setequal(tableset[[tolambdas1[i]]], names(table)[j]))   #1 is the offspring species
    }
    
    sympatry_list[[j]]=lambdas[sympatry_logical_list[[j]]]
    
  }  
  
  #vicariance speciation
  vicariance_logical=rep(F,length(lambdas))
  for (i in 1:length(lambdas)){
    vicariance_logical[i]=(length(intersect(tableset[[tolambdas1[i]]],tableset[[tolambdas2[i]]]))==0 && #it is allopatry
                         setequal(tableset[[fromlambdas[i]]],union(tableset[[tolambdas1[i]]],tableset[[tolambdas2[i]]]))) #union of offsprings is equal to ancestor
  }

  vicariance=lambdas[vicariance_logical]

  #founder event speciations (sensu BioGeoBears)
  founder_logical_list=list()
  founder_list=list()  
  
  for (j in 1:length(names(table))){
    
    founder_logical_list[[j]]=rep(F,length(lambdas))
    
    for (i in 1:length(lambdas)){
      founder_logical_list[[j]][i]=(length(intersect(tableset[[tolambdas1[i]]],tableset[[tolambdas2[i]]]))==0 && #it is allopatry
                                   ((setequal(tableset[[fromlambdas[i]]],tableset[[tolambdas1[i]]]) && setequal(tableset[[tolambdas2[i]]], names(table)[j])) | #2 is offspring species
                                      (setequal(tableset[[fromlambdas[i]]],tableset[[tolambdas2[i]]]) && setequal(tableset[[tolambdas1[i]]], names(table)[j])))) #1 is offspring species
    }
    
    founder_list[[j]]=lambdas[founder_logical_list[[j]]]
    
  }
  
  #extinction
  extinction_q_logical_list=list()
  extinction_mu_logical_list=list()
  extinction_list=list()
  
  for (j in 1:length(names(table))){
    
    #local extinction
    extinction_q_logical_list[[j]]=rep(F,length(qs))
    
    for (i in 1:length(qs)){
      extinction_q_logical_list[[j]][i]=is.subset(tableset[[toqs[i]]],tableset[[fromqs[i]]]) && #offspring is subset of ancestor
                                        setequal(setdiff(tableset[[fromqs[i]]],tableset[[toqs[i]]]),names(table)[j]) #their difference is focal area
    }
    
    #global extinction
    extinction_mu_logical_list[[j]]=rep(F,length(mus))
    
    for (i in 1:length(mus)){
      extinction_mu_logical_list[[j]][i]=setequal(tableset[[frommus[i]]],names(table)[j]) #it is extinction in focal area
    }
    
    extinction_list[[j]]=c(qs[extinction_q_logical_list[[j]]],mus[extinction_mu_logical_list[[j]]])
    
  }
  
  

  #migration
  migration_logical_list=list()
  migration_list=list()
  
  for (j in 1:length(names(table))){
    
    migration_logical_list[[j]]=rep(F,length(qs))
    
    for (i in 1:length(qs)){
      migration_logical_list[[j]][i]=is.subset(tableset[[fromqs[i]]],tableset[[toqs[i]]]) && #ancestor is a subset of offspring
                                     setequal(setdiff(tableset[[toqs[i]]],tableset[[fromqs[i]]]), names(table)[j]) #their difference is focal area
    }
    
    migration_list[[j]]=qs[migration_logical_list[[j]]]
    
  }
  
  #everything else and the forbidden parameters is a zero combination
  zerocombinations=c(setdiff(names, c(unlist(sympatry_list), vicariance, unlist(founder_list), unlist(extinction_list), unlist(migration_list))),forbidden)


  #build formulas 
  ##############################
  
  #sympatry
  if (no.sympatry) {
    sympatry_formulas=zero_builder(unlist(sympatry_list))
  } else if (single.sympatry) {
    sympatry_formulas=formula_builder(unlist(sympatry_list))
  } else {
    sympatry_formulas=unlist(lapply(sympatry_list, formula_builder)) 
  }
  
  #vicariance
  if (no.vicariance) {
    vicariance_formulas=zero_builder(vicariance)
  } else if (single.vicariance) {
    vicariance_formulas=formula_builder(vicariance)
  } else {
    vicariance_formulas=NULL
  }
  
  #founder
  if (no.founder) {
    founder_formulas=zero_builder(unlist(founder_list))
  } else if (single.founder) {
    founder_formulas=formula_builder(unlist(founder_list))
  } else if (pooled.founder){
    founder_formulas=unlist(lapply(founder_list, formula_builder))
  } else {
    founder_formulas=NULL
  }
  
  #extinction
  if (no.extinction) {
    extinction_formulas=zero_builder(unlist(extinction_list))
  } else if (single.extinction) {
    extinction_formulas=formula_builder(unlist(extinction_list))
  } else {
    extinction_formulas=unlist(lapply(extinction_list, formula_builder))
  }
    
  #migration
  if (no.migration) {
    migration_formulas=zero_builder(unlist(migration_list))
  } else if (single.migration) {
    migration_formulas=formula_builder(unlist(migration_list))
  } else if (pooled.migration) {
    migration_formulas=unlist(lapply(migration_list, formula_builder))
  } else {
    migration_formulas=NULL
  }
  
  #not meaningful parameters
  zerocombinations_formulas=zero_builder(zerocombinations)
  
  #stick them in one list
  ##############################
  formulas=as.list(c( sympatry_formulas,
                      vicariance_formulas,
                      founder_formulas,
                      extinction_formulas,
                      migration_formulas,
                      zerocombinations_formulas))
  return(formulas)
}

##usage
#regionstable=data.frame(A=c(1,0,0,1,1,0,1),B=c(0,1,0,1,0,1,1),C=c(0,0,1,0,1,1,1))
#formulas=restrict_classe2geosse(regionstable, no.founder = T, single.extinction = T, no.vicariance=T) 


