library(robustDist)


regular.ls.loss = function(log.branch.length, my.topology, seq.dist){
  
  ## first bring branch length to the absolute scale
  branch.length = exp(log.branch.length)
  
  num.taxa = dim(seq.dist)[1]
  
  my.ls = 0
  
  my.phylo = my.topology
  my.phylo$edge.length = branch.length
  
  phylo.dist = cophenetic.phylo(my.phylo)
  
  for (i in 2:num.taxa){
    for (j in 1:(i-1)){
      my.ls = my.ls + (phylo.dist[my.phylo$tip.label[i],my.phylo$tip.label[j]] - seq.dist[my.phylo$tip.label[i],my.phylo$tip.label[j]])^2
    }
  }
  
  return(my.ls)
}


regular.ls.fit = function(init.brlen, my.topology, seq.dist){
  
  return.val = NULL
  par.num = length(init.brlen)
  
  optim.out = optim(
    log(init.brlen),
    regular.ls.loss,
    lower = rep(-10,par.num),	# this should constrain br len to be min of 4.5e-5
    upper = rep(0,par.num),	# this should constrain br len to be max of 1
    my.topology = my.topology,
    seq.dist = seq.dist,
    control = list(fnscale=+1),
    method = "L-BFGS-B"
  )
  
  if (optim.out$convergence != 0){
    stop("loss function minimization routine did not converge")
  }else{
    return.val = exp(optim.out$par)
  }
  
  return(return.val)
}




new.ls.loss = function(log.branch.length, my.topology, seq.table, regist.matrix){
  
  ## first bring branch length to the absolute scale
  branch.length = exp(log.branch.length)
  
  num.taxa = dim(seq.table)[1]
  
  my.ls = 0
  
  my.phylo = my.topology
  my.phylo$edge.length = branch.length
  
  phylo.dist = cophenetic.phylo(my.phylo)

  ## define JC69 model with its eigen decomposition
  regist.mat = matrix(1, nrow = 4, ncol = 4) - diag(1, 4)
  
  for (i in 2:num.taxa){
    for (j in 1:(i-1)){
      my.jc <- as.eigen(jc.mc(4/3,c(0.25,0.25,0.25,0.25)))
      my.jc = rescale.mc(my.jc, phylo.dist[my.phylo$tip.label[i],my.phylo$tip.label[j]])
      my.ls = my.ls + (phylo.dist[my.phylo$tip.label[i],my.phylo$tip.label[j]] - pair.robust.dist(my.jc, regist.matrix, seq.table[c(my.phylo$tip.label[i],my.phylo$tip.label[j]),]))^2
    }
  }
  
  return(my.ls)
}



new.ls.fit.optimx <- function(init.brlen, my.topology, seq.table, method="nlminb", low=-100){
  
  return.val = NULL
  par.num = length(init.brlen)
  regist.mat = matrix(1,4,4) - diag(rep(1,4))
  optim.out<-list(par=1,convcode=1)
  count<-0
  
  # try up to 10 times to reach convergence
  while((optim.out$convcode != 0) & count<10){
    
    optim.out <- optimx(
      log(init.brlen),
      new.ls.loss,
      lower = rep(low,par.num),	# this should constrain br len to be min of 4.5e-5 (-100 gives 3.72e-44)
      upper = rep(0,par.num),		# this should constrain br len to be max of 1
      method = method,
      my.topology = my.topology,
      seq.table = seq.table,
      regist.matrix = regist.mat
    )
    count<-count+1
    init.brlen<-runif(par.num,0,0.5)
  }
  
#  if(!any(unlist(optim.out$par)==0) & any(unlist(optim.out$par)==low)){
#    optim.out$conv<-6		# so, 6 means that it hit the lower boundary
#  }
#  if(any(unlist(optim.out$par)==0) & !any(unlist(optim.out$par)==low)){
#    optim.out$conv<-7		# so, 7 means that it hit the upper boundary
#  }
#  if(any(unlist(optim.out$par)==0) & any(unlist(optim.out$par)==low)){
#    optim.out$conv<-8		# so, 8 means that it hit both boundaries
#  }
  
  
  return.val<-list(par.est=exp(unlist(optim.out[1:par.num])),ls=unlist(optim.out$value),conv=unlist(optim.out$convcode),count=count)
  return(return.val)
}



read.phylosim.nuc<-function(alignment){
  n.taxa<-nrow(alignment)
  n.sites<-ncol(alignment)
  matrixed.align<-t(matrix(as.vector(unlist(apply(alignment,1,strsplit,""))),ncol=n.taxa))
  
  not.in.alphabet <- ((matrixed.align != "A")*(matrixed.align != "G")*(matrixed.align != "C")*(matrixed.align != "T")*
                        (matrixed.align != "a")*(matrixed.align != "g")*(matrixed.align != "c")*(matrixed.align != "t"))
  not.in.alphabet <- not.in.alphabet == 1
  
  align.table = matrix(0, nrow = n.taxa, ncol = n.sites)
  rownames(align.table) = rownames(alignment)
  ##print(not.in.alphabet)
  
  align.table[not.in.alphabet] = 0
  align.table[matrixed.align == "A" | matrixed.align == "a"] = 1
  align.table[matrixed.align == "G" | matrixed.align == "g"] = 2
  align.table[matrixed.align == "C" | matrixed.align == "c"] = 3
  align.table[matrixed.align == "T" | matrixed.align == "t"] = 4
  return(align.table)
}



