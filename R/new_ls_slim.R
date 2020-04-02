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


# order of arguments changed 4/1/20, might have implications on old code...
new.ls.fit.optimx <- function(my.topology, seq.table, init.brlen = NULL, method="nlminb", low=-100, high=2){

  if(is.null(init.brlen)){
    init.brlen <- rep(0.1, dim(my.topology$edge)[1])
  }

  # change A,C,T,G to numeric if needed
  if(typeof(seq.table) == "character"){
    seq.table <- read.phylosim.nuc(seq.table)
  }

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
      upper = rep(high,par.num),		# this should constrain br len to be max of e^2 ~ 7.389
      method = method,
      my.topology = my.topology,
      seq.table = seq.table,
      regist.matrix = regist.mat
    )

    # Keep track of iterations
    count<-count+1

    # If convergence was not reached, re-initialize starting values to random values
    init.brlen<-runif(par.num,0,0.5)
  }


# 3/24/20 keeping this out for now as it doesn't seem necessary...

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


#' Least squares phylogenetic inference
#'
#' This function performs phylogenetic inference via ordinary least squares. It was written taking elements
#' from Liam Revell's optim.phylo.ls. The main difference is that it allows for an exhaustive search among
#' all possible topologies (if not, it will do an NNI search, starting from the NJ tree). This function infers an unrooted tree.
#'
#' @param seq.table a nucleotide sequence alignment, formatted as an n x s character matrix where n = # of taxa, s = # of sites
#' @param set.neg.to.zero if TRUE, negative branch lengths will be converted to 0
#' @param search.all if TRUE, an exhaustive search across all topologies will be performed. Otherwise, an NNI search will be performed.
#' @param model substitution model for which to calculate the distance matrix
#' @param tol in NNI search, keep searching if improvement is at least this amount
#'
#' @keywords phylogeny, OLS
#' @export
#' @examples
#' phylo.ls(seq.table)
phylo.ls <- function(alignment, set.neg.to.zero = TRUE, search.all = FALSE, model="JC69", tol = 1e-10){
  seq.table <- as.character(alignment)
  data.bin<-as.DNAbin(as.alignment(seq.table))
  D <- as.matrix(dist.dna(data.bin,model=model))
  n <- nrow(D)

#    if(class(D) == "dist"){
#      D <- as.matrix(D)
#    }
#    n <- nrow(D)

  if(search.all){
    # Do search through all possible topologies
    all.trees <- allTrees(n, tip.label = row.names(as.matrix(D))) # change this (as.matrix(D) is redundant)
    allQ <- vector()
    bestQ <- Inf

    # Benchmarked [for loop] speed against [lapply] up to n.tips=8, basically the same either way (as I should've known).
    for (i in 1:length(all.trees)) {
      all.trees[[i]]$edge.length <- rep(dim(all.trees[[i]]$edge)[1], 1)
      all.trees[[i]] <- ls.tree(tree=all.trees[[i]], D=D)	# this used to be ls.tree. 3/30/20 back to ls.tree because nnls.tree produces ties
      allQ[i] <- attr(all.trees[[i]], "Q-score")
    }
    best <- which(allQ == min(allQ))
    best.tree <- all.trees[[best]]

  } else {
    # Do nni search, starting from NJ tree
    tree <- nj(D)
    best.tree <- ls.tree(tree = tree, D = D)

    Q <- Inf
    bestQ <- 0
    while((Q - bestQ) > tol){
      Q <- attr(best.tree, "Q-score")

      nni.trees <- lapply(nni(best.tree), ls.tree, D=D)
      nniQ <- sapply(nni.trees, function(x) attr(x, "Q-score"))
      best <- which(nniQ == min(nniQ))
      bestQ <- nniQ[best]

      if(bestQ > Q){
        bestQ <- Q
      } else {
        best.tree <- nni.trees[[best]]
      }
    }
  }
  if(set.neg.to.zero){
    best.tree$edge.length <- ifelse(best.tree$edge.length>0, best.tree$edge.length, 0)
  }
  return(best.tree)
}




#' Least squares phylogenetic inference without distances
#'
#' This function performs phylogenetic inference via least squares with our new loss function. It allows for an exhaustive search among
#' all possible topologies (if not, it will do an NNI search, starting from the NJ tree). This function infers an unrooted tree.
#'
#' @param alignment a nucleotide sequence alignment, of class phyDat
#' @param initvals starting values for each branch
#' @param search.all if TRUE, an exhaustive search across all topologies will be performed. Otherwise, an NNI search will be performed.
#' @param method optimization method for optimx function
#' @param low lower bound for numeric optimization
#' @param high upper bound for numeric optimization
#' @param tol in NNI search, keep searching if improvement is at least this amount
#' @return An unrooted phylogeny
#'
#' @keywords phylogeny, OLS
#' @export
#' @examples
#' phylo.ls.nodist(seq.table)
phylo.ls.nodist <- function(alignment, initvals = NULL, search.all = FALSE, method="nlminb", low=-100, high=2, tol = 1e-10){
  seq.table <- as.character(alignment)
  n <- nrow(seq.table)

  if(search.all){
    all.trees <- allTrees(n,tip.label=row.names(seq.table))
    allQ <- vector()

    # for loop is same speed as using lapply. Not sure if I can speed this up at all.
    for (i in 1:length(all.trees)) {
      output <- new.ls.fit.optimx(all.trees[[i]], seq.table, initvals, method, low, high) # 3/30/20 check this
      all.trees[[i]]$edge.length <- output$par.est
      all.trees[[i]]$ls <- output$ls

      if(output$conv!=0){
        cat('Warning: convergence not reached on at least one potential tree','\n')
      }
      all.trees[[i]]$convergence <- output$conv
      allQ[i] <- output$ls
    }

    best<-which(allQ==min(allQ))
    best.tree <- all.trees[[best]]

  } else {
    # Do nni search

    # Start from NJ tree using distance matrix according to JC69
    data.bin<-as.DNAbin(as.alignment(seq.table))
    D <- as.matrix(dist.dna(data.bin,model="JC69"))
    best.tree <- nj(D)

    output <- new.ls.fit.optimx(my.topology = best.tree, seq.table = seq.table)
    best.tree$edge.length <- output$par.est
    best.tree$ls <- output$ls

    Q <- Inf
    bestQ <- 0
    while((Q - bestQ) > tol){
      Q <- best.tree$ls

      nni.trees <- nni(best.tree)  # not sure if the ordering of nni output is deterministic, so store it first
      nni.output <- lapply(nni.trees, new.ls.fit.optimx, seq.table = seq.table)
      nniQ <- sapply(nni.output, `[[`, "ls")
      best <- which(nniQ == min(nniQ))
      bestQ <- nniQ[best]

      if(bestQ > Q){
        bestQ <- Q
      } else {
        best.tree <- nni.trees[[best]]
        best.tree$edge.length <- nni.output[[best]]$par.est
        best.tree$ls <- nni.output[[best]]$ls
      }
    }
  }
  return(best.tree)
}


#' Maximum likelihood phylogenetic inference
#'
#' This is basically a wrapper function around the pml function from the phangorn package, made to suit our purposes here.
#'
#' @param alignment a nucleotide sequence alignment, of class phyDat
#' @param search.all if TRUE, an exhaustive search across all topologies will be performed. Otherwise, an NNI search will be performed.
#' @param tol in NNI search, keep searching if improvement is at least this amount. Used 1e-8 as default value consistent with phangorn.
#' @return An unrooted phylogeny
#'
#' @keywords phylogeny, OLS
#' @export
#' @examples
#' phylo.ML(seq.table)
phylo.ML <- function(alignment, search.all = FALSE, tol = 1e-8){

  n <- length(alignment)

  if(search.all){
    all.trees <- allTrees(n, tip.label=names(alignment))
    allQ <- vector()

    # 4/2/20 5PM TEST THIS, HASN'T BEEN RUN AT ALL YET
    for (i in 1:length(all.trees)) {
      output <- optim.pml(pml(all.trees[[i]], alignment))
      all.trees[[i]] <- output$tree
      allQ[i] <- output$logLik
    }

    best<-which(allQ==min(allQ))
    best.tree <- all.trees[[best]]

  } else {
    # Do nni search

    # Start from NJ tree using distance matrix according to JC69
    data.bin<-as.DNAbin(as.alignment(seq.table))
    D <- as.matrix(dist.dna(data.bin,model="JC69"))
    best.tree <- nj(D)

    output <- new.ls.fit.optimx(my.topology = best.tree, seq.table = seq.table)
    best.tree$edge.length <- output$par.est
    best.tree$ls <- output$ls

    Q <- Inf
    bestQ <- 0
    while((Q - bestQ) > tol){
      Q <- best.tree$ls

      nni.trees <- nni(best.tree)  # not sure if the ordering of nni output is deterministic, so store it first
      nni.output <- lapply(nni.trees, new.ls.fit.optimx, seq.table = seq.table)
      nniQ <- sapply(nni.output, `[[`, "ls")
      best <- which(nniQ == min(nniQ))
      bestQ <- nniQ[best]

      if(bestQ > Q){
        bestQ <- Q
      } else {
        best.tree <- nni.trees[[best]]
        best.tree$edge.length <- nni.output[[best]]$par.est
        best.tree$ls <- nni.output[[best]]$ls
      }
    }
  }
  return(best.tree)
}
