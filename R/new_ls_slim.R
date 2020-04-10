
#' A phyloLSnoDist function
#'
#' This might not be needed. Revisit later.
#'
#' @param log.branch.length Natural log of the branch lengths at which to calculate the loss
#' @param my.topology A phylogeny in ape format
#' @param seq.dist not sure
#'
#' @keywords phylogeny, OLS
#' @export
#' @examples
#' regular.ls.loss(log.branch.length, my.topology, seq.dist)

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


#' A phyloLSnoDist function
#'
#' Not sure if I need this one anymore either.
#'
#' @param init.brlen Initial branch lengths at which to start the search
#' @param my.topology a phylogeny in ape format
#' @param seq.dist not sure
#'
#' @keywords phylogeny, OLS
#' @export
#' @examples
#' regular.ls.fit(init.brlen, my.topology, seq.dist)


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


#' L2 loss function
#'
#' This function calculates the L2 loss as described in...
#'
#' @param log.branch.length natural log of the branch lengths of the phylogeny
#' @param my.topology a phylogeny in ape format
#' @param seq.dist 4/6/20 not needed?
#' @param regist.matrix 4/6/20 I don't know why this is an argument if it is defined in the function...
#'
#' @keywords phylogeny, OLS
#' @export
#' @examples
#' phylo.ML(seq.table)


# 4/6/20 timecheck this and see if we can speed it up
new.ls.loss = function(log.branch.length, my.topology, seq.table){

  ## first bring branch length to the absolute scale
  branch.length = exp(log.branch.length)

  num.taxa = dim(seq.table)[1]

  my.ls = matrix(0,nrow=num.taxa,ncol=num.taxa)

  my.phylo = my.topology
  my.phylo$edge.length = branch.length

  phylo.dist = cophenetic.phylo(my.phylo)

  ## define JC69 model with its eigen decomposition
  regist.matrix = matrix(1, nrow = 4, ncol = 4) - diag(1, 4)
  jc.69 <- as.eigen(jc.mc(4/3,c(0.25,0.25,0.25,0.25)))

  for (i in 2:num.taxa){
    for (j in 1:(i-1)){
      my.jc = rescale.mc(jc.69, phylo.dist[my.phylo$tip.label[i],my.phylo$tip.label[j]])
      my.ls[i,j] = (phylo.dist[my.phylo$tip.label[i],my.phylo$tip.label[j]] - pair.robust.dist(my.jc, regist.matrix, seq.table[c(my.phylo$tip.label[i],my.phylo$tip.label[j]),]))^2
    }
  }

  return(sum(my.ls))
}


#' Find optimal branch lengths according to L2 loss
#'
#' This function takes a topology and nucleotide sequence alignment, and finds the optimal branch lengths
#' according to the new loss function
#'
#' @param my.topology a phylogenetic topology on which to optimize the branch lengths
#' @param seq.table a nucleotide sequence alignment, of class \code{phyDat}. Use \code{phyDat}
#' function from \code{phangorn} package to transform data if necessary.
#' @param initvals starting values for each branch
#' @param init.brlen initial branch lengths to use in the optimization routine
#' @param method option for optimx function to choose which optimization method to use
#' @param low contraint on optimization. Defaults to -100 which will constrain br len to be min of 3.72e-44
#' @param high constraint on optimization. Defaults to 2 which will constrain br len to be max of e^2 ~ 7.389
#' @param starttests control parameter for optimx, default to FALSE for speed
#' @param kkt control parameter for optimx, default to FALSE for speed
#'
#' @return An unrooted phylogeny
#'
#' @keywords phylogeny, OLS
#' @export
#' @examples
#' phylo.ML(seq.table)


# order of arguments changed 4/1/20, might have implications on old code...
new.ls.fit.optimx <- function(my.topology, seq.table, init.brlen=NULL, method="nlminb", low=-100, high=NULL, starttests=TRUE, kkt=TRUE){

  if(class(seq.table)!="phyDat"){
    cat('Error: alignment must be of class phyDat')
  } else{

    # If not otherwise specified, use 0.1 as starting value for all branches
    if(is.null(init.brlen)){
      init.brlen <- rep(0.1, dim(my.topology$edge)[1])
    }

    # If not otherwise specified, use max of pairwise distances as the upper limit
    if(is.null(high)){
      pair.dists <- dist.dna(as.DNAbin(seq.table))
      high <- log(max(pair.dists))
    }

    seq.table <- read.phylosim.nuc(as.character(seq.table))

    return.val = NULL
    par.num = length(init.brlen)
    #  regist.mat = matrix(1,4,4) - diag(rep(1,4))
    optim.out<-list(par=1,convcode=1)
    count<-0

    # try up to 10 times to reach convergence
    while((optim.out$convcode != 0) & count<10){

      optim.out <- optimx(
        log(init.brlen),
        new.ls.loss,
        lower = rep(low,par.num),
        upper = rep(high,par.num),
        method = method,
        my.topology = my.topology,
        seq.table = seq.table,
        control=list(starttests=starttests, kkt=kkt)
      )

      # Keep track of iterations
      count<-count+1

      # If convergence was not reached, re-initialize starting values to random values
      init.brlen<-runif(par.num,0,0.1)  # 4/10/20 changed from 0.5 to 0.1 because I suspect that smaller starting values will be better
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

}


#' Maximum likelihood phylogenetic inference
#'
#' This is basically a wrapper around the pml function from the phangorn package, made to suit our purposes here
#' and perform either an exhaustive search through all topologies or an NNI search with pml as the workhorse.
#' I don't really recommend that anyone should use this function instead of optim.pml from phangorn unless you see a strong reason to.
#'
#' @param alignment a nucleotide sequence alignment, of class \code{phyDat}. Use \code{phyDat}
#' function from \code{phangorn} package to transform data if necessary.
#' @param initvals starting values for each branch
#' @param search.all if TRUE, an exhaustive search across all topologies will be performed. Otherwise, an NNI search will be performed.
#' @param tol in NNI search, keep searching if improvement is at least this amount. Used 1e-8 as default value consistent with phangorn.
#' @return An unrooted phylogeny
#'
#' @keywords phylogeny, OLS
#' @export
#' @examples
#' phylo.ML(seq.table)


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
#' @param alignment a nucleotide sequence alignment, of class \code{phyDat}. Use \code{phyDat}
#' function from \code{phangorn} package to transform data if necessary.
#' @param initvals starting values for each branch
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
  data.bin<-as.DNAbin(alignment)
  D <- as.matrix(dist.dna(data.bin,model=model))
  n <- nrow(D)

#    if(class(D) == "dist"){
#      D <- as.matrix(D)
#    }
#    n <- nrow(D)

  if(search.all){
    # Do search through all possible topologies
    all.trees <- allTrees(n, tip.label = row.names(D)) # change this (as.matrix(D) is redundant)
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
#' @param alignment a nucleotide sequence alignment, of class \code{phyDat}. Use \code{phyDat}
#' function from \code{phangorn} package to transform data if necessary.
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
#'

# 4/9/20 not sure if bound_rm is a good idea. try the JC69 constraint idea and see if that does it.
phylo.ls.nodist <- function(alignment, initvals = NULL, search.all = FALSE, method="nlminb", low=-100, high=NULL, tol = 1e-10, bound_rm = TRUE){
#  seq.table <- as.character(alignment)   # keep it as phyDat

  if(class(alignment) != "phyDat"){
    cat('Error: alignment must be of class phyDat')
  } else {

    # If not otherwise specified, use max of pairwise distances as the upper limit
    if(is.null(high)){
      pair.dists <- dist.dna(as.DNAbin(alignment))
      high <- log(max(pair.dists))
    }

    if(search.all){
      all.trees <- allTrees(length(my.align),tip.label=names(alignment))
      allQ <- vector()

      # for loop is same speed as using lapply. Not sure if I can speed this up at all.
      for (i in 1:length(all.trees)) {
        output <- new.ls.fit.optimx(all.trees[[i]], alignment, initvals, method, low, high) # 3/30/20 check this
        all.trees[[i]]$edge.length <- output$par.est
        all.trees[[i]]$ls <- output$ls

        if(output$conv!=0){
          cat('Warning: convergence not reached on at least one potential tree','\n')
        }
        all.trees[[i]]$convergence <- output$conv
        allQ[i] <- output$ls
      }

      if(bound_rm){
        hit.bound <- apply(sapply(all.trees, `[[`, "edge.length")==exp(2), 2, sum) + apply(sapply(all.trees, `[[`, "edge.length")==exp(-100), 2, sum)
        all.trees <- all.trees[hit.bound==0]
        allQ <- allQ[hit.bound==0]
      }

      best<-which(allQ==min(allQ))
      best.tree <- all.trees[[best]]
      if(best.tree$convergence!=0){
        cat('Warning: convergence not reached on best tree.')
      }

    } else {
      # Do nni search

      # Start from NJ tree using distance matrix according to JC69
      data.bin<-as.DNAbin(alignment)
      D <- as.matrix(dist.dna(data.bin,model="JC69"))
      best.tree <- nj(D)

      output <- new.ls.fit.optimx(my.topology = best.tree, seq.table = alignment)
      best.tree$edge.length <- output$par.est
      best.tree$ls <- output$ls

      Q <- Inf
      bestQ <- 0
      while((Q - bestQ) > tol){
        Q <- best.tree$ls

        nni.trees <- nni(best.tree)  # not sure if the ordering of nni output is deterministic, so store it first
        nni.output <- lapply(nni.trees, new.ls.fit.optimx, seq.table = alignment)
        nniQ <- sapply(nni.output, `[[`, "ls")


        if(bound_rm){
          hit.bound <- apply(sapply(nni.trees, `[[`, "edge.length")==exp(2), 2, sum) + apply(sapply(nni.trees, `[[`, "edge.length")==exp(-100), 2, sum)
          nni.trees <- nni.trees[hit.bound==0]
          nniQ <- nniQ[hit.bound==0]
        }

        if(length(nniQ)==0){
          bestQ <- Q
        } else {
          best <- which(nniQ == min(nniQ))
          bestQ <- nniQ[best]
          if(nni.output[[best]]$conv!=0){
            cat('Warning: convergence not reached on best tree.')
          }

          if(bestQ > Q){
            bestQ <- Q
          } else {
            best.tree <- nni.trees[[best]]
            best.tree$edge.length <- nni.output[[best]]$par.est
            best.tree$ls <- nni.output[[best]]$ls

          }
        }
      }
    }
    return(best.tree)

  }
}




#' Maximum likelihood phylogenetic inference
#'
#' This is basically a wrapper around the pml function from the phangorn package, made to suit our purposes here
#' and perform either an exhaustive search through all topologies or an NNI search with pml as the workhorse.
#' I don't really recommend that anyone should use this function instead of optim.pml from phangorn unless you see a strong reason to.
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

    for (i in 1:length(all.trees)) {
      all.trees[[i]]$edge.length <- rep(0.1, dim(all.trees[[i]]$edge)[1])
      output <- optim.pml(pml(all.trees[[i]], alignment), control=pml.control(trace=0))
      all.trees[[i]] <- output$tree
      all.trees[[i]]$logLik <- output$logLik
      allQ[i] <- output$logLik
    }

    best<-which(allQ==max(allQ))
    best.tree <- all.trees[[best]]

  } else {
    # Do nni search

    # Start from NJ tree using distance matrix according to JC69
    data.bin<-as.DNAbin(alignment)
    D <- as.matrix(dist.dna(data.bin,model="JC69"))
    best.tree <- nj(D)

    output <- optim.pml(pml(best.tree, alignment), control=pml.control(trace=0)) # this might be what messes things up...
    best.tree <- output$tree
    best.tree$logLik <- output$logLik

    Q <- Inf
    bestQ <- 0
    while((Q - bestQ) > tol){
      Q <- best.tree$logLik

      nni.trees <- nni(best.tree)  # not sure if the ordering of nni output is deterministic, so store it first
      nni.output <- lapply(nni.trees, pml, data = alignment)
      nniQ <- sapply(nni.output, `[[`, "logLik")
      best <- which(nniQ == max(nniQ))
      bestQ <- nniQ[best]

      if(bestQ < Q){
        bestQ <- Q
      } else {
        best.tree <- nni.trees[[best]]
        best.tree <- nni.output[[best]]$tree
        best.tree$logLik <- nni.output[[best]]$logLik
      }
    }
  }
  return(best.tree)
  # 4/2/20 I think it's done but could use a little more testing
}


#' Least Squares K80 without distances
#'
#' This function calculates the distance free loss function under the Kimura 2 parameter model (K80).
#'
#' @param log.br.len natural log of the branch lengths
#' @param log.kappa natural log of the transition/transversion ratio
#' @param my.topology the tree on which to calculate the loss function
#' @param seq.table a nucleotide sequence alignment
#'
#' @keywords phylogeny, OLS
#' @export
#' @examples
#' new.loss.K80(log.br.len, kappa, my.topology, seq.table)

new.loss.K80 = function(log.params, my.topology, seq.table){

  ## first bring branch length to the absolute scale
  n.br<-length(log.params) - 1
  branch.length = exp(log.params[1:n.br])
  ts.tv.ratio <- exp(log.params[n.br+1])

  ts.matrix = matrix(c(0,1,0,0,1,0,0,0,0,0,0,1,0,0,1,0), 4,4,byrow=TRUE)
  tv.matrix = matrix(1,4,4) - ts.matrix - diag(rep(1,4))

  num.taxa = dim(seq.table)[1]

  my.ls = 0

  my.phylo = my.topology
  my.phylo$edge.length = branch.length

  phylo.dist = cophenetic.phylo(my.phylo)

  for (i in 2:num.taxa){
    for (j in 1:(i-1)){
      ## define K80 model with its eigen decomposition
      my.K80 = as.eigen.hky(c(ts.tv.ratio,1), c(0.25,0.25,0.25,0.25), scale=T)
      my.K80 = rescale.mc(my.K80, phylo.dist[my.phylo$tip.label[i],my.phylo$tip.label[j]])

      tree.ts.dist = pair.conv.dist(my.K80, ts.matrix)
      tree.tv.dist = pair.conv.dist(my.K80, tv.matrix)

      robust.ts.dist = pair.robust.dist(my.K80, ts.matrix, seq.table[c(my.phylo$tip.label[i],my.phylo$tip.label[j]),])
      robust.tv.dist = pair.robust.dist(my.K80, tv.matrix, seq.table[c(my.phylo$tip.label[i],my.phylo$tip.label[j]),])

      my.ls = my.ls + (tree.ts.dist - robust.ts.dist)^2 + (tree.tv.dist - robust.tv.dist)^2

    }
  }
  return(my.ls)
}


#' Least Squares K80 without distances
#'
#' This function calculates the distance free loss function under the Kimura 2 parameter model (K80).
#'
#' @param log.br.len natural log of the branch lengths
#' @param log.kappa natural log of the transition/transversion ratio
#' @param my.topology the tree on which to calculate the loss function
#' @param seq.table a nucleotide sequence alignment
#'
#' @keywords phylogeny, OLS
#' @export
#' @examples
#' new.loss.K80(log.br.len, kappa, my.topology, seq.table)

new.ls.fit.K80 <- function(my.topology, seq.table, init.brlen = NULL, init.kappa = NULL, method="nlminb", low=-100, high=2, high.k = 100, starttests=FALSE, kkt=FALSE){

  # change from simSeq output to character matrix if needed
  if(typeof(seq.table) == "list"){
    seq.table <- as.character(seq.table)
  }

  # change A,C,T,G to numeric if needed
  if(typeof(seq.table) == "character"){
    seq.table <- read.phylosim.nuc(seq.table)
  }


  n.br <- dim(seq.table)[1]*2 - 3

  if(is.null(init.brlen)){
    init.brlen <- rep(0.1, n.br)
  }

  if(is.null(init.kappa)){
    init.kappa <- 1
  }

  return.val = NULL
  par.num = n.br + 1
  optim.out<-list(par=1,convcode=1)
  count<-0

  # try up to 10 times to reach convergence
  while((optim.out$convcode != 0) & count<10){

    optim.out <- optimx(
      log(c(init.brlen, init.kappa)),
      new.loss.K80,
      lower = rep(low,par.num),
      upper = c(rep(high,n.br),high.k),
      method = method,
      my.topology = my.topology,
      seq.table = seq.table,
      control=list(starttests=starttests, kkt=kkt)
    )

    # Keep track of iterations
    count<-count+1

    # If convergence was not reached, re-initialize starting values to random values
    init.brlen <- runif(n.br,0,0.5)
    init.kappa <- runif(1, 0, 20)
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


#' K80 eigen decomposition
#'
#' This function defines the K80 model via its eigen decomposition.
#'
#' @param hky.rates transition and transversion rates
#' @param mc.stat stationary distribution
#' @param scale if TRUE then scale the transition matrix. Defaults to FALSE.
#'
#' @keywords phylogeny, OLS
#' @export
#' @examples
#' as.eigen.hky(hky.rates, mc.stat, scale = F)

as.eigen.hky = function(hky.rates, mc.stat, scale = F){

  transition.rate = hky.rates[1]
  transversion.rate = hky.rates[2]

  left.eigen = matrix(rep(0,16), nrow = 4, ncol = 4)
  right.eigen = matrix(rep(0,16), nrow = 4, ncol = 4)

  pi.Y = mc.stat[3] + mc.stat[4]
  pi.R = mc.stat[1] + mc.stat[2]

  left.eigen[,1] = rep(1,4)
  left.eigen[,2] = c(-1/pi.R,-1/pi.R, 1/pi.Y, 1/pi.Y)
  left.eigen[,3] = c(mc.stat[2]/pi.R, -mc.stat[1]/pi.R, 0, 0)
  left.eigen[,4] = c(0, 0, -mc.stat[4]/pi.Y, mc.stat[3]/pi.Y)

  right.eigen[,1] = mc.stat
  right.eigen[,2] = c(-pi.Y*mc.stat[1],-pi.Y*mc.stat[2], pi.R*mc.stat[3], pi.R*mc.stat[4])
  right.eigen[,3] = c(1, -1, 0, 0)
  right.eigen[,4] = c(0, 0, -1, 1)

  ##  right.eigen[,1] = c(mc.stat[4], mc.stat[3],mc.stat[1], mc.stat[2])
  ##  right.eigen[,2] = c(pi.R*mc.stat[4],pi.R*mc.stat[3], -pi.Y*mc.stat[1], -pi.Y*mc.stat[2])
  ##  right.eigen[,3] = c(0,0,-1, -1)
  ##  right.eigen[,4] = c(1,-1,0, 0)

  eigen.val = c(0, -transversion.rate, -(pi.Y*transversion.rate + pi.R*transition.rate),
                -(pi.Y*transition.rate + pi.R*transversion.rate))

  revmc.eigen = list(hky.mc(hky.rates, mc.stat, scale = F)[[1]],
                     mc.stat, eigen.val,left.eigen,
                     t(right.eigen))
  names(revmc.eigen) = c("rate.matrix", "stationary",
                         "values", "vectors", "invvectors")

  class(revmc.eigen) = c("revmc", "eigen")

  if (scale){
    revmc.eigen = rescale.mc(revmc.eigen, 1.0/sum(-mc.stat*diag(revmc.eigen$rate.matrix)))
  }

  return(revmc.eigen)

}



# This is a problem...
#Error: $ operator is invalid for atomic vectors
#In addition: Warning messages:
# 1: In if (bestQ > Q) { :
#    the condition has length > 1 and only the first element will be used
# 2: In if (bestQ > Q) { :
#    the condition has length > 1 and only the first element will be used
# 3: In phy$tip.label <- attr(x, "TipLabel") : Coercing LHS to a list



