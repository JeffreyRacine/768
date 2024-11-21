  ## This function returns the weight matrix for a local polynomial,
  ## and was rewritten 14/1/15 in Toulouse while visiting JP. It
  ## supports mixed data types. Basic error checking is
  ## undertaken. deriv = 0, strips off weights for mean, = p partials
  ## up to order 2. No cross-partials in this one. Basically useful
  ## for univariate case when deriv > 0 though could be refined - the
  ## old function was slower but had more capability (that basically
  ## went unused).

  ## Update - from ?npksum, "The option permutation.operator= can be
  ## used to ‘mix and match’ operator strings to create a ‘hybrid’
  ## kernel, in addition to the kernel sum with no operators applied,
  ## one for each continuous dimension in the data. For example, for a
  ## two-dimensional data frame of numeric datatypes,
  ## permutation.operator=c("derivative") will return the usual kernel
  ## sum as if operator = c("normal","normal") in the ksum member, and
  ## in the p.ksum member, it will return kernel sums for operator =
  ## c("derivative","normal"), and operator =
  ## c("normal","derivative"). This makes the computation of gradients
  ## much easier."

  ## So, the upshot is that I could, for the multivariate case, add
  ## the derivative stuff.

  Hmat.lp <- function(deriv=0,
                      mydata.train=NULL,
                      mydata.eval=NULL,
                      bws=NULL,
                      p=0,
                      shrink=TRUE,
                      warning.immediate=TRUE,
                      ...) {

      ## 14/1/15, Toulouse - note that the weights herein ** DO NOT **
      ## shrink towards the lc estimator (neither for the function nor
      ## derivatives), unlike the function returned in
      ## glpreg(). However, they all appear to agree with the previous
      ## Hmat.lp with ** also ** did not shrink towards the lc
      ## estimator. This is noticeably faster, which ought to render
      ## Tikhonov faster as well.

      ## Basic error checking...

      if(is.null(mydata.train)) stop("You must provide training data")
      if(is.null(mydata.eval)) mydata.eval <- mydata.train
      if(is.null(bws)) stop("You must provide bandwidths")

      n.train=nrow(mydata.train)
      n.eval=nrow(mydata.eval)

      X.train <- as.data.frame(mydata.train)
      X.eval <- as.data.frame(mydata.eval)

      ## Check whether it appears that training and evaluation data are
      ## conformable...

      if(ncol(X.train)!=ncol(X.eval))
          stop("Error: training and evaluation data have unequal number of columns\n")

      X.col.numeric <- sapply(1:ncol(X.train),function(i){is.numeric(X.train[,i])})

      ## k represents the number of numeric regressors, this will return
      ## zero if there are none

      k <- ncol(as.data.frame(X.train[,X.col.numeric]))

      if(k > 0) {
          X.train.numeric <- as.data.frame(X.train[,X.col.numeric])
          X.eval.numeric <- as.data.frame(X.eval[,X.col.numeric])
      }

      if(deriv<0||deriv>2)
          stop(paste("Error: deriv= (integer) is invalid\n[min = ", 0, ", max = ",  p, "]\n",sep=""))

      if(p < 0)
          stop(paste("Error: p (order of polynomial) must be a non-negative integer\np is (", p, ")\n",sep=""))

      K.x <- npksum(txdat=X.train,
                    exdat=X.eval,
                    bws=bws,
                    return.kernel.weights=TRUE,
                    ...)$kw

      if(p==0) {

          ## No shrinking necessary for local constant estimator

          if(deriv==0) {

              Hmat <- t(K.x)/np:::NZD(rowSums(t(K.x)))

          } else if(deriv==1) {

              ## Note this is not general XXX Feb 25 2015, for
              ## univariate z only

              K.x.deriv <- npksum(txdat=X.train,
                                  exdat=X.eval,
                                  bws=bws,
                                  return.kernel.weights=TRUE,
                                  operator="derivative",
                                  ...)$kw/np:::NZD(bws)

              rSk <- np:::NZD(rowSums(t(K.x)))

              Hmat <- t(K.x.deriv)/rSk-t(K.x)/rSk*(rowSums(t(K.x.deriv))/rSk)

          }

      }

      if(p > 0) {

          ## Re-use this matrix, shrinking occurs here

          W.z <- crs:::W.glp(xdat=X.train.numeric,
                             degree=rep(p,NCOL(X.train.numeric)))

          if(is.null(mydata.eval)) {
              ## Guess we could avoid copy with conditional statement below using either W.z or W.z.eval
              W.z.eval <- W.z
          } else {
              W.z.eval <- crs:::W.glp(xdat=X.train.numeric,
                                      exdat=as.data.frame(X.eval.numeric),
                                      degree=rep(p,NCOL(X.train.numeric)))
          }

          nc <- ncol(W.z)

          WzkWz.inv <- list()

          for(i in 1:ncol(K.x)) {

              if(tryCatch(WzkWz.inv[[i]] <- as.matrix(chol2inv(chol(t(W.z)%*%(K.x[,i]*W.z)))),
                          error = function(e){
                              return(matrix(FALSE,nc,nc))
                          })[1,1]!=FALSE) {

              } else {

                  if(shrink==FALSE) {

                      ## If we do not explicitly engage ridging then we do not fail
                      ## and terminate, rather, we return NA when Wmat.sum is
                      ## singular

                      Hmat <- NA

                  } else {

                      ## Ridging

                      epsilon <- 1/n.train
                      ridge <- 0

                      while(tryCatch(as.matrix(chol2inv(chol((chol(t(W.z)%*%(K.x[,i]*W.z)+diag(rep(ridge,nc))))))),
                                     error = function(e){
                                         return(matrix(FALSE,nc,nc))
                                     })[1,1]==FALSE) {
                          ridge <- ridge + epsilon
                      }

                      WzkWz.inv[[i]] <- as.matrix(chol2inv(chol(t(W.z)%*%(K.x[,i]*W.z)+diag(rep(ridge,nc)))))

                      warning(paste("Ridging obs. ", i, ", ridge = ", signif(ridge,6),sep=""),
                              immediate.=warning.immediate,
                              call.=!warning.immediate)

                  }

              }

          }
      }

      if(p==1) {

          if(deriv==0) Hmat <- t(sapply(1:ncol(K.x),function(i){W.z.eval[i,,drop=FALSE]%*%WzkWz.inv[[i]]%*%t(W.z)*K.x[,i]}))
          if(deriv==1) {
              W.z.deriv.1 <- crs:::W.glp(xdat=X.train.numeric,
                                         exdat=as.matrix(X.eval.numeric),
                                         degree=rep(p,NCOL(X.train.numeric)),
                                         gradient.vec = 1)

              Hmat <- t(sapply(1:ncol(K.x),function(i){W.z.deriv.1[i,,drop=FALSE]%*%WzkWz.inv[[i]]%*%t(W.z)*K.x[,i]}))
          }

      }

      if(p >= 2) {

          if(deriv==0) {

              Hmat <- t(sapply(1:ncol(K.x),function(i){W.z.eval[i,,drop=FALSE]%*%WzkWz.inv[[i]]%*%t(W.z)*K.x[,i]}))

          }


          if(deriv==1) {

              W.z.deriv.1 <- crs:::W.glp(xdat=X.train.numeric,
                                         exdat=as.matrix(X.eval.numeric),
                                         degree=rep(p,NCOL(X.train.numeric)),
                                         gradient.vec = 1)

              Hmat <- t(sapply(1:ncol(K.x),function(i){W.z.deriv.1[i,,drop=FALSE]%*%WzkWz.inv[[i]]%*%t(W.z)*K.x[,i]}))

          }

          if(deriv==2) {

              W.z.deriv.2 <- crs:::W.glp(xdat=X.train.numeric,
                                         exdat=as.matrix(X.eval.numeric),
                                         degree=rep(p,NCOL(X.train.numeric)),
                                         gradient.vec = 2)

              Hmat <- t(sapply(1:ncol(K.x),function(i){W.z.deriv.2[i,,drop=FALSE]%*%WzkWz.inv[[i]]%*%t(W.z)*K.x[,i]}))

          }

      }

      return(Hmat)
  }

