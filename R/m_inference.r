
#' Title
#'
#' @param x 
#' @param pfn 
#' @param qfn 
#' @param rt.lv 
#'
#' @return
#' @export
#'
#' @examples
rtlv.basic <- function(x, pfn, qfn, rt.lv=c(2,5,10,20,50,100)){
  tmp.rtlv <- t(apply(X = x, MARGIN = 2, FUN = function(data){
    total.size <- length(data)
    s.x <- data[!is.na(data)]
    qfn(1-1/rt.lv, pfn(samlmu(s.x, 3)))
  }))
  colnames(tmp.rtlv) <- paste0("rt",rt.lv)
  tmp.rtlv
}

exceedence.prob.basic <- function(x, pfn, cdfn, z=seq(50, 500, 50)){
  tmp <- t(apply(X = x, MARGIN = 2, FUN = function(data){
    total.size <- length(data)
    s.x <- data[!is.na(data)]
    param <- pfn(samlmu(s.x, 3))
    1-cdfn(z, param)
  }))
  colnames(tmp) <- paste0("z",z)
  tmp
}

#' Title
#'
#' @param data.obj variable 이름으로 묶인 리스트, 리스트의 요소는 컬럼이 station으로 구분되는 data.frame 이어야 한다.
#' @param pfns Parameter estimation Functions, data.obj의 variable 개수와 동일해야 함.
#' @param qfns Quantile Functions, data.obj의 variable 개수와 동일해야 함.
#' @param rt.lv 
#'
#' @return
#' @export
#'
#' @examples
rtlv.multi.var <- function(data.obj, pfns, qfns, rt.lv=c(2,5,10,20,30,50,100)){
  # usethis::use_package("pbapply")
  
  if( typeof(pfns) !="list" ) stop("'pfns' must be 'list' type")
  if( typeof(qfns) !="list" ) stop("'qfns' must be 'list' type")
  
  if( (length(data.obj) != length(pfns)) || (length(pfns) != length(qfns)) ){
    stop('Must satisfy the following: "length of variables == length of pfns == length of qfns"')
  }
  
  vars <- names(data.obj)
  names(pfns) <- vars
  names(qfns) <- vars
  
  res <- 
    pbapply::pblapply(vars, function(v_){
      x_ <- data.obj[[v_]]
      pfn_ <- pfns[[v_]]
      qfn_ <- qfns[[v_]]
      if(!is.null(dim(x_))){ # x_ 의 dim이 null 이 아니면 더이상 하위 리스트가 없음. (ex. obs 자료)
        rtlv.basic(x=x_, pfn=pfn_, qfn=qfn_, rt.lv=rt.lv)
        
      }else{ # x_ 의 dim이 null 이면 하위 리스트가 있음. (ex. model scenario 자료)
        lapply(x_, function(model.x_){
          rtlv.basic(x=model.x_, pfn=pfn_, qfn=qfn_, rt.lv=rt.lv)
        })  
      }
    })
  names(res) <- vars
  res
}


#' Title
#'
#' @param ref 
#' @param model 
#' @param target.rtlv 
#'
#' @return
#' @export
#'
#' @examples
tilde.normalization <- function(ref, model, target.rtlv){
  # usethis::use_package("pbapply")
  
  if(!all(names(ref) == names(model))) stop('The variables order of "ref" and "model" should be the same.')  
  
  # Calc 'Normalization Factors'
  cat('Calc "Normalization Factors"\n')
  norm.factors <- 
    pbapply::pbmapply(function(ref_, model_){
      total.df <- rbind(ref_, do.call(rbind, model_))
      normalize.factor <- 
        data.frame(med=apply(total.df, 2, median, na.rm=T),
                   max=apply(total.df, 2, max, na.rm=T),
                   min=apply(total.df, 2, min, na.rm=T))
      normalize.factor
    }, ref, model, SIMPLIFY = FALSE)
  
  cat('Normalize Target\n')
  
  # target.rtlv[[1]] 의 dim이 null 이 아니면 더이상 하위 리스트가 없음. (ex. obs 자료)
  if(!is.null(dim(target.rtlv[[1]]))){ 
    mapply(function(target.rtlv_, norm.factors_){
      tmp <- 
        pbapply::pblapply(X = 1:nrow(target.rtlv_),
                          FUN = function(i.loc){
                            rt.T <- target.rtlv_[i.loc,]
                            med_ <- norm.factors_$med[i.loc]
                            min_ <- norm.factors_$min[i.loc]
                            max_ <- norm.factors_$max[i.loc]
                            
                            norm.value <- rt.T
                            norm.value[rt.T >= med_] <- (norm.value[rt.T >= med_] - med_)/(max_ - med_)
                            norm.value[rt.T <  med_] <- (norm.value[rt.T < med_] - med_)/(med_ - min_)
                            norm.value
                          })
      tmp_ <- do.call(rbind, tmp)
      rownames(tmp_) <- paste0("G",1:nrow(tmp_))
      tmp_
    }, target.rtlv, norm.factors, SIMPLIFY = FALSE)
  }else{ # target.rtlv[[1]] 의 dim이 null 이면 하위 리스트가 있음. (ex. model scenario 자료)
    mapply(function(multi.target.rtlv_, norm.factors_){
      pbapply::pblapply(multi.target.rtlv_, function(target.rtlv_){
        tmp <-
          lapply(X = 1:nrow(target.rtlv_),
                 FUN = function(i.loc){
                   rt.T <- target.rtlv_[i.loc,]
                   med_ <- norm.factors_$med[i.loc]
                   min_ <- norm.factors_$min[i.loc]
                   max_ <- norm.factors_$max[i.loc]
                   
                   norm.value <- rt.T
                   norm.value[rt.T >= med_] <- (norm.value[rt.T >= med_] - med_)/(max_ - med_)
                   norm.value[rt.T <  med_] <- (norm.value[rt.T < med_] - med_)/(med_ - min_)
                   norm.value
                 })
        tmp_ <- do.call(rbind, tmp)
        rownames(tmp_) <- paste0("G",1:nrow(tmp_))
        tmp_
      })
    }, target.rtlv, norm.factors, SIMPLIFY = FALSE)
  }
}


#' Title
#'
#' @param ref 
#' @param model 
#' @param simul.N 
#' @param random.seed 
#' @param candi.sigma.d 
#'
#' @return
#' @export
#'
#' @examples
calc.sigma.d <- function(ref, model, simul.N=10000, random.seed=123456, significance.level=0.05, candi.sigma.d=seq(1.0, 0.01, -0.01)){
  
  if(!all(names(ref) == names(model))) stop('The variables order of "ref" and "model" should be the same.')
  if(length(candi.sigma.d)<2) stop("'candi.sigma.d' : The number of candidates for sigma.D is too small.")
  
  vars <- names(ref)
  
  # usethis::use_package("pbapply")
  # usethis::use_package("MCMCpack")
  # usethis::use_package("dplyr")
  require(dplyr)
  
  set.seed(random.seed)
  N <- simul.N
  model.len <- get.model.length(model)
  
  equal.weight <- rep(1/model.len, model.len)
  simul.weight <- MCMCpack::rdirichlet(N, alpha=rep(2, model.len))
  simul.t_ <- apply((simul.weight - replicate(N, equal.weight) %>% t())^2 / replicate(N, equal.weight) %>% t(), 1, sum)
  
  
  chisq.result.with.sigmad.res <-
    pbapply::pblapply(candi.sigma.d, function(sigma.d){
      
      perf.weight.avg <- 
        mapply(function(ref_, model_){
         
          p.k <-
            lapply(model_, function(model__){
              dk.square <- apply((model__-ref_)^2, 1, sum)
              exp(-dk.square/(sigma.d)^2)
            })
          
          pi.weight.nume <- do.call(rbind, p.k)
          pi.weight.nume <- data.frame(model=rownames(pi.weight.nume), pi.weight.nume) %>% dplyr::arrange(model)
          
          perf.weight <- 
            apply(pi.weight.nume[,-1], 2, function(loc.data){
              loc.data/sum(loc.data)
            })
          perf.weight.df <- data.frame(model=pi.weight.nume$model, perf.weight)
          
          perf.weight.df.mean <- apply(perf.weight.df[,-1], 1, median)
          
          t_ <- sum((perf.weight.df.mean - equal.weight)^2 / equal.weight)
          
          p.value <- sum(simul.t_ >= t_)/N
          
          data.frame(stats=t_, p.value=p.value)
        }, ref, model, SIMPLIFY = FALSE)
      perf.weight.avg
      
    })
  var.agg <- 
    lapply(vars, function(v_){
      
      chisq.result.with.sigmad <- do.call(rbind, lapply(chisq.result.with.sigmad.res, function(x) x[[v_]]))
      chisq.result.with.sigmad <- 
        data.frame(sigma.d=candi.sigma.d, chisq.result.with.sigmad) %>% 
        mutate(result=ifelse(p.value>=significance.level, "H0 accept", "H0 reject"), group="A")
      
      optim.sigma.d.chisq <- 
        chisq.result.with.sigmad[tail(which(chisq.result.with.sigmad$result=="H0 accept"),1)+1,]
    }) %>% do.call(what=rbind)
  print(var.agg)
  sigma.D <- round(mean(var.agg$sigma.d),2)
  print(paste0("Avg of sigma.D = ", sigma.D))
  return(list(simul.res=chisq.result.with.sigmad.res, var.agg=var.agg, sigma.D=sigma.D))
}


#' Title
#'
#' @param ref 
#' @param model 
#' @param sigma.D 
#'
#' @return
#' @export
#'
#' @examples
calc.performance.weights <- function(ref, model, sigma.D){
  
  if(!all(names(ref) == names(model))) stop('The variables order of "ref" and "model" should be the same.')
  
  # usethis::use_package("abind")
  # usethis::use_package("dplyr")
  # usethis::use_package("pbapply")
  # 
  
  model.names <- get.model.names(model)
  
  perf.weight.avg <- 
    pbapply::pbmapply(function(ref_, model_){
      
      p.k <- 
        lapply(model_, function(model__){
          dk.square <- apply((model__-ref_)^2, 1, sum)
          exp(-dk.square/sigma.D^2)
        })
      pi.weight.nume <- do.call(rbind, p.k)
      pi.weight.nume <- data.frame(model=rownames(pi.weight.nume), pi.weight.nume) %>% dplyr::arrange(model)
      
      perf.weight <- 
        apply(pi.weight.nume[,-1], 2, function(loc.data){
          loc.data/sum(loc.data)
        })
      perf.weight
    }, ref, model, SIMPLIFY = FALSE) %>% abind::abind(along = 3) %>% apply(MARGIN = 1:2, mean)
  
  
  perf.weight.df <- data.frame(model=sort(model.names), perf.weight.avg)
  
  perf.weight.list <- 
    lapply(seq_along(model.names), function(k){
      t(perf.weight.df[k,-1])
    })
  
  names(perf.weight.list) <- model.names
  
  return(list(weight.avg=perf.weight.avg, weight.df=perf.weight.df, weight.list=perf.weight.list))
}

#' Title
#'
#' @param perf.weight_ 
#' @param indep.weight_ 
#'
#' @return
#' @export
#'
#' @examples
calc.pi.weight <- function(perf.weight_, indep.weight_){
  
  if(length(perf.weight_$model) != length(indep.weight_$model)){
    stop("The number of models must be the same.")
  }
  
  if(!all(perf.weight_$model==indep.weight_$model)){
    stop("The order of models must be the same.")
  }
  
  # model.len <- length(perf.weight_$model)
  # if( sum(mapply(function(a,b) a==b , perf.weight_$model, indep.weight_$model, SIMPLIFY = T)) != model.len){
  #   stop("You must check the order of models")
  # }
  
  # Multiplication of Perf.weights and Indep.Weights
  pi.weight_ <- perf.weight_[,-1] * replicate(ncol(perf.weight_)-1, indep.weight_[,2])
  
  # Normalize so that the sum of weights is 1.
  pi.weight <- apply(pi.weight_, 2, function(x){
    x/sum(x)
  })
  
  data.frame(model=indep.weight_$model, pi.weight)
}

#' Title
#'
#' @param x 
#' @param pfn 
#' @param ref.rtlv 
#' @param rt.lv 
#'
#' @return
#' @export
#'
#' @examples
return.period.basic <- function(x, pfn, ref.rtlv, rt.lv){
  require(dplyr)
  tmp <- # model 자료의 gev 파라메터 추정 
    apply(X = x, MARGIN = 2, FUN = function(g){
      pfn(lmom::samlmu(x = g, nmom = 3))
    })
  
  tmp2 <- # obs 자료 merge
    data.frame(loc=rownames(t(tmp)), t(tmp)) %>% 
    dplyr::arrange(loc) %>%
    dplyr::left_join(ref.rtlv, by = c("loc"="loc"))
  
  tmp3 <- # return period 추정
    apply(tmp2, 1, function(r_){
      loc_ <- r_[1]
      gev.para_ <- r_[2:4]
      ref.rtlvs <- as.numeric(r_[-1:-4])
      1/(1-lmom::cdfgev(ref.rtlvs, as.numeric(gev.para_)))
    })
  tmp.rs <- data.frame(loc=tmp2$loc, t(tmp3))
  colnames(tmp.rs) <- c("loc", paste0("rt",rt.lv))
  tmp.rs
}


exceedence.prob.multi.var <- function(data.obj, pfns, cdfns, z){
  # usethis::use_package("pbapply")
  
  if( typeof(pfns) !="list" ) stop("'pfns' must be 'list' type")
  if( typeof(cdfns) !="list" ) stop("'cdfns' must be 'list' type")
  
  if( (length(data.obj) != length(pfns)) || (length(pfns) != length(cdfns)) ){
    stop('Must satisfy the following: "length of variables == length of pfns == length of cdfns"')
  }
  
  vars <- names(data.obj)
  names(pfns) <- vars
  names(cdfns) <- vars
  
  res <- 
    pbapply::pblapply(vars, function(v_){
      x_ <- data.obj[[v_]]
      pfn_ <- pfns[[v_]]
      cdfn_ <- cdfns[[v_]]
      if(!is.null(dim(x_))){ # x_ 의 dim이 null 이 아니면 더이상 하위 리스트가 없음. (ex. obs 자료)
        exceedence.prob.basic(x=x_, pfn=pfn_, cdfn=cdfn_, z=z)
        
      }else{ # x_ 의 dim이 null 이면 하위 리스트가 있음. (ex. model scenario 자료)
        lapply(x_, function(model.x_){
          exceedence.prob.basic(x=model.x_, pfn=pfn_, cdfn=cdfn_, z=z)
        })  
      }
    })
  names(res) <- vars
  res
}
