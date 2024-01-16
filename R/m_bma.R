##################################################
## 예전 BMA 모듈 모아놓은것
## by Yonggwan, 2021.05.13
##################################################

# 시나리오 자료에 bma 가중치 적용하여 리턴레벨 계산
#' Title
#'
#' @param obj 
#' @param weight 
#' @param rt.lv 
#'
#' @return
#' @export
#'
#' @examples
applied_scena_for_weight <- function(obj, weight, rt.lv=c(2, 5, 10, 20, 50)) {
  obj.list <- split(obj %>% dplyr::select(-model), obj$model)
  tmp.rtlv <- 
    lapply(X = obj.list, FUN = function(x){
      tmp.table <- 
        t(apply(X = x, MARGIN = 2, FUN = function(data){
          require(lmom)
          s.x <- data[!is.na(data)]
          quagev(1-1/rt.lv, pelgev(samlmu(s.x, 3)))
        }))
      colnames(tmp.table) <- paste0("rt",rt.lv)
      tmp.table
    })
  # print(tmp.rtlv)
  before.sum <- lapply(X = seq_along(tmp.rtlv), FUN = function(k, m.rtlv, weight){
    ncol_ = ncol(m.rtlv[[k]])
    #m.rtlv[[k]]*weight[[k]]
    # print('-----')
    # print(matrix(rep(weight[[k]], ncol_), ncol=ncol_))
    # print(apply(matrix(rep(weight[[k]], ncol_), ncol=ncol_), 2, sum))
    m.rtlv[[k]]*matrix(rep(weight[[k]], ncol_), ncol=ncol_)
  }, m.rtlv = tmp.rtlv, weight = weight)
  
  # summation
  apply(simplify2array(before.sum), 1:2, sum)
  
}


# rfa 가중합으로 바뀌면서 rtlv 를 파라메터로 받는걸로 바뀜
#' 입력 arg 두 리스트 모두 모델이름으로 네이밍이 되어 있어야함.
#'
#' @param rtlv.list 
#' @param weights 
#'
#' @return
#' @export
#'
#' @examples
weighted.sum <- function(rtlv.list, weights) {
  # usethis::use_package("abind")
  
  model.list <- names(rtlv.list)
  
  if( (is.character(weights)) && (weights=="equal")){
    weights <- lapply(rtlv.list, function(x) matrix(1/length(model.list), ncol=1, nrow=nrow(x)))
  }
  
  if(length(rtlv.list) != length(weights)) stop("Length of 'rtlv.list' and 'weights' should be the same.")
  if(!all(names(rtlv.list) == names(weights))) stop("Names of 'rtlv.list' and 'weights' should be the same.")
  
  before.sum <- lapply(X = model.list, FUN = function(model_){
    ncol_ <- ncol(rtlv.list[[model_]])
    rtlv.list[[model_]]*matrix(rep(weights[[model_]], ncol_), ncol=ncol_)
  })
  
  # summation
  apply(abind::abind(before.sum, along = 3), 1:2, sum)
  
}

#' Title
#'
#' @param x.rp 
#' @param weight.list 
#'
#' @return
#' @export
#'
#' @examples
return.period.weighted.sum <- function(x.rp, weight.list, rt.lv){
  
  rp.3d.array <- # 모델 자료 -> 3d-array
    abind::abind(
      lapply(x.rp, function(x) x %>% dplyr::select(-loc)),
      along = 3)
  
  mask.3d.array <- # 이상값 masking (3d-array)
    apply(rp.3d.array, 2:3, function(x){
      is.finite(x) & x<100 # Stisfy both 'finite value' and 'x<100'
    })
  
  weight.3d.array <- # weight.list -> 3d-array
    abind::abind(
      lapply(weight.list, function(pw){
        replicate(length(rt.lv), pw, simplify = T)
      }), along = 3)
  
  rp.rs <- matrix(NA, ncol = ncol(rp.3d.array), nrow=nrow(rp.3d.array))
  for(y in 1:ncol(rp.3d.array)){
    for(x in 1:nrow(rp.3d.array)){
      rp.values <- rp.3d.array[x,y,]
      m.values <- mask.3d.array[x,y,]
      w.values <- weight.3d.array[x,y,]
      
      rp.values[is.infinite(rp.values)] <- 0
      w.values <- w.values/sum(w.values[m.values]) # re-weight
      w.values[!m.values] <- 0
      rp.rs[x,y] <- sum(rp.values*w.values) # weigthed sum
    }
  }
  
  colnames(rp.rs) <- colnames(x.rp[[1]])[-1] # except for the name of first column : 'loc' 
  rownames(rp.rs) <- x.rp[[1]]$loc
  return(as.matrix(rp.rs))
}
