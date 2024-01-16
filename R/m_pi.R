#' Title
#'
#' @param sceDir 
#' @param sceName 
#' @param scePeriod 
#' @param variableNnames 
#' @param pfn 
#' @param qfn 
#' @param weight 
#' @param rt.lv 
#' @param filter.gird 
#' @param rfa9 
#'
#' @return
#' @export
#'
#' @examples
ensemble.sce <- 
  function(refDir, refPeriod, sceDir, sceName, scePeriod, variableNames, pfn, qfn, weight, 
           rt.lv=c(2, 5, 10, 20, 30, 50), 
           ref.filter.grid=NULL, sce.filter.grid=NULL, 
           weight.filter.grid=NULL,
           rfa9=FALSE){
    
    require(dplyr)
    
    names(sceDir) <- sceName
    names(pfn) <- variableNames
    names(qfn) <- variableNames
    
    ref <- readata.fromFilePath(refDir, refPeriod, variable_names = variableNames, filter.grid = ref.filter.grid, rfa9 = rfa9)
    ref.rtlv <- rtlv.multi.var(ref, pfns = pfn, qfns = qfn, rt.lv = rt.lv)
    
    ref.rtlv <- 
      lapply(ref.rtlv, function(var.data){
        ref.rtlv_ <- as.data.frame(var.data)
        ref.rtlv_$loc <- rownames(ref.rtlv_)
        ref.rtlv_ %>% dplyr::arrange(loc)
      })
    
    ## TODO :: 시나리오별 모두 모델 개수가 같은지, 모델 순서가 같은지 확인하는 코드 삽입
    ## TODO :: 그러기 위해서는 데이터를 읽는 것과 rtlv 구하는 것을 분리해서 작성해야 함. 나중에 하는 걸로..
    ## 반드시 모든 시나리오 디렉토리의 모델 개수와 모델 목록이 동일해야 함. 
    ## 만약 다를 경우 올바른 결과를 보장하지 못함.
    cat('Load Datas..\n')    
    sce.data <- 
      ## 미래 시나리오 디렉토리별
      lapply(sceDir, function(sce){
        ## 미래 시나리오 기간별
        lapply(scePeriod, function(p){
          readata.fromDirPath(sce, p, variable_names = variableNames, filter.grid = sce.filter.grid, rfa9 = rfa9)
        })
      })
    
    
    cat('Calc ReturnLevel & ReturnPeriod..\n')
    sce.stats <- 
      ## 미래 시나리오 디렉토리별
      lapply(sce.data, function(scePeriod){ 
        ## 미래 시나리오 기간별
        lapply(scePeriod, function(var.data){
          ## 변수별
          var.res <- 
            lapply(variableNames, function(var.name){
              ## 모델별      
              ### return level
              rtlv <- 
                lapply(var.data[[var.name]], function(data_){
                  rtlv.basic(data_, pfn[[var.name]], qfn[[var.name]], rt.lv = rt.lv) 
                })  
              
              ### return period
              rtpd <- 
                lapply(var.data[[var.name]], function(data_){
                  return.period.basic(x = data_, pfn = pfn[[var.name]], ref.rtlv = ref.rtlv[[var.name]], rt.lv = rt.lv)
                })  
              list(return.level = rtlv, return.period=rtpd)        
            })
          names(var.res) <- variableNames
          var.res
        })
      })
    
    weight.list <-
      split(x = weight %>% dplyr::select(-model), f = weight$model) %>%
      lapply(function(x){
        if(is.null(weight.filter.grid)){
          x %>% t() 
        }else{
          x %>% dplyr::select(weight.filter.grid) %>% t()
        }
      })
    
    cat('Do Ensemble..')
    
    sce.rtlv.ensembled <- 
      ## 미래 시나리오별
      lapply(sce.stats, function(sce.p){  ## grid 적용 안되는 듯.
        ## 미래 시나리오 기간별
        lapply(sce.p, function(sce.p.var){
          ## 변수별
          lapply(sce.p.var, function(sce.p.var.model){
            weighted.sum(sce.p.var.model$return.level, weight.list)
          })
        })
      })
    
    sce.rtpd.ensembled <- 
      ## 미래 시나리오별
      lapply(sce.stats, function(sce.p){  ## grid 적용 안되는 듯.
        ## 미래 시나리오 기간별
        lapply(sce.p, function(sce.p.var){
          ## 변수별
          lapply(sce.p.var, function(sce.p.var.model){
            return.period.weighted.sum(sce.p.var.model$return.period, weight.list, rt.lv)
          })
        })
      })
    
    cat('OK\n')
    list(sce.data      = sce.data,
         sce.stats     = sce.stats, # 앙상블 하기 전 자료
         return.level  = sce.rtlv.ensembled,
         return.period = sce.rtpd.ensembled,
         weight.list   = weight.list)
    
  }
#' #' Title
#' #'
#' #' @param sceDir 
#' #' @param sceName 
#' #' @param scePeriod 
#' #' @param variableNnames 
#' #' @param pfn 
#' #' @param qfn 
#' #' @param weight 
#' #' @param rt.lv 
#' #' @param filter.gird 
#' #' @param rfa9 
#' #'
#' #' @return
#' #' @export
#' #'
#' #' @examples
#' ensemble.sce <- 
#'   function(sceDir, sceName, scePeriod, variableNames, pfn, qfn, weight, rt.lv=c(2, 5, 10, 20, 30, 50), filter.grid=NULL, rfa9=FALSE){
#'     usethis::use_package("dplyr")
#'     require(dplyr)
#'     
#'     names(sceDir) <- sceName
#'     
#'     ## TODO :: 시나리오별 모두 모델 개수가 같은지, 모델 순서가 같은지 확인하는 코드 삽입
#'     ## TODO :: 그러기 위해서는 데이터를 읽는 것과 rtlv 구하는 것을 분리해서 작성해야 함. 나중에 하는 걸로..
#'     ## 반드시 모든 시나리오 디렉토리의 모델 개수와 모델 목록이 동일해야 함. 
#'     ## 만약 다를 경우 올바른 결과를 보장하지 못함.
#'     sce.rtlv <- 
#'       ## 미래 시나리오 디렉토리별
#'       lapply(sceDir, function(sce){  ## grid 적용 안되는 듯.
#'         ## 미래 시나리오 기간별
#'         lapply(scePeriod, function(p){
#'           var.data <- readata.fromDirPath(sce, p, variable_names = variableNames, filter.grid = filter.grid, rfa9 = rfa9)
#'           ## 변수별
#'           lapply(var.data, function(model.data){
#'             ## 모델별      
#'             lapply(model.data, function(data_){
#'               rtlv.basic(data_, pfn, qfn, rt.lv = rt.lv)      
#'             })  
#'           })
#'           
#'         })
#'       })
#'     
#'     weight.list <-
#'       split(x = weight %>% dplyr::select(-model), f = weight$model) %>%
#'       lapply(function(x){
#'         if(is.null(filter.grid)){
#'           x %>% t() 
#'         }else{
#'           x %>% dplyr::select(filter.grid) %>% t()
#'         }
#'       })
#'     
#'     sce.rtlv.ensembled <- 
#'       ## 미래 시나리오별
#'       lapply(sce.rtlv, function(sce.p){  ## grid 적용 안되는 듯.
#'         ## 미래 시나리오 기간별
#'         lapply(sce.p, function(sce.p.var){
#'           ## 변수별
#'           lapply(sce.p.var, function(sce.p.var.model){
#'             weighted.sum(sce.p.var.model, weight.list)
#'           })
#'         })
#'       })
#'   }
#'   

#


#' Title
#'
#' @param sce.data 
#' @param sce.stats 
#' @param ensembled.rtlv 
#' @param weight.list 
#'
#' @return
#' @export
#'
#' @examples
var.rtlv <- function(ensemble.sce.obj, n.boots=100, random.seed=121212){
  # var.rtlv <- function(sce.data, sce.stats, ensembled.rtlv, weight.list){
  require(dplyr)
  set.seed(random.seed)
  
  # ensemble.sce.obj <- sce.ensembled.bc
  
  sce.data = ensemble.sce.obj$sce.data
  sce.stats = ensemble.sce.obj$sce.stats
  ensembled.rtlv = ensemble.sce.obj$return.level
  weight.list = ensemble.sce.obj$weight.list
  
  eq.weight.list <- 
    lapply(weight.list, function(x){
      x[] <- 1/length(weight.list)
      x
    })
  
  model.value.df <- 
    reshape2::melt(sce.stats, id.vars="loc") %>% 
    dplyr::filter(L4=="return.level") %>% 
    dplyr::select(-loc, -variable, -L4) %>% 
    dplyr::group_by(Var2, Var1, L1, L2, L3) %>% tidyr::nest()
  colnames(model.value.df) <- c("loc","rt","variable","period","sce","model.value")
  
  ensembled.rtlv.df <- 
    reshape2::melt(ensembled.rtlv, id.vars="loc") %>%
    dplyr::group_by(Var2, Var1, L1, L2, L3)
  colnames(ensembled.rtlv.df) <- c("loc","rt","ensembled.value","variable","period","sce")
  
  weight.df <- 
    reshape2::melt(weight.list) %>% 
    dplyr::select(-L1) %>% 
    dplyr::group_by(Var1) %>% tidyr::nest()
  colnames(weight.df) <- c("loc","model.weight")
  
  eq.weight.df <- 
    reshape2::melt(eq.weight.list) %>% 
    dplyr::select(-L1) %>% 
    dplyr::group_by(Var1) %>% tidyr::nest()
  colnames(eq.weight.df) <- c("loc","model.eq.weight")
  
  joined_.df <- dplyr::left_join(model.value.df, ensembled.rtlv.df, by=c("variable"="variable","loc"="loc", "rt"="rt","period"="period","sce"="sce"))
  joined_.df <- dplyr::left_join(joined_.df, weight.df, by=c("loc"="loc"))
  full.df    <- dplyr::left_join(joined_.df, eq.weight.df, by=c("loc"="loc"))
  
  calc.Vbm <- function(xk, xbar, w){
    Xk = xk %>% dplyr::arrange(L5)
    W = w %>% dplyr::arrange(Var2)
    
    tibble::tibble(model=W$Var2, value=(Xk$value - xbar)^2 * W$value)
  }
  
  full.df <-   
    full.df %>% 
    dplyr::mutate(simple.mean.value=purrr::map_dbl(model.value, ~mean(.$value))) # simple mean
  
  cat("Calc Vbm..")
  res_ <- 
    full.df %>% 
    dplyr::mutate(Vbm=purrr::pmap(.l = list(model.value, ensembled.value, model.weight), .f = calc.Vbm)) %>%
    dplyr::mutate(Vbm.eq=purrr::pmap(.l = list(model.value, simple.mean.value, model.eq.weight), .f = calc.Vbm))
  
  res_ <- res_ %>% 
    dplyr::mutate(Vbm.sum=purrr::map_dbl(Vbm, function(x) sum(x$value))) %>%
    dplyr::mutate(Vbm.eq.sum=purrr::map_dbl(Vbm.eq, function(x) sum(x$value)))
  cat("Ok\n")
  
  #### Vwm
  #### Vbm에서 계산된 rt.lv 목록 가져옴(Vwm 에서 쓰려고)
  rt.lv_ <- as.numeric(gsub(x = unique(res_$rt), pattern = "rt", replacement = ""))
  
  sce.data <- reshape2::melt(sce.data) %>% 
    dplyr::select(-Var1)
  colnames(sce.data) <- c("loc","value","model","variable","period","sce")
  
  sce.data <- 
    sce.data %>% 
    dplyr::group_by(variable, loc, sce, period, model) %>% 
    tidyr::nest() %>%
    dplyr::arrange(variable, loc, sce, period)
  
  boots.gev <- function(data, n.boots=100){
    s.x <- data[!is.na(data)]
    n.gen.data <- length(s.x)
    gev.param <- lmom::pelgev(lmom::samlmu(s.x, 3))
    
    boots.res <- 
      sapply(1:n.boots, function(bi){
        rgev.data <- 
          SpatialExtremes::rgev(n.gen.data, loc=gev.param[1], scale = gev.param[2], shape = gev.param[3]*-1)
        
        lmom::quagev(1-1/rt.lv_, lmom::pelgev(lmom::samlmu(rgev.data, 3)))
        
      })
    rownames(boots.res) <- paste0('rt',rt.lv_)
    reshape2::melt(t(apply(boots.res, MARGIN = 1, FUN = stats::var))) %>% 
      dplyr::select(-Var1)
  }
  
  ## Vwm 계산
  cat("Bootstrap for var[rk(T)]..(it would take a few minutes)..")
  sce.boots.data <- sce.data %>% dplyr::mutate(bdata.variance=purrr::map(data, boots.gev, n.boots=n.boots))
  sce.boots.data <- sce.boots.data %>% tidyr::unnest(col=bdata.variance)
  cat("Ok\n")
  
  cat("Calc Vwm..")
  tmp__ <- 
    sce.boots.data %>%
    dplyr::left_join(weight.df %>% tidyr::unnest(col=model.weight),
                     by=c("loc"="loc", "model"="Var2")) %>%
    dplyr::left_join(eq.weight.df %>% tidyr::unnest(col=model.eq.weight),
                     by=c("loc"="loc", "model"="Var2")) %>%
    dplyr::mutate(Vwm=purrr::map2_dbl(value.x, value.y, function(var.value,weight.value){
      var.value*weight.value
    })) %>%
    dplyr::mutate(Vwm.eq=purrr::map2_dbl(value.x, value, function(var.value,weight.value){
      var.value*weight.value
    })) %>%
    dplyr::group_by(variable, Var2, loc, sce, period) %>%
    tidyr::nest() %>%
    dplyr::mutate(Vwm.sum=purrr::map_dbl(data, ~sum(.$Vwm))) %>%
    dplyr::mutate(Vwm.eq.sum=purrr::map_dbl(data, ~sum(.$Vwm.eq)))
  cat("Ok\n")
  
  cat("Final..")
  result <- 
    dplyr::left_join(res_, tmp__, 
                     by = c("variable"="variable","loc"="loc","rt"="Var2","period"="period","sce"="sce")) %>%
    dplyr::mutate(total.var=purrr::map2_dbl(Vbm.sum, Vwm.sum, function(vbm_, vwm_) vbm_+vwm_)) %>%
    dplyr::mutate(total.var.eq=purrr::map2_dbl(Vbm.eq.sum, Vwm.eq.sum, function(vbm_, vwm_) vbm_+vwm_))
  cat("Ok\n")
  
  result
}


BSS <- function(ref.value, new.value){(1-(new.value/ref.value))*100}

DV <- function(ref.value, new.value){ref.value - new.value}

RI <- function(ref.value, new.value){100*(ref.value-new.value)/new.value}

bss.table <- function(var.rtlv.obj, var.name, rtlv){
  require(dplyr)
  var.rtlv.obj %>% 
    dplyr::filter(variable==var.name) %>% 
    dplyr::mutate(bss=purrr::map2_dbl(total.var.eq, total.var, BSS)) %>%
    dplyr::mutate(bss.bm=purrr::map2_dbl(Vbm.eq.sum, Vbm.sum, BSS)) %>% 
    dplyr::mutate(bss.wm=purrr::map2_dbl(Vwm.eq.sum, Vwm.sum, BSS)) %>%
    dplyr::mutate(dv=purrr::map2_dbl(total.var.eq, total.var, DV)) %>%
    dplyr::mutate(dv.bm=purrr::map2_dbl(Vbm.eq.sum, Vbm.sum, DV)) %>%
    dplyr::mutate(dv.wm=purrr::map2_dbl(Vwm.eq.sum, Vwm.sum, DV)) %>%
    dplyr::mutate(ri=purrr::map2_dbl(total.var.eq, total.var, RI)) %>%
    dplyr::select(bss:ri) %>% 
    dplyr::ungroup() %>% 
    dplyr::group_by(rt, variable, period, sce) %>% 
    tidyr::nest() %>% 
    dplyr::mutate(bss.mean    =purrr::map_dbl(data, ~mean(.$bss)),
                  bss.bm.mean =purrr::map_dbl(data, ~mean(.$bss.bm)),
                  bss.wm.mean =purrr::map_dbl(data, ~mean(.$bss.wm)),
                  dv.mean     =purrr::map_dbl(data, ~round(mean(.$dv))),
                  dv.bm.mean  =purrr::map_dbl(data, ~round(mean(.$dv.bm))),
                  dv.wm.mean  =purrr::map_dbl(data, ~round(mean(.$dv.wm))),
                  ri.mean     =purrr::map_dbl(data, ~round(mean(.$ri))),
    ) %>% 
    dplyr::filter(rt==rtlv) %>% 
    dplyr::select(-data)
}


rtlv.table <- function(return.level.obj, rt.lv){
  # return.level.obj <- sce.ensembled.bc$return.level
  # rt.lv <- c('rt20', 'rt50')
  require(dplyr)
  tmp.res <- 
    reshape2::melt(return.level.obj) %>% 
    `colnames<-`(c('loc','rt','value','var','period','sce')) %>% 
    dplyr::group_by(var, sce, period, rt) %>% 
    tidyr::nest() %>% 
    dplyr::mutate(mean   = purrr::map_dbl(data, ~round(mean(.$value))),
                  q1     = purrr::map_dbl(data, ~round(stats::quantile(.$value, 0.25))),
                  median = purrr::map_dbl(data, ~round(mean(.$value))),
                  q3     = purrr::map_dbl(data, ~round(stats::quantile(.$value, 0.75))),
                  sce.p  = paste0(sce,period)) %>% 
    dplyr::filter(rt %in% rt.lv) %>% 
    dplyr::arrange(rt) %>% 
    dplyr::ungroup() %>% 
    dplyr::select(-data, -period, -sce) %>% 
    tidyr::gather(stats, value, -rt, -var, -sce.p) %>% 
    tidyr::spread(sce.p, value)
  tmp.res$stats <- factor(tmp.res$stats, levels = c('mean', 'q1', 'median', 'q3'))
  
  tmp.res %>% dplyr::arrange(var, rt, stats)
}
