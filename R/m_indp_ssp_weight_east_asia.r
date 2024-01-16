#' indpendence weight 
#'
#' 기본적으로 RFA 를 적용하지 않은 자료로부터 계산함.
#'
#' by Yire Shin
#'
#' @param histDir 
#' @param sspDir 
#' @param year 
#' @param candi.sigma.s 
#'
#' @return
#' @export
#'
#' @examples
indp_ssp_east_asia <- function(histDir, sceDir, ref.period, sce.period, variableNames, candi.sigma.s, filter.grid=NULL){
  
  require(dplyr)
  
  ########## DEBUG ##########
  # histDir <- "/home/rstudio/workspace/nas/stat8039/git_projects/nimr2020_rproj/pi_weight/data/2020.10.07/CMIP6_EastAsia_variables_nbc2/historical"
  # s245Dir <- "/home/rstudio/workspace/nas/stat8039/git_projects/nimr2020_rproj/pi_weight/data/2020.10.07/CMIP6_EastAsia_variables_nbc2/ssp245"
  # s370Dir <- "/home/rstudio/workspace/nas/stat8039/git_projects/nimr2020_rproj/pi_weight/data/2020.10.07/CMIP6_EastAsia_variables_nbc2/ssp370"
  # s585Dir <- "/home/rstudio/workspace/nas/stat8039/git_projects/nimr2020_rproj/pi_weight/data/2020.10.07/CMIP6_EastAsia_variables_nbc2/ssp585"
  # ref.period <- 1850:2010
  # sce.period <- 2015:2100
  # sceDir <- c(s245Dir, s370Dir, s585Dir)
  # filter.grid <- NULL
  # variableNames <- c("AMP1", "AMP5", "ATP", "AMCWD", "AMCDD")
  # candi.sigma.s = seq(0.1, 1, 0.05)
  ########## DEBUG ##########
  
  cat("Read Datas..\n")
  ## Read Ref-m.datas
  hist.tmp.tb <- 
    readata.fromDirPath(histDir, ref.period, variable_names = variableNames, filter.grid = filter.grid, rfa9 = F) %>%
    tibble::enframe(name = "variable", value = "data") %>%
    dplyr::mutate(data=purrr::map(data, tibble::enframe)) %>% 
    tidyr::unnest(col=data)
  
  ## add 'sce' col to Ref-m.datas
  hist.tmp.tb <- tibble::add_column(tibble::tibble(sce=rep(0, nrow(hist.tmp.tb))), hist.tmp.tb)
  
  ## Read Sce-m.datas
  sce.tmp.tb <- 
    ## 미래 시나리오 디렉토리별
    lapply(sceDir, function(sce){
      tibble::enframe(readata.fromDirPath(sce, sce.period, variable_names = variableNames, filter.grid = filter.grid, rfa9 = F), 
                      name = "variable", value="data")
    }) %>% tibble::enframe(name="sce",value = "data") %>% tidyr::unnest(col=data) %>% 
    dplyr::mutate(data=purrr::map(data, tibble::enframe)) %>% 
    tidyr::unnest(col=data)
  
  ## row bind
  data.tb <- tibble::add_row(hist.tmp.tb, sce.tmp.tb)
  
  model.list <- data.tb$name %>% unique
  
  ## 자료 변환 ##
  # 지점별 7년 이동평균 수행 
  # moving average function
  
  ma <- function(x, n = 7){
    tmp <- tibble::as_tibble(stats::na.omit(stats::filter(x, rep(1 / n, n), sides = 2)))
    colnames(tmp) <- paste0("G",1:ncol(tmp))
    tmp
  }
  
  cat('Moving Average..')
  data.ma.tb <- data.tb %>% dplyr::mutate(ma_data=purrr::map(value, ma))
  cat('OK\n')
  
  ## 자료 변환 ##
  # 이동평균된 자료 합치기(hist + 시나리오 자료)_rbind 형태로 합치기
  # timeseries(hist+ssp245+ssp370+ssp585)
  
  cat('Binding Scenario Datas..')
  data.tb.tmp <- data.ma.tb %>% dplyr::group_by(variable, name) %>% tidyr::nest() %>% 
    dplyr::mutate(data=purrr::map(data, function(x){
      x %>% dplyr::select(-value) %>% # 필요 없는 컬럼 제외
        dplyr::arrange() %>% # sce 를 기준으로 오름차순 정렬 (0: ref 데이터, 1부터는 sce 데이터 차례대로임)
        tidyr::unnest(ma_data, names_repair="minimal") %>% # 각각 tibble 로 wrapping 되어있는 것을 unnest 처리
        dplyr::select(-sce)
    })) %>% dplyr::mutate(data = purrr::map(data, function(x){
      x %>% dplyr::mutate(t=1:n())
    }))
  cat('OK\n')
  
  nm <-length(model.list)
  
  #' loading distance 계산함수 ##
  #'
  #' by Yire Shin.
  #' @param pca 
  #' @param tnum 
  #'
  #' @return
  #' @export
  #'
  #' @examples
  calc_delta <-function(pca,tnum){
    delta <-matrix(rep(0,nm*nm),ncol=nm)#
    for(i in 1:nm){
      for(j in 1:nm){
        delta[i,j]<-dist(pca$rotation[c(i,j),1:tnum]) # pca$rotation # the matrix of variable loadings
      }        
    }
    return(delta)
  }
  
  cat('Shaping Datas..')
  ### 데이터 모양 변경
  data.stn.tb <- 
    data.tb.tmp %>% 
    tidyr::unnest(data) %>% 
    tidyr::gather(key="stn", value="value", -variable, -name, -t) %>% 
    tidyr::spread(name, value) %>% 
    dplyr::mutate(stn2=as.numeric(gsub(pattern = "G",replacement = "", x = stn))) %>%
    dplyr::select(-stn) %>%
    dplyr::arrange(variable, stn2, t) %>% 
    dplyr::group_by(variable, stn2) %>% 
    tidyr::nest()
  cat('Ok\n')
  
  cat('CORR, PCA, Delta.IJ..')
  ## 상관분석 수행 ##
  # prcomp()
  # 공분산행렬의 eigenvalue을 이용하지 않고, 
  # 원 데이터에 대해 SVD(Singular Value Decomposition : 특이값분해)를 
  # 수행하여 계산 
  # spearman 수행  
  data.stn.tb.tmp <- 
    data.stn.tb %>% 
    dplyr::mutate(corr=purrr::map(data, function(x){
      cor(x %>% dplyr::select(-t), method="spearman") 
    })) %>% 
    dplyr::mutate(pca=purrr::map(corr, function(x){ ## 주성분분석 수행 ##
      stats::prcomp(x)
    })) %>% 
    dplyr::mutate(tnum=purrr::map(corr, function(x){
      which.min(base::eigen(x)$value > 1) # Eigen value가 1보다 큰 주성분 선택 
    })) %>% 
    dplyr::mutate(delta_ij=purrr::map2(pca, tnum, calc_delta)) 
  cat('OK..')
  
  # 모델 i와 j간의 거리 모든 지역 평균 #
  d_mean.tb <- 
    data.stn.tb.tmp %>% 
    dplyr::group_by(variable) %>% 
    tidyr::nest() %>% 
    dplyr::mutate(data=purrr::map(data, ~abind::abind(.$delta_ij, along = 3))) %>% 
    dplyr::mutate(data=purrr::map(data, ~apply(., 1:2, mean)))
  
  d_mean_total.tmp <- apply(abind::abind(d_mean.tb$data, along = 3), 1:2, mean)
  
  
  #' 시간 단축을 위해 이레코드 중 가장 수정하기 쉬운 일부를 약간 변경함.
  do__ <- function(sigma.s){
    # independence 수식 적용
    indep1 <- 
      data.stn.tb.tmp %>% 
      dplyr::mutate(v_sta=purrr::map(delta_ij, ~apply(., 1, function(x) sum(exp(-abs(x^2)/sigma.s^2))))) %>% 
      dplyr::mutate(s_i=purrr::map(v_sta, ~1/.))
    
    indep2 <- 
      indep1 %>% 
      dplyr::group_by(stn2) %>% 
      tidyr::nest() %>% 
      dplyr::mutate(s_i_total=purrr::map(data, ~apply(do.call(rbind, .$s_i), 2, mean)))
    
    s_i_total.stn <- do.call(rbind, indep2$s_i_total)
    s_i.mean <- colMeans(s_i_total.stn) # mean of each station 1/sigma  (for all station)
    
    # independence 가중치 
    indep3 <- 
      indep1 %>% 
      dplyr::mutate(indp_i=purrr::map(v_sta, ~(1/.)/(sum(1/.)))) %>% 
      dplyr::arrange(variable, stn2)  # 다시 한번 정렬
    
    # mean of each station Indep Weight (for all station)
    indp_i.mean <- 
      indep3 %>% 
      dplyr::ungroup() %>% 
      dplyr::select(stn2, variable, indp_i) %>% 
      dplyr::group_by(stn2) %>% 
      tidyr::nest() %>% 
      dplyr::mutate(value=purrr::map(data, ~apply(do.call(rbind, .$indp_i), 2, mean))) %>% 
      dplyr::pull(value) %>% 
      do.call(what = rbind) %>% 
      colMeans
    
    # model.list <- data.tb$name %>% unique
    
    return(list(indp_i_mean=indp_i.mean,model.list=model.list))
  }
  
  cat('Calc Idp-weights..')
  res_ <- pbapply::pblapply(candi.sigma.s, do__)
  names(res_) <- candi.sigma.s
  cat('Ok\n')
  return(res_)
}


#' Entropy 를 이용하여 sigma.S 선택
#'
#' @param model.his.dirpath 
#' @param model.sce.dirpath 
#' @param year 
#' @param candi.sigma.s 
#' @param mc.cores mc.cores를 2 이상으로 주면, parallel 로 실행, 리눅스에서만 실행할것.
#'
#' @return
#' @export
#'
#' @examples
calc.sigma.s <- function(hisDir, sceDir, ref.period, sce.period, variableNames, candi.sigma.s, filter.grid=NULL){#, mc.cores=1){
  
  ###### parallel 기능 없애면서 주석 처리
  # indp_ <-
  # if(mc.cores>1){
  #   if((Sys.info())["sysname"] == "Windows") stop("Sorry, Parallel execution on 'Windows systems' is not implemented. Set 'mc.cores=1' and run again")
  #   usethis::use_package("parallel")
  #   usethis::use_package("pbmcapply")
  #   if(mc.cores>length(candi.sigma.s)){
  #     cat("mc.cores =",mc.cores, ",length(candi.sigma.s) =", length(candi.sigma.s),"\n")
  #     print("Since mc.cores is larger than candy.sigma.s, the value of mc.cores has been automatically modified.")
  #     mc.cores <- length(candi.sigma.s)
  #   }
  #   
  #   pbmcapply::pbmclapply(candi.sigma.s, function(ss){
  #     indp_ssp_east_asia(model.his.dirpath, model.sce.dirpath, year, ss)
  #   }, mc.cores=mc.cores)
  #   
  # }else{
  #   lapply(candi.sigma.s, function(ss){
  #     indp_ssp_east_asia(model.his.dirpath, model.sce.dirpath, year, ss)
  #   })
  # }
  # indp_ <- indp_ssp_east_asia_org(hisDir, sceDir, year=c(1850,2010, 2015, 2100), candi.sigma.s) # by Yire Shin.
  
  indp_ <- indp_ssp_east_asia(histDir, sceDir, ref.period, sce.period, variableNames, candi.sigma.s, filter.grid)
  
  tmp.sigma.s <- 
    lapply(indp_, function(indp__){
      E.pi <- indp__$indp_i_mean
      -sum(E.pi * log(E.pi))
    }) %>% do.call(what = rbind) %>% data.frame(candi.sigma.s)  
  
  colnames(tmp.sigma.s) <- c("value", "sigma.s")
  optim.sigma.s <- tmp.sigma.s %>% filter(value==min(value))
  
  min.idx <- which.min(tmp.sigma.s$value)
  optim.res <- indp_[[min.idx]]
  
  # Weights
  indep.weight.df <- data.frame(model=optim.res$model.list, w=optim.res$indp_i_mean) %>% dplyr::arrange(model)
  
  return(list(res=indp_, optim.res=optim.res, tmp.sigma.s=tmp.sigma.s, sigma.S=optim.sigma.s, weights=indep.weight.df))
}








