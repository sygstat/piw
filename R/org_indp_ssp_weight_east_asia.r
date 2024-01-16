#' indpendence weight 
#'
#' 기본적으로 RFA 를 적용하지 않은 자료로부터 계산함.
#'
#' by Yire Shin
#'
indp_ssp_east_asia_org <- function(histDir,sspDir,year,candi.sigma.s){
  
  
  ########## DEBUG ##########
  histDir <- "/home/rstudio/workspace/nas/stat8039/git_projects/nimr2020_rproj/pi_weight/data/2020.10.07/CMIP6_EastAsia_variables_nbc2/historical"
  s245Dir <- "/home/rstudio/workspace/nas/stat8039/git_projects/nimr2020_rproj/pi_weight/data/2020.10.07/CMIP6_EastAsia_variables_nbc2/ssp245"
  s370Dir <- "/home/rstudio/workspace/nas/stat8039/git_projects/nimr2020_rproj/pi_weight/data/2020.10.07/CMIP6_EastAsia_variables_nbc2/ssp370"
  s585Dir <- "/home/rstudio/workspace/nas/stat8039/git_projects/nimr2020_rproj/pi_weight/data/2020.10.07/CMIP6_EastAsia_variables_nbc2/ssp585"
  year <- c(1850, 2010, 2015, 2100)
  sspDir <- c(s245Dir, s370Dir, s585Dir)
  filter.grid <- NULL
  variableNames <- c("AMP1", "AMP5", "ATP", "AMCWD", "AMCDD")
  candi.sigma.s = seq(0.1, 1, 0.05)
  ########## DEBUG ##########
  
  lf_hist <-list.files(histDir)
  lf_s245 <-list.files(sspDir[1])
  lf_s370 <-list.files(sspDir[2])
  lf_s585 <-list.files(sspDir[3])
  
  # ssp 데이터 목록 합치기 #
  
  lf_ssp  <-c(lf_s245,lf_s370,lf_s585)
  
  # hist 모델 목록중 ssp 모델 목록만 선택 #
  
  model_selc <-lf_hist[lf_hist %in% unique(lf_ssp)]
  
  # 선택된 모델 목록 # 
  
  lf_hist_selc <-which(lf_hist %in% model_selc)
  lf_s245_selc <-which(lf_s245 %in% model_selc)
  lf_s370_selc <-which(lf_s370 %in% model_selc)
  lf_s585_selc <-which(lf_s585 %in% model_selc)
  
  # 선택된 모델 이름 #
  
  model_selc <-model_selc[lf_hist_selc]
  
  ## 모델 선택2 : 사용자가 정의한 reference 기간       ##
  ##              & 시나리오 기간이 일치하는 모델 선택 ##
  
  ref_year  <-c(year[1]:year[2])
  ssp_year  <-c(year[3]:year[4])
  
  
  ## 기간이 다른 모델 선택 #
  
  
  # 최종 입력 모델 선택 #
  
  model_selc <-model_selc[union(union(union(lf_hist_selc,lf_s245_selc), lf_s370_selc),lf_s585_selc)]
  
  model.list <- gsub(x=model_selc, pattern = ".rds", replacement = "")
  
  nm <-length(model.list)
  
  
  hist <- list()
  s245 <- list()
  s370 <- list()
  s585 <- list()
  
  system.time({
    for(i in 1:nm){
      
      hist[[i]]<-readRDS(file.path(histDir,lf_hist[lf_hist_selc[i]]))
      s245[[i]]<-readRDS(file.path(s245Dir,lf_s245[lf_s245_selc[i]]))
      s370[[i]]<-readRDS(file.path(s370Dir,lf_s370[lf_s370_selc[i]]))
      s585[[i]]<-readRDS(file.path(s585Dir,lf_s585[lf_s585_selc[i]]))
      
      #print(lf_hist[lf_hist_selc[i]])  
    }
    
    
    ## 자료 변환1 ##
    
    # raw data form information #
    # row = latitude(37), column = longitude(78), array = timeseries(161)
    # list_num = model(22), real value = 41, number of NA = 99-41 =58
    
    
    # transform data  #
    # 모델별 AMP1, AMP5, ATP, AMCWD, AMcdd
    
    # row = total station(2886), column = time series(161) # 
    # list number = model(22)
    
    
    var_name <- c("AMP1", "AMP5", "ATP", "AMCWD", "AMDWD")
    
    
    hist_amp1  <- lapply(hist, function(x) apply(x$AMP1, 3,function(x) c(x)))
    hist_amp5  <- lapply(hist, function(x) apply(x$AMP5, 3,function(x) c(x)))  
    hist_atp   <- lapply(hist, function(x) apply(x$ATP,  3,function(x) c(x))) 
    hist_amcwd <- lapply(hist, function(x) apply(x$AMCWD,3,function(x) c(x)))   
    hist_amcdd <- lapply(hist, function(x) apply(x$AMCDD,3,function(x) c(x)))  
    
    s245_amp1  <- lapply(s245, function(x) apply(x$AMP1, 3,function(x) c(x)))
    s245_amp5  <- lapply(s245, function(x) apply(x$AMP5, 3,function(x) c(x)))  
    s245_atp   <- lapply(s245, function(x) apply(x$ATP,  3,function(x) c(x))) 
    s245_amcwd <- lapply(s245, function(x) apply(x$AMCWD,3,function(x) c(x)))   
    s245_amcdd <- lapply(s245, function(x) apply(x$AMCDD,3,function(x) c(x)))  
    
    s370_amp1  <- lapply(s370, function(x) apply(x$AMP1, 3,function(x) c(x)))
    s370_amp5  <- lapply(s370, function(x) apply(x$AMP5, 3,function(x) c(x)))  
    s370_atp   <- lapply(s370, function(x) apply(x$ATP,  3,function(x) c(x))) 
    s370_amcwd <- lapply(s370, function(x) apply(x$AMCWD,3,function(x) c(x)))   
    s370_amcdd <- lapply(s370, function(x) apply(x$AMCDD,3,function(x) c(x)))  
    
    s585_amp1  <- lapply(s585, function(x) apply(x$AMP1, 3,function(x) c(x)))
    s585_amp5  <- lapply(s585, function(x) apply(x$AMP5, 3,function(x) c(x)))  
    s585_atp   <- lapply(s585, function(x) apply(x$ATP,  3,function(x) c(x))) 
    s585_amcwd <- lapply(s585, function(x) apply(x$AMCWD,3,function(x) c(x)))   
    s585_amcdd <- lapply(s585, function(x) apply(x$AMCDD,3,function(x) c(x)))  
    
  })
  
  # number of stations area which have value #
  
  tot_num  <-seq(1:length(hist_amp1[[1]][,1])) # total station number = 1,...2886
  sta_num  <-which(!is.na(hist_amp1[[1]][,1])) # observed station number = 3,4,12,..98
  sta_len  <-length(sta_num)                  # number of station =1180
  na_num   <-which(!tot_num %in% sta_num)     # non observed station number/ length(na_num) = 1706
  
  
  ## 자료 변환2 ##
  
  # 지점별 모델별 강수자료
  # row = time series, col = model number, list number = station number#
  
  hist_amp1_stn <-lapply(1:sta_len, function(s) sapply(hist_amp1, function(x,num=s) x[sta_num[num],]))
  hist_amp5_stn <-lapply(1:sta_len, function(s) sapply(hist_amp5, function(x,num=s) x[sta_num[num],]))
  hist_atp_stn <-lapply(1:sta_len, function(s) sapply(hist_atp, function(x,num=s) x[sta_num[num],]))
  hist_amcwd_stn <-lapply(1:sta_len, function(s) sapply(hist_amcwd, function(x,num=s) x[sta_num[num],]))
  hist_amcdd_stn <-lapply(1:sta_len, function(s) sapply(hist_amcdd, function(x,num=s) x[sta_num[num],]))
  
  s245_amp1_stn <-lapply(1:sta_len, function(s) sapply(s245_amp1, function(x,num=s) x[sta_num[num],]))
  s245_amp5_stn <-lapply(1:sta_len, function(s) sapply(s245_amp5, function(x,num=s) x[sta_num[num],]))
  s245_atp_stn <-lapply(1:sta_len, function(s) sapply(s245_atp, function(x,num=s) x[sta_num[num],]))
  s245_amcwd_stn <-lapply(1:sta_len, function(s) sapply(s245_amcwd, function(x,num=s) x[sta_num[num],]))
  s245_amcdd_stn <-lapply(1:sta_len, function(s) sapply(s245_amcdd, function(x,num=s) x[sta_num[num],]))
  
  s370_amp1_stn <-lapply(1:sta_len, function(s) sapply(s370_amp1, function(x,num=s) x[sta_num[num],]))
  s370_amp5_stn <-lapply(1:sta_len, function(s) sapply(s370_amp5, function(x,num=s) x[sta_num[num],]))
  s370_atp_stn <-lapply(1:sta_len, function(s) sapply(s370_atp, function(x,num=s) x[sta_num[num],]))
  s370_amcwd_stn <-lapply(1:sta_len, function(s) sapply(s370_amcwd, function(x,num=s) x[sta_num[num],]))
  s370_amcdd_stn <-lapply(1:sta_len, function(s) sapply(s370_amcdd, function(x,num=s) x[sta_num[num],]))
  
  s585_amp1_stn <-lapply(1:sta_len, function(s) sapply(s585_amp1, function(x,num=s) x[sta_num[num],]))
  s585_amp5_stn <-lapply(1:sta_len, function(s) sapply(s585_amp5, function(x,num=s) x[sta_num[num],]))
  s585_atp_stn <-lapply(1:sta_len, function(s) sapply(s585_atp, function(x,num=s) x[sta_num[num],]))
  s585_amcwd_stn <-lapply(1:sta_len, function(s) sapply(s585_amcwd, function(x,num=s) x[sta_num[num],]))
  s585_amcdd_stn <-lapply(1:sta_len, function(s) sapply(s585_amcdd, function(x,num=s) x[sta_num[num],]))
  
  
  
  
  ## 자료 변환3 ##
  # 지점별 7년 이동평균 수행 
  # moving average function
  # row = time series, col = model number, list number = station number
  
  ma <- function(x, n = 7){stats::filter(x, rep(1 / n, n), sides = 2)}
  
  hist_ma_amp1_stn<-lapply(hist_amp1_stn, function(x) na.omit(ma(x)))
  hist_ma_amp5_stn<-lapply(hist_amp5_stn, function(x) na.omit(ma(x)))
  hist_ma_atp_stn<-lapply(hist_atp_stn, function(x) na.omit(ma(x)))
  hist_ma_amcwd_stn<-lapply(hist_amcwd_stn, function(x) na.omit(ma(x)))
  hist_ma_amcdd_stn<-lapply(hist_amcdd_stn, function(x) na.omit(ma(x)))
  
  s245_ma_amp1_stn<-lapply(s245_amp1_stn, function(x) na.omit(ma(x)))
  s245_ma_amp5_stn<-lapply(s245_amp5_stn, function(x) na.omit(ma(x)))
  s245_ma_atp_stn<-lapply(s245_atp_stn, function(x) na.omit(ma(x)))
  s245_ma_amcwd_stn<-lapply(s245_amcwd_stn, function(x) na.omit(ma(x)))
  s245_ma_amcdd_stn<-lapply(s245_amcdd_stn, function(x) na.omit(ma(x)))
  
  s370_ma_amp1_stn<-lapply(s370_amp1_stn, function(x) na.omit(ma(x)))
  s370_ma_amp5_stn<-lapply(s370_amp5_stn, function(x) na.omit(ma(x)))
  s370_ma_atp_stn<-lapply(s370_atp_stn, function(x) na.omit(ma(x)))
  s370_ma_amcwd_stn<-lapply(s370_amcwd_stn, function(x) na.omit(ma(x)))
  s370_ma_amcdd_stn<-lapply(s370_amcdd_stn, function(x) na.omit(ma(x)))
  
  s585_ma_amp1_stn<-lapply(s585_amp1_stn, function(x) na.omit(ma(x)))
  s585_ma_amp5_stn<-lapply(s585_amp5_stn, function(x) na.omit(ma(x)))
  s585_ma_atp_stn<-lapply(s585_atp_stn, function(x) na.omit(ma(x)))
  s585_ma_amcwd_stn<-lapply(s585_amcwd_stn, function(x) na.omit(ma(x)))
  s585_ma_amcdd_stn<-lapply(s585_amcdd_stn, function(x) na.omit(ma(x)))
  
  
  
  ## 자료 변환4 ##
  # 이동평균된 자료 합치기(hist + 시나리오 자료)_rbind 형태로 합치기
  # row = timeseries(hist+ssp245+ssp370+ssp585), col = model number, list number = station number
  
  data_amp1_stn <-lapply(1:sta_len, 
                         function(x) rbind(hist_ma_amp1_stn[[x]],s245_ma_amp1_stn[[x]],
                                           s370_ma_amp1_stn[[x]],s585_ma_amp1_stn[[x]]))
  
  data_amp5_stn <-lapply(1:sta_len, 
                         function(x) rbind(hist_ma_amp5_stn[[x]],s245_ma_amp5_stn[[x]],
                                           s370_ma_amp5_stn[[x]],s585_ma_amp5_stn[[x]]))
  
  data_atp_stn <-lapply(1:sta_len, 
                        function(x) rbind(hist_ma_atp_stn[[x]],s245_ma_atp_stn[[x]],
                                          s370_ma_atp_stn[[x]],s585_ma_atp_stn[[x]]))
  
  data_amcwd_stn <-lapply(1:sta_len, 
                          function(x) rbind(hist_ma_amcwd_stn[[x]],s245_ma_amcwd_stn[[x]],
                                            s370_ma_amcwd_stn[[x]],s585_ma_amcwd_stn[[x]]))
  
  data_amcdd_stn <-lapply(1:sta_len, 
                          function(x) rbind(hist_ma_amcdd_stn[[x]],s245_ma_amcdd_stn[[x]],
                                            s370_ma_amcdd_stn[[x]],s585_ma_amcdd_stn[[x]]))
  
  
  
  ## 상관분석 수행 ##
  # prcomp()
  # 공분산행렬의 eigenvalue을 이용하지 않고, 
  # 원 데이터에 대해 SVD(Singular Value Decomposition : 특이값분해)를 
  # 수행하여 계산 
  # spearman 수행  
  
  corr_amp1_stn <-lapply(data_amp1_stn, function(x) cor(x,method="spearman"))
  corr_amp5_stn <-lapply(data_amp5_stn, function(x) cor(x,method="spearman"))
  corr_atp_stn <-lapply(data_atp_stn, function(x) cor(x,method="spearman"))
  corr_amcwd_stn <-lapply(data_amcwd_stn, function(x) cor(x,method="spearman"))
  corr_amcdd_stn <-lapply(data_amcdd_stn, function(x) cor(x,method="spearman"))
  
  
  
  ## 주성분분석 수행 ##
  # pca$rotation # the matrix of variable loadings
  # Eigen value가 1보다 큰 주성분 선택 
  
  pca_amp1   <-lapply(corr_amp1_stn,   function(x) prcomp(x)) # list number # 41
  pca_amp5   <-lapply(corr_amp5_stn,   function(x) prcomp(x)) # list number # 41
  pca_atp    <-lapply(corr_atp_stn,    function(x) prcomp(x)) # list number # 41
  pca_amcwd  <-lapply(corr_amcwd_stn,  function(x) prcomp(x)) # list number # 41
  pca_amcdd  <-lapply(corr_amcdd_stn,  function(x) prcomp(x)) # list number # 41
  
  
  tnum_amp1  <-sapply(1:sta_len, function(x) which.min(eigen(corr_amp1_stn[[x]])$value > 1))
  tnum_amp5  <-sapply(1:sta_len, function(x) which.min(eigen(corr_amp5_stn[[x]])$value > 1))
  tnum_atp   <-sapply(1:sta_len, function(x) which.min(eigen(corr_atp_stn[[x]])$value > 1))
  tnum_amcwd <-sapply(1:sta_len, function(x) which.min(eigen(corr_amcwd_stn[[x]])$value > 1))
  tnum_amcdd <-sapply(1:sta_len, function(x) which.min(eigen(corr_amcdd_stn[[x]])$value > 1))
  
  
  #tnum <-sapply(pca, function(x) which.max(summary(x)$importance[3,] > 0.8))
  
  
  
  ## loading distance 계산함수 ##
  
  calc_delta <-function(pca,tnum){
    
    delta <-matrix(rep(0,nm*nm),ncol=nm)#
    
    for(i in 1:nm){
      
      for(j in 1:nm){
        
        delta[i,j]<-dist(pca$rotation[c(i,j),1:tnum])
        
      }        
      
    }
    
    return(delta)
    
  }
  
  
  ## loading distance 계산 ## 
  
  delta_ij_amp1  <-list()
  delta_ij_amp5  <-list()
  delta_ij_atp   <-list()
  delta_ij_amcwd <-list()
  delta_ij_amcdd <-list()
  
  for(i in 1:sta_len){
    
    delta_ij_amp1[[i]]  <-calc_delta(pca_amp1[[i]],tnum=tnum_amp1[i])
    delta_ij_amp5[[i]]  <-calc_delta(pca_amp5[[i]],tnum=tnum_amp5[i])
    delta_ij_atp[[i]]   <-calc_delta(pca_atp[[i]],tnum=tnum_atp[i])
    delta_ij_amcwd[[i]] <-calc_delta(pca_amcwd[[i]],tnum=tnum_amcwd[i])
    delta_ij_amcdd[[i]] <-calc_delta(pca_amcdd[[i]],tnum=tnum_amcdd[i])
    
  }
  
  
  
  ## 함수 출력값 ##
  
  
  # 모델 i와 j간의 loading distance
  # 가까우면 dependence, 멀면 independence
  # row = model i, col = model j, list number = station number
  
  d_dist_amp1    <-delta_ij_amp1
  d_dist_amp5    <-delta_ij_amp5
  d_dist_atp     <-delta_ij_atp
  d_dist_amcwd   <-delta_ij_amcwd
  d_dist_amcdd   <-delta_ij_amcdd
  
  
  
  # d_dist list에서 array 변환 <- 평균구하기 위해서 
  d_array_amp1   <-array(as.numeric(unlist(d_dist_amp1)), dim=c(nm, nm, sta_len))
  d_array_amp5   <-array(as.numeric(unlist(d_dist_amp5)), dim=c(nm, nm, sta_len))
  d_array_atp    <-array(as.numeric(unlist(d_dist_atp)), dim=c(nm, nm, sta_len))
  d_array_amcwd  <-array(as.numeric(unlist(d_dist_amcwd)), dim=c(nm, nm, sta_len))
  d_array_amcdd  <-array(as.numeric(unlist(d_dist_amcdd)), dim=c(nm, nm, sta_len))
  
  
  
  # 모델 i와 j간의 거리 모든 지역 평균 #
  d_mean_amp1    <-apply(d_array_amp1,  c(1,2), function(x) mean(x))
  d_mean_amp5    <-apply(d_array_amp5,  c(1,2), function(x) mean(x))
  d_mean_atp     <-apply(d_array_atp,   c(1,2), function(x) mean(x))
  d_mean_amcwd   <-apply(d_array_amcwd, c(1,2), function(x) mean(x))
  d_mean_amcdd   <-apply(d_array_amcdd, c(1,2), function(x) mean(x))
  d_mean_total   <-(d_mean_amp1+d_mean_amp5+d_mean_atp+d_mean_amcwd+d_mean_amcdd)/5
  
  
  #' 시간 단축을 위해 이레코드 중 가장 수정하기 쉬운 일부를 약간 변경함.
  do__ <- function(sigma.s){
    # independence 수식 적용  
    v_sta_amp1     <-lapply(delta_ij_amp1, function(x) apply(x, 1, function(x) sum(exp(-abs(x^2)/sigma.s^2))))
    v_sta_amp5     <-lapply(delta_ij_amp5, function(x) apply(x, 1, function(x) sum(exp(-abs(x^2)/sigma.s^2))))
    v_sta_atp      <-lapply(delta_ij_atp, function(x) apply(x, 1, function(x) sum(exp(-abs(x^2)/sigma.s^2))))
    v_sta_amcwd    <-lapply(delta_ij_amcwd, function(x) apply(x, 1, function(x) sum(exp(-abs(x^2)/sigma.s^2))))
    v_sta_amcdd    <-lapply(delta_ij_amcdd, function(x) apply(x, 1, function(x) sum(exp(-abs(x^2)/sigma.s^2))))
    
    
    s_i_amp1    <-lapply(v_sta_amp1, function(x) (1/x))
    s_i_amp5    <-lapply(v_sta_amp5, function(x) (1/x))
    s_i_atp     <-lapply(v_sta_atp, function(x) (1/x))
    s_i_amcwd   <-lapply(v_sta_amcwd, function(x) (1/x))
    s_i_amcdd   <-lapply(v_sta_amcdd, function(x) (1/x))
    
    
    s_i_total    <-list()
    
    # mean of each 5 variable Indep Weight (for each station)
    
    for(i in 1:length(s_i_amp1)){
      
      s_i_total[[i]] <-(s_i_amp1[[i]]+s_i_amp5[[i]]+s_i_atp[[i]]+s_i_amcwd[[i]]+s_i_amcdd[[i]])/5
      
    }
    
    s_i_total_stn <- matrix(as.numeric(unlist(s_i_total)),ncol=nm, byrow = T)
    
    
    # mean of each station 1/sigma  (for all station)
    
    s_i_mean <- colMeans(s_i_total_stn)
    
    
    # independence 가중치 
    indp_i_amp1    <-lapply(v_sta_amp1, function(x) (1/x)/(sum(1/x)))
    indp_i_amp5    <-lapply(v_sta_amp5, function(x) (1/x)/(sum(1/x)))
    indp_i_atp     <-lapply(v_sta_atp, function(x) (1/x)/(sum(1/x)))
    indp_i_amcwd   <-lapply(v_sta_amcwd, function(x) (1/x)/(sum(1/x)))
    indp_i_amcdd   <-lapply(v_sta_amcdd, function(x) (1/x)/(sum(1/x)))
    
    indp_i_total    <-list()
    
    # mean of each 5 variable Indep Weight (for each station)
    
    for(i in 1:length(indp_i_amp1)){
      
      indp_i_total[[i]] <-(indp_i_amp1[[i]]+indp_i_amp5[[i]]+indp_i_atp[[i]]+indp_i_amcwd[[i]]+indp_i_amcdd[[i]])/5
      
    }
    
    indp_i_total_stn <- matrix(as.numeric(unlist(indp_i_total)),ncol=nm, byrow = T)
    
    
    # mean of each station Indep Weight (for all station)
    
    indp_i_mean <- colMeans(indp_i_total_stn)
    
    
    return(list(d_dist_amp1   =d_dist_amp1,
                d_dist_amp5   =d_dist_amp5,
                d_dist_atp    =d_dist_atp,
                d_dist_amcwd  =d_dist_amcwd,
                d_dist_amcdd  =d_dist_amcdd,
                d_mean_amp1   =d_mean_amp1,        
                d_mean_amp5   =d_mean_amp5,    
                d_mean_atp    =d_mean_atp,    
                d_mean_amcwd  =d_mean_amcwd,      
                d_mean_amcdd  =d_mean_amcdd,
                d_mean_total  =d_mean_total,
                v_sta_amp1    =v_sta_amp1,          
                v_sta_amp5    =v_sta_amp5,          
                v_sta_atp     =v_sta_atp,      
                v_sta_amcwd   =v_sta_amcwd,        
                v_sta_amcdd   =v_sta_amcdd,        
                indp_i_amp1   =indp_i_amp1,
                indp_i_amp5   =indp_i_amp5,      
                indp_i_atp    =indp_i_atp,        
                indp_i_amcwd  =indp_i_amcwd,    
                indp_i_amcdd  =indp_i_amcdd,  
                indp_i_total  =indp_i_total,
                indp_i_mean   =indp_i_mean,
                s_i_mean    =s_i_mean,
                model.list=model.list))
  }
  
  res_ <- lapply(candi.sigma.s, do__)
  names(res_) <- candi.sigma.s
  return(res_)
}

