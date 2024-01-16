
if(FALSE){
  # Setting -----------------------------------------------------------------
  library(piw)
  
  setwd("~/workspace/projects/pi-weights-rpackage/")
  
  obs.file.path <- 
    "not_importing_data/2020.10.05/OBS.rds"
  
  # model 자료 경로
  histDir.bc <- "not_importing_data/2020.10.05/historical"
  s245Dir.bc <- "not_importing_data/2020.10.05/ssp245"
  s370Dir.bc <- "not_importing_data/2020.10.05/ssp370"
  s585Dir.bc <- "not_importing_data/2020.10.05/ssp585"
  
  histDir <- "not_importing_data/2020.10.07/CMIP6_EastAsia_variables_nbc2/historical"
  s245Dir <- "not_importing_data/2020.10.07/CMIP6_EastAsia_variables_nbc2/ssp245"
  s370Dir <- "not_importing_data/2020.10.07/CMIP6_EastAsia_variables_nbc2/ssp370"
  s585Dir <- "not_importing_data/2020.10.07/CMIP6_EastAsia_variables_nbc2/ssp585"
  
  
  # 그림 저장 경로
  eps.save.path <- "not_importing_output/eps/"
  pdf.save.path <- "not_importing_output/pdf/"
  png.save.path <- "not_importing_output/png/"
  
  
  # Set period : RFA
  filter.period <- 1973:2010
  
  # set Return Level to be considered
  rt.lv <- c(2, 5, 10, 20, 30, 50, 100)
  
  vars <- c("AMP1","AMP5","ATP","AMCWD","AMCDD")
  
  # Load Datas --------------------------------------------------------------
  
  obs <- readata.fromFilePath(obs.file.path, filter.period, variable_names = vars, rfa9 = T)
  his    <- readata.fromDirPath(histDir   , filter.period, variable_names = vars, rfa9 = T)
  his.bc <- readata.fromDirPath(histDir.bc, filter.period, variable_names = vars, rfa9 = T)
  
  # Calc Return Lv. ---------------------------------------------------------
  
  library(lmom)
  
  # Parameter estimation Functions by L-moments
  pfns <- list(pelgev, pelgev, pelnor, pelgev, pelgev)
  
  # Quantile Functions
  qfns <- list(quagev, quagev, quanor, quagev, quagev)
  
  obs.rtlv <- rtlv.multi.var(obs, pfns, qfns)
  his.rtlv <- rtlv.multi.var(his, pfns, qfns)
  his.rtlv.bc <- rtlv.multi.var(his, pfns, qfns)
  
  
  # Normalization -----------------------------------------------------------
  
  #'  기존의 연구에서는 rt.lv <- c(2, 5, 10, 20, 30, 50, 100) 을 이용하여 정규화 수행.
  
  obs.tilde <- tilde.normalization(obs, his, target.rtlv = obs.rtlv)
  his.tilde <- tilde.normalization(obs, his, target.rtlv = his.rtlv)
  
  
  # Performance Weight ------------------------------------------------------
  
  #'
  #'  기존 연구에서는 sigma.D 를 동아시아 전체 자료를 이용해서 결정함.
  #'  어떤 자료를 이용하는냐에 따라 sigma.D 의 값이 달라질것임.
  #'  예를 들면 한반도 자료만을 이용해서 sigma.D 를 계산한것과 
  #'  동아시아 전체 자료를 이용해서 sigma.D를 계산한 것은 다를 수 있음.
  
  
  ## Sigma.D
  res.sigma.d <- calc.sigma.d(ref = obs.tilde, 
                              model = his.tilde, 
                              simul.N = 10000, 
                              random.seed = 197079,
                              significance.level = 0.05,
                              candi.sigma.d = seq(1.0, 0.01, -0.01))
  
  sigma.D <- res.sigma.d$sigma.D
  
  ## 모델별 Perf-weights 값은 지역별로 다름.
  
  ## Weights
  perf.weight.res <- calc.performance.weights(ref=obs.tilde, model = his.tilde, sigma.D)
  perf.weight <- perf.weight.res$weight.df
  
  # Independence Weight -----------------------------------------------------
  
  ## 기존 연구에서는 sigma.S 를 동아시아 전체 자료를 이용해서 결정함.
  ## 어떤 자료를 이용하는냐에 따라 sigma.S 의 값이 달라질것임.
  ## 예를 들면 한반도 자료만을 이용해서 sigma.S 를 계산한것과 
  ## 동아시아 전체 자료를 이용해서 sigma.S를 계산한 것은 다를 수 있음.
  
  # Find Optimal Sigma.S from Entropy Approach
  system.time({
    res.sigma.s <- 
      calc.sigma.s(hisDir = histDir, 
                   sceDir = c(s245Dir,s370Dir,s585Dir),
                   ref.period = 1850:2010, 
                   sce.period = 2015:2100,
                   variableNames = vars,
                   candi.sigma.s = seq(0.1, 1, 0.05),
                   filter.grid = NULL)
  })
  # user  system elapsed 
  # 184.231  39.811 195.221 
  
  # user  system elapsed 
  # 106.124  18.840 106.374
  
  sigma.s <- res.sigma.s$sigma.S$sigma.s
  

  # >> Main : Figure 2 ---------------------------------------------------------
  # Visualization
  plot.sigma.s.res(res.sigma.s)
  
  ## perf-weights 와는 다르게 모델별 ind-weights 는 전 지역에 대해 동일한 값임.
  ## Weights
  indep.weight <- res.sigma.s$weights
  
  
  # PI Weight ---------------------------------------------------------------
  
  ## PI-weight 를 계산하기 위해 사용된 자료는 모두 bias-correction을 하지 않은 자료임.
  ## bias-correction을 적용한 자료를 사용할 경우, 모든 모델이 등가중치를 가지게 되어 
  
  pi.weight <- calc.pi.weight(perf.weight_ = perf.weight, indep.weight_ = indep.weight)
  

  # >> Main : Figure 5 ------------------------------------------------------
  ## plot
  plot.weights(perf.weight, indep.weight, pi.weight) +
    ggplot2::ggtitle(as.expression(bquote(sigma[S] ~ "=" ~ .(sprintf("%.2f",sigma.s)) ~ "," ~ sigma[D] ~ "=" ~ .(sprintf("%.2f",sigma.D )))))
  
  
  #### 여기까지 pi-weight 계산 끝
  
  #### 이제부터 한반도 격자에 대해서 pi-weight 앙상블 적용

  # Ensemble RtLv & RtPd ----------------------------------------------------
  
  kr.grids <- c(393:394,443:444,494:497,552:556,616:620,672:676,
                729:732,786:791,844:849,900:904,950:951)
  
  sce.ensembled.bc <- 
    ensemble.sce(refDir = obs.file.path,
                 refPeriod = 1973:2010,
                 sceDir = c(s245Dir.bc, s370Dir.bc, s585Dir.bc), 
                 sceName = c("s245","s370","s585"), 
                 scePeriod = list(p1=2021:2050,p2=2046:2075,p3=2071:2100),
                 variableNames = c("AMP1"), 
                 pfn = list(lmom::pelgev),
                 qfn = list(lmom::quagev),
                 weight = pi.weight,
                 rt.lv = c(2, 5, 10, 20, 30, 50, 100, 150, 200),
                 ref.filter.grid = kr.grids,
                 sce.filter.grid = kr.grids,
                 weight.filter.grid = kr.grids,
                 rfa9 = TRUE
    )
  
  
  
  # Return Level boxplot ----------------------------------------------------
  
  # 한반도 지역의 obs : return level 계산
  obs <- readata.fromFilePath(obs.file.path, 1973:2010, "AMP1", kr.grids, rfa9=TRUE)
  obs.rtlv <- rtlv.multi.var(obs, list(lmom::pelgev), list(lmom::quagev))
  
  # 한반도 지역 his : return level 계산
  his <- readata.fromDirPath(histDir, 1973:2010, "AMP1", kr.grids, rfa9 = TRUE)
  his.rtlv <- rtlv.multi.var(his, list(lmom::pelgev), list(lmom::quagev))
  
  # 상자그림에 보여줄 his : return level의 평균 계산
  his.ensembled <- 
    lapply(his.rtlv, function(x){
      weighted.sum(rtlv.list = x, weights = "equal")
    })
  
  # >> Main : Figure 7 ------------------------------------
  # 20y return level 상자그림 그리기
  rtlv.boxplot(ref.rtlv = obs.rtlv, 
               his.rtlv.ensembled = his.ensembled,
               sce.rtlv.ensembled = sce.ensembled.bc$return.level,
               variable.name = "AMP1", 
               rt.name = "rt20",
               ylim=c(50, 500),
               label=c("SSP2-4.5","SSP3-7.0","SSP5-8.5"),
               label.y.coord = 500,
               label.x.coord = c(4,7,10),
               axis.x.label=c("OBS","HIS(NBC)","P1","P2","P3","P1","P2","P3","P1","P2","P3"),
               y.label = "Return Level", 
               calibrate.text.size=-3)
  

  # >> SM : Figure S 3 ------------------------------------
  # 50y return level 상자그림 그리기
  rtlv.boxplot(ref.rtlv = obs.rtlv, 
               his.rtlv.ensembled = his.ensembled,
               sce.rtlv.ensembled = sce.ensembled.bc$return.level,
               variable.name = "AMP1", 
               rt.name = "rt50",
               ylim=c(50, 500),
               label=c("SSP2-4.5","SSP3-7.0","SSP5-8.5"),
               label.y.coord = 500,
               label.x.coord = c(4,7,10),
               axis.x.label=c("OBS","HIS(NBC)","P1","P2","P3","P1","P2","P3","P1","P2","P3"),
               y.label = "Return Level(50y)", 
               calibrate.text.size=-3)
  

  # Return Period Boxplot ----------------------------------------------------

  rtpd.bplot <- 
    rtpd.boxplot(ensemble.sce.obj = sce.ensembled.bc,
                 var.name = "AMP1",
                 rtlv = c("rt20","rt50"), 
                 facet.label = c("rt20"="20-year","rt50"="50-year"),
                 axis.x.names = c("P1","P2","P3","P1","P2","P3","P1","P2","P3"),
                 calibrate.text.size=-1)
  
  # 라벨 추가
  lb <- data.frame(x=c(2,5,8,2,5,8), 
                   y=c(rep(65,3), rep(65,3)),
                   label=c("SSP2-4.5","SSP3-7.0","SSP5-8.5","SSP2-4.5","SSP3-7.0","SSP5-8.5"),
                   variable = c(rep("rt20",3), rep("rt50", 3)))

  # >> Main : Figure 9 ------------------------------------------------------
  rtpd.bplot +
    ggplot2::scale_y_continuous(limits = c(3, 65), breaks = seq(0, 65, by = 10))+
    ggplot2::geom_label(data=lb, ggplot2::aes(x=x, y=y, label=label, group=variable), inherit.aes = F, size=7/2.2)
  
  
  

  # Variance based on BMA ---------------------------------------------------
  
  # 분산 계산, 논문에서는 n.boots = 500 으로 했음.
  var.rtlv.pi <- var.rtlv(sce.ensembled.bc, n.boots = 500, random.seed = 121212)
  

  # >> Main : Figure 14 ------------------------------------------------------
  library(ggpubr)
  my.col <- RColorBrewer::brewer.pal(3, 'Set1')
  my.col <- c(my.col[1], 'black', my.col[3])
  ggpubr::ggarrange(
    var.scatter.plot(var.rtlv.pi, "AMP1", "rt20", "total")+ ggplot2::scale_color_manual(values = my.col),
    var.scatter.plot(var.rtlv.pi, "AMP1", "rt20", "vbm")  + ggplot2::scale_color_manual(values = my.col),
    var.scatter.plot(var.rtlv.pi, "AMP1", "rt20", "vwm")  + ggplot2::scale_color_manual(values = my.col),
    nrow = 3, common.legend = TRUE, legend='top'
  )
  

  # >> Main : Table 4 -------------------------------------------------------
  bss.table(var.rtlv.obj = var.rtlv.pi, var.name="AMP1", rtlv = 'rt20')
  
  # Exceedence Probability --------------------------------------------------
  kr.grids <- c(393:394,443:444,494:497,552:556,616:620,672:676,
                729:732,786:791,844:849,900:904,950:951)
  
  z_ = c(seq(50, 500, 50), 750, 1000)
  
  # 한국 격자에 대한 가중값만 뽑아서 list 형태로 변환
  pi.weight.list <-
    split(x = pi.weight %>% dplyr::select(-model), f = pi.weight$model) %>%
    lapply(function(x){
      x %>% dplyr::select(kr.grids) %>% t()
    })
  
  # AMP1 자료만 사용, RFA 적용하지 않은 채로, 한국 격자에 대해서만
  obs <- readata.fromFilePath(obs.file.path, 
                              filter.period = 1973:2010, 
                              variable_names = 'AMP1', 
                              filter.grid = kr.grids, 
                              rfa9 = F)
  # 확률 계산
  obs.exc.prob <- 
    exceedence.prob.multi.var(data.obj = obs, 
                              pfns = list(lmom::pelgev),
                              cdfns = list(lmom::cdfgev), 
                              z = z_)
  
  # AMP1 자료만 사용, RFA 적용하지 않은 채로, 한국 격자에 대해서만
  his.bc <- readata.fromDirPath(histDir.bc, 
                                filter.period = 1973:2010, 
                                variable_names = 'AMP1', 
                                filter.grid = kr.grids, 
                                rfa9 = F)
  # historical 자료 확률 계산하고 앙상블 
  his.bc.exc.prob.ensembled <-  
    exceedence.prob.multi.var(data.obj = his.bc, 
                              pfns = list(lmom::pelgev),
                              cdfns = list(lmom::cdfgev), z = z_) %>% 
    lapply(function(p.var.model){
      weighted.sum(p.var.model, pi.weight.list)
    })
  
  # scenario 자료,  한국 격자에 대해서만, RFA 하지 않음.
  # 자료 읽고 바로 확률 계산해서 앙상블 적용 
  sce.bc.exc.prob.ensembled <- 
    ## 미래 시나리오 디렉토리별
    lapply(c(s245Dir.bc, s370Dir.bc, s585Dir.bc), function(sce){
      ## 미래 시나리오 기간별
      lapply(list(p1=2021:2050,p2=2046:2075,p3=2071:2100), function(p){
        
        # 데이터 읽어서
        data <- readata.fromDirPath(sce, p, variable_names = 'AMP1', filter.grid = kr.grids, rfa9 = F)
        
        # exc prob 계산하고
        exc.prob <- 
          exceedence.prob.multi.var(data.obj = data, 
                                    pfns = list(lmom::pelgev),
                                    cdfns = list(lmom::cdfgev),
                                    z = z_)
        # 앙상블 
        lapply(exc.prob, function(p.var.model){
          weighted.sum(p.var.model, pi.weight.list)
        })
      })
    }) %>% `names<-`(c('s245','s370','s585'))
  
  # exceedence probability 결과 모두 합치기
  exc.prob.all <- 
    rbind(
      # obs 그림을 p1, p2, p3 에 똑같이 넣기 위해..
      reshape2::melt(obs.exc.prob) %>% 
        `colnames<-`(c('loc','z','value','var')) %>% 
        dplyr::mutate(period='p1', sce='obs'),
      reshape2::melt(obs.exc.prob) %>% 
        `colnames<-`(c('loc','z','value','var')) %>% 
        dplyr::mutate(period='p2', sce='obs'),
      reshape2::melt(obs.exc.prob) %>% 
        `colnames<-`(c('loc','z','value','var')) %>% 
        dplyr::mutate(period='p3', sce='obs'),
      
      # his 결과
      reshape2::melt(his.bc.exc.prob.ensembled) %>% 
        `colnames<-`(c('loc','z','value','var')) %>% 
        dplyr::mutate(period='p0', sce='hist'),
      
      # sce 결과 
      reshape2::melt(sce.bc.exc.prob.ensembled) %>% 
        `colnames<-`(c('loc','z','value','var','period','sce'))
    ) %>% dplyr::mutate(z.n = as.numeric(gsub(pattern = 'z', replacement = '', x = z)))
  
  # 모든 격자에 대해 median 값 계산
  exc.prob.all.med <- exc.prob.all %>% 
    dplyr::select(-loc, -z) %>% 
    dplyr::group_by(var, sce, period, z.n) %>% 
    tidyr::nest() %>% 
    dplyr::mutate(med=purrr::map(data, ~median(.$value))) %>% 
    tidyr::unnest(med)
 

  # >> SM : Figure S 6 ------------------------------------
  exc.prob.all.med %>% 
    dplyr::filter(var == "AMP1", sce %in% c("obs", "s245", "s370", "s585"), z.n <= 300) %>%
    ggplot2::ggplot(mapping = ggplot2::aes(x = z.n, y = med, colour=sce)) +
    ggplot2::geom_line() +
    ggplot2::scale_x_continuous(breaks = z_[z_<=300])+
    ggplot2::facet_grid(~period) +
    ggplot2::theme_bw() + 
    ggplot2::ylab("Probability") + 
    ggplot2::xlab("AMP1")+
    ggplot2::theme(legend.position = "top", 
                   legend.title    = ggplot2::element_blank())
  

  # >> SM : Table S3 -----------------------------------------------------------
  # Return Level stats
  rtlv.table(return.level.obj = sce.ensembled.bc$return.level, 
                  rt.lv = c('rt20', 'rt50'))
  
  
  
}