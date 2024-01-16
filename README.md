# PI-Weights-Rpackage

Future Projections and Uncertainty Assessment of Precipitation Extremes in the Korean Peninsula from the CMIP6 Ensemble with a Statistical Framework.
[https://www.mdpi.com/2073-4433/12/1/97]

## Installation

Install the latest development version (on GitHub) via `{remotes}`:

``` r
remotes::install_github("sygstat/piw")
```

## Getting started

#### Load Datasets

If you would like to download the dataset used in this example, please contact us via email.
(jspark@chonnam.ac.kr or syg.stat@gmail.com)

```r
library(piw)

setwd("~/workspace/projects/pi-weights-rpackage/")

obs.file.path <- 
  "not_importing_data/2020.10.05/OBS.rds"

# model data path
histDir.bc <- "not_importing_data/2020.10.05/historical"
s245Dir.bc <- "not_importing_data/2020.10.05/ssp245"
s370Dir.bc <- "not_importing_data/2020.10.05/ssp370"
s585Dir.bc <- "not_importing_data/2020.10.05/ssp585"

histDir <- "not_importing_data/2020.10.07/CMIP6_EastAsia_variables_nbc2/historical"
s245Dir <- "not_importing_data/2020.10.07/CMIP6_EastAsia_variables_nbc2/ssp245"
s370Dir <- "not_importing_data/2020.10.07/CMIP6_EastAsia_variables_nbc2/ssp370"
s585Dir <- "not_importing_data/2020.10.07/CMIP6_EastAsia_variables_nbc2/ssp585"
  
  
# save path for figures
eps.save.path <- "not_importing_output/eps/"
pdf.save.path <- "not_importing_output/pdf/"
png.save.path <- "not_importing_output/png/"


# Set period : RFA
filter.period <- 1973:2010

# set Return Level to be considered
rt.lv <- c(2, 5, 10, 20, 30, 50, 100)

vars <- c("AMP1","AMP5","ATP","AMCWD","AMCDD")

obs <- readata.fromFilePath(obs.file.path, filter.period, variable_names = vars, rfa9 = T)
his    <- readata.fromDirPath(histDir   , filter.period, variable_names = vars, rfa9 = T)
his.bc <- readata.fromDirPath(histDir.bc, filter.period, variable_names = vars, rfa9 = T)
```

#### Calc rtlv
```r
library(lmom)

# Parameter estimation Functions by L-moments
pfns <- list(pelgev, pelgev, pelnor, pelgev, pelgev)

# Quantile Functions
qfns <- list(quagev, quagev, quanor, quagev, quagev)

obs.rtlv <- rtlv.multi.var(obs, pfns, qfns)
his.rtlv <- rtlv.multi.var(his, pfns, qfns)
his.rtlv.bc <- rtlv.multi.var(his, pfns, qfns)
```

#### Normalization
```r

obs.tilde <- tilde.normalization(obs, his, target.rtlv = obs.rtlv)
his.tilde <- tilde.normalization(obs, his, target.rtlv = his.rtlv)
```

#### Get the Performance Weights
```r

## Sigma.D
res.sigma.d <- calc.sigma.d(ref = obs.tilde, 
                            model = his.tilde, 
                            simul.N = 10000, 
                            random.seed = 197079,
                            significance.level = 0.05,
                            candi.sigma.d = seq(1.0, 0.01, -0.01))

sigma.D <- res.sigma.d$sigma.D

## Weights
perf.weight.res <- calc.performance.weights(ref=obs.tilde, model = his.tilde, sigma.D)
perf.weight <- perf.weight.res$weight.df
```

#### Get the Independence Weights
```r
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
```

#### Figure 2
```r
# Visualization
plot.sigma.s.res(res.sigma.s)
```

#### Figure 5
```r

## Weights
indep.weight <- res.sigma.s$weights

pi.weight <- calc.pi.weight(perf.weight_ = perf.weight, indep.weight_ = indep.weight)

## plot
plot.weights(perf.weight, indep.weight, pi.weight) +
  ggplot2::ggtitle(as.expression(bquote(sigma[S] ~ "=" ~ .(sprintf("%.2f",sigma.s)) ~ "," ~ sigma[D] ~ "=" ~ .(sprintf("%.2f",sigma.D )))))
```

#### Ensemble rtlv & rtpd
```r
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
```

#### Boxplot
```r
# obs : return level 
obs <- readata.fromFilePath(obs.file.path, 1973:2010, "AMP1", kr.grids, rfa9=TRUE)
obs.rtlv <- rtlv.multi.var(obs, list(lmom::pelgev), list(lmom::quagev))

# his : return level 
his <- readata.fromDirPath(histDir, 1973:2010, "AMP1", kr.grids, rfa9 = TRUE)
his.rtlv <- rtlv.multi.var(his, list(lmom::pelgev), list(lmom::quagev))

# his : mean of return level for boxplot
his.ensembled <- 
  lapply(his.rtlv, function(x){
    weighted.sum(rtlv.list = x, weights = "equal")
  })

# >> Main : Figure 7 ------------------------------------
# 20y return level boxplot
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
# 50y return level boxplot
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
```


  
  
```r
# Variance based on BMA ---------------------------------------------------

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

pi.weight.list <-
  split(x = pi.weight %>% dplyr::select(-model), f = pi.weight$model) %>%
  lapply(function(x){
    x %>% dplyr::select(kr.grids) %>% t()
  })

obs <- readata.fromFilePath(obs.file.path, 
                            filter.period = 1973:2010, 
                            variable_names = 'AMP1', 
                            filter.grid = kr.grids, 
                            rfa9 = F)
# clac prob
obs.exc.prob <- 
  exceedence.prob.multi.var(data.obj = obs, 
                            pfns = list(lmom::pelgev),
                            cdfns = list(lmom::cdfgev), 
                            z = z_)


his.bc <- readata.fromDirPath(histDir.bc, 
                              filter.period = 1973:2010, 
                              variable_names = 'AMP1', 
                              filter.grid = kr.grids, 
                              rfa9 = F)

his.bc.exc.prob.ensembled <-  
  exceedence.prob.multi.var(data.obj = his.bc, 
                            pfns = list(lmom::pelgev),
                            cdfns = list(lmom::cdfgev), z = z_) %>% 
  lapply(function(p.var.model){
    weighted.sum(p.var.model, pi.weight.list)
  })


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

# exceedence probability 
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

# median
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
  
```
## Learn more

...

## Citation

...

-----
