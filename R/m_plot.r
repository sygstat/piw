##################################################
## plot 관련 모듈들
## by Yonggwan, 2021.05.13
##################################################

#' Title
#'
#' @param perf.weight_ 
#' @param indep.weight_ 
#' @param pi.weight_ 
#' @param arrange.by.pi 
#'
#' @return
#' @export
#'
#' @examples
plot.weights <- function(perf.weight_, indep.weight_, pi.weight_, arrange.by.pi=TRUE){
  require(dplyr)
  # usethis::use_package("reshape2")
  # usethis::use_package("ggplot2")
  # usethis::use_package("dplyr")
  
  melted.total.weight <- 
    reshape2::melt(
      list(
        reshape2::melt(perf.weight_) %>% dplyr::group_by(model) %>% dplyr::summarise(w=mean(value)),
        indep.weight_,
        reshape2::melt(pi.weight_) %>% dplyr::group_by(model) %>% dplyr::summarise(w=mean(value))
      ), id.vars = "model"
    ) %>% dplyr::mutate(type=c("perf","indep","PI")[L1]) %>% dplyr::select(-L1, -variable)
  
  if(arrange.by.pi){
    melted.total.weight$model <- 
      factor(melted.total.weight$model,
             levels = pi.weight$model[order((melted.total.weight %>% filter(type=="PI"))$value, decreasing = T)])
  }
  
  pi.weight.plot <- 
    ggplot2::ggplot(melted.total.weight %>% filter(type=="PI"), ggplot2::aes(x=model, y=value)) +
    ggplot2::geom_line(ggplot2::aes(group=1)) + 
    ggplot2::geom_point(ggplot2::aes(color="PI-weight"), size=7, shape=16) + 
    ggplot2::geom_point(data=melted.total.weight %>% filter(type=="perf"), ggplot2::aes(color="performance"), size=7, shape=8) +
    ggplot2::geom_point(data=melted.total.weight %>% filter(type=="indep"), ggplot2::aes(color="independence"), size=7, shape=17) +
    ggplot2::geom_hline(ggplot2::aes(yintercept=1/length(perf.weight_$model)), color="gray") + 
    ggplot2::scale_color_brewer(palette = "Dark2") +
    ggplot2::theme_classic() + 
    ggplot2::theme(plot.title = ggplot2::element_text(size=22),
                   axis.text.x = ggplot2::element_text(size=22),
                   axis.text.y = ggplot2::element_text(size=22),
                   axis.title.x = ggplot2::element_text(size=22),
                   axis.title.y = ggplot2::element_text(size=22),
                   axis.ticks.length= grid::unit(.25, "cm"),
                   strip.text.x = ggplot2::element_text(size=25),
                   strip.text.y = ggplot2::element_text(size=25),
                   legend.title = ggplot2::element_blank(),
                   legend.text = ggplot2::element_text(size=18),
                   legend.position = c(0.85,0.9),
                   legend.key = ggplot2::element_rect(fill = "white", colour = "black"))+
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(shape = c(17, 8, 16))))+
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust=0.1, hjust=0.99)) + 
    ggplot2::ylab("PI-weight") + 
    ggplot2::xlab('')
  
  pi.weight.plot
}


#' sigma.s 마다 entropy 계산한 것을 시각화
#'
#' @param res.sigma.s.obj 
#'
#' @return
#' @export
#'
#' @examples
plot.sigma.s.res <- function(res.sigma.s.obj){
  
  # usethis::use_package("ggplot2")
  
  tmp.sigma.s <- res.sigma.s.obj$tmp.sigma.s
  optim.sigma.s <- res.sigma.s.obj$sigma.S
  sigma.s <- res.sigma.s.obj$sigma.S$sigma.s
  
  sigma.s.plot <- 
    ggplot2::ggplot(tmp.sigma.s, ggplot2::aes(x=sigma.s, y=value)) +
    ggplot2::geom_line() +
    ggplot2::geom_segment(ggplot2::aes(x=optim.sigma.s$sigma.s, xend=optim.sigma.s$sigma.s,
                                       y=optim.sigma.s$value-0.001, yend=optim.sigma.s$value), color="gray", size=0.3)+
    ggplot2::xlab(expression(sigma[S])) + 
    ggplot2::ylab("Entropy") +
    ggplot2::scale_x_continuous(breaks = c(0.1, optim.sigma.s$sigma.s, 1.0)) +
    ggplot2::scale_y_continuous(expand = c(0,0)) +
    ## Theme
    ggplot2::theme_classic() + 
    ggplot2::theme(plot.title = ggplot2::element_text(size=22),
                   axis.text.x = ggplot2::element_text(size=22),
                   axis.text.y = ggplot2::element_text(size=22),
                   axis.title.x = ggplot2::element_text(size=22),
                   axis.title.y = ggplot2::element_text(size=22),
                   axis.ticks.length=grid::unit(.25, "cm"),
                   strip.text.x = ggplot2::element_text(size=25),
                   strip.text.y = ggplot2::element_text(size=25),
                   legend.title = ggplot2::element_blank(),
                   legend.text = ggplot2::element_text(size=18),
                   legend.position = c(0.85,0.9),
                   legend.key = ggplot2::element_rect(fill = "white", colour = "black"))
  sigma.s.plot
}


#' Title
#'
#' @param ref.rtlv 
#' @param his.rtlv.ensembled 
#' @param sce.rtlv.ensembled 
#' @param variable.name 
#' @param rt.name 
#' @param ylim 
#' @param label 
#' @param label.y.coord 
#' @param label.x.coord 
#' @param axis.label 
#' @param y.label 
#' @param calibrate.text.size 
#'
#' @return
#' @export
#'
#' @examples
rtlv.boxplot <- function(ref.rtlv, his.rtlv.ensembled, sce.rtlv.ensembled, variable.name, rt.name, 
                         ylim=c(50, 400),
                         label=c("SSP2-4.5","SSP3-7.0","SSP5-8.5"),
                         label.y.coord = 400,
                         label.x.coord = c(4,7,10),
                         axis.x.label=c("OBS","HIS(NBC)","P1","P2","P3","P1","P2","P3","P1","P2","P3"),
                         y.label = "Return Level",
                         calibrate.text.size=0){

  ref.melted <- reshape2::melt(ref.rtlv) %>% dplyr::mutate(period="p0", sce="ref0")
  his.melted <- reshape2::melt(his.rtlv.ensembled) %>% dplyr::mutate(period="p0", sce="ref1")
  sce.melted <- reshape2::melt(sce.rtlv.ensembled)
  
  colnames(ref.melted) <- c("loc","rt","value","var","period","sce")
  colnames(his.melted) <- c("loc","rt","value","var","period","sce")
  colnames(sce.melted) <- c("loc","rt","value","var","period","sce")
  
  all.melted <- rbind(ref.melted, his.melted, sce.melted) %>% mutate(plot.order=paste0(toupper(sce),toupper(period)))
  
  lb <- data.frame(x=label.x.coord, 
                   y=c(rep(label.y.coord,length(label))),
                   label=label,
                   variable = c(rep(rt.name,length(label))))
  
  rtlv.plot <-
    ggplot2::ggplot(all.melted %>%
                      dplyr::filter(rt==rt.name, var==variable.name), 
                    ggplot2::aes(x=plot.order, y=value, fill=sce)) +
    ggplot2::geom_hline(yintercept = c(100,200,300,400,500), color="lightgray")+
    ggplot2::geom_boxplot() +
    ggplot2::scale_fill_brewer(palette="Set2") +
    ggplot2::scale_x_discrete(labels=axis.x.label)+
    ggplot2::geom_vline(xintercept = c(3.5,6.5,9.5)-1, linetype="dashed", color="gray")+
    ggplot2::geom_label(data = lb, ggplot2::aes(x=x, y=y, label=label), inherit.aes = F, size=7)+
    ggplot2::theme_classic() + 
    ggplot2::theme(plot.title       = ggplot2::element_text(size=22+calibrate.text.size),
                   axis.text.x      = ggplot2::element_text(size=22+calibrate.text.size),
                   axis.text.y      = ggplot2::element_text(size=22+calibrate.text.size),
                   axis.title.x     = ggplot2::element_text(size=22+calibrate.text.size),
                   axis.title.y     = ggplot2::element_text(size=22+calibrate.text.size),
                   axis.ticks.length= ggplot2::unit(.25, "cm"),
                   strip.text.x     = ggplot2::element_text(size=25+calibrate.text.size),
                   strip.text.y     = ggplot2::element_text(size=25+calibrate.text.size),
                   legend.title     = ggplot2::element_blank(),
                   legend.text      = ggplot2::element_text(size=18+calibrate.text.size),
                   legend.position  = c(0.85,0.9),
                   legend.key        = ggplot2::element_rect(fill = "white", colour = "black")) +
    ggplot2::theme(legend.position = "none", axis.text.x = ggplot2::element_text(angle=0)) +
    ggplot2::ylim(ylim[1],ylim[2])+
    ggplot2::xlab("") +
    ggplot2::ylab(y.label)
  rtlv.plot
}


#' Title
#'
#' @param var.rtlv.obj 
#' @param variable 
#' @param rt 
#' @param variance 
#' @param calibrate.text.size 
#'
#' @return
#' @export
#'
#' @examples
var.scatter.plot <- function(var.rtlv.obj, var.name, rtlv, variance, calibrate.text.size=0){
  require(dplyr)
  
  if( nrow( var.rtlv.obj %>% dplyr::filter(variable==var.name) ) == 0) 
    stop("'var.name' is something wrong, please check it.")
  if( nrow( var.rtlv.obj %>% dplyr::filter(rt==rtlv) ) == 0) 
    stop("'rtlv' is something wrong, please check it.")
  
  base.plot <- 
    if(variance=="vbm"){
      ggplot2::ggplot(var.rtlv.obj %>% dplyr::filter(variable==var.name, rt==rtlv), 
                      ggplot2::aes(x=Vbm.eq.sum, y=Vbm.sum))+
        ggplot2::xlab(expression(Var[paste(BM,',',Eq)])) + 
        ggplot2::ylab(expression(Var[paste(BM,',PI')]))
    }else if(variance=="vwm"){
      ggplot2::ggplot(var.rtlv.obj %>% dplyr::filter(variable==var.name, rt==rtlv), 
                      ggplot2::aes(x=Vwm.eq.sum, y=Vwm.sum))+
        ggplot2::xlab(expression(Var[paste(WM,',',Eq)])) + 
        ggplot2::ylab(expression(Var[paste(WM,',PI')]))
    }else if(variance=="total"){
      ggplot2::ggplot(var.rtlv.obj %>% dplyr::filter(variable==var.name, rt==rtlv), 
                      ggplot2::aes(x=total.var.eq, y=total.var))+
        ggplot2::xlab(expression(Var[Eq])) + 
        ggplot2::ylab(expression(Var['PI']))
    }else{
      stop("'variance' must be the one of ('vbm','vwm','total')")
    }
  extra.plot <- 
    if(variance=="vbm"){
      ggplot2::theme(strip.text.x = ggplot2::element_blank())
    }else if(variance=="vwm"){
      ggplot2::theme(strip.text.x = ggplot2::element_blank())
    }else if(variance=="total"){
      ggplot2::theme(strip.text.x = ggplot2::element_text(size=15+calibrate.text.size))
    }
  
  base.plot+
    ggplot2::geom_point(ggplot2::aes(shape=period, color=period), size=1.3, alpha=0.7) +
    ggplot2::geom_abline(slope = 1, intercept = 0) +
    ggplot2::facet_wrap(sce~., scales = 'free')+
    ggplot2::theme(plot.title       = ggplot2::element_text(size=22+calibrate.text.size),
                   axis.text.x      = ggplot2::element_text(size=22+calibrate.text.size),
                   axis.text.y      = ggplot2::element_text(size=22+calibrate.text.size),
                   axis.title.x     = ggplot2::element_text(size=22+calibrate.text.size),
                   axis.title.y     = ggplot2::element_text(size=22+calibrate.text.size),
                   axis.ticks.length= ggplot2::unit(.25, "cm"),
                   strip.text.y     = ggplot2::element_text(size=25+calibrate.text.size),
                   legend.title     = ggplot2::element_blank(),
                   legend.text      = ggplot2::element_text(size=18+calibrate.text.size),
                   legend.position  = c(0.85,0.9),
                   legend.key       = ggplot2::element_rect(fill = "white", colour = "black"))+
    ggplot2::theme_bw()+
    ggplot2::theme(legend.position = 'right', 
                   legend.key   = ggplot2::element_rect(colour = "black")) +
    extra.plot
}


#' Title
#'
#' @param ensemble.sce.obj 
#' @param var.name 
#' @param rtlv 
#' @param facet.label 
#' @param axis.x.names 
#'
#' @return
#' @export
#'
#' @examples
rtpd.boxplot <- function(ensemble.sce.obj, var.name, rtlv, facet.label, axis.x.names,
                         calibrate.text.size=0){
  require(dplyr)
  
  rtpd.obj <- ensemble.sce.obj$return.period
  rtpd.gather <- reshape2::melt(rtpd.obj) %>% mutate(L4=paste0(L1,L2))
  
  rtpd.pi.plot <- 
    ggplot2::ggplot(rtpd.gather %>% dplyr::filter(Var2 %in% rtlv, L3==var.name),
                    ggplot2::aes(x=L4, y=value, fill=L2)) +
    ggplot2::geom_boxplot()+ 
    ggplot2::facet_wrap(~Var2, scales = "fix",
                        labeller = ggplot2::labeller(Var2=facet.label))+
    ggplot2::geom_vline(xintercept = c(4.5,7.5,10.5)-1, linetype="dashed", color="gray")+
    ggplot2::geom_point(color="blue", alpha=0.1)
  
  sce.names <- rtpd.gather %>% dplyr::filter(Var2 %in% rtlv) %>% dplyr::pull(L1) %>% unique
  
  for(i in 1:length(sce.names)){
    rtpd.pi.plot <- 
      rtpd.pi.plot +   
      ggplot2::geom_line(data=rtpd.gather %>% dplyr::filter(Var2 %in% rtlv, L1==sce.names[i]), 
                         ggplot2::aes(group=Var1), color="blue", alpha=0.1)
  }
  
  rtpd.pi.plot <- 
    rtpd.pi.plot +
    ggplot2::scale_fill_brewer(palette="Set3") +
    ggplot2::theme_classic() + 
    ggplot2::theme(plot.title       = ggplot2::element_text(size=22+calibrate.text.size),
                   axis.text.x      = ggplot2::element_text(size=22+calibrate.text.size),
                   axis.text.y      = ggplot2::element_text(size=22+calibrate.text.size),
                   axis.title.x     = ggplot2::element_text(size=22+calibrate.text.size),
                   axis.title.y     = ggplot2::element_text(size=22+calibrate.text.size),
                   axis.ticks.length= ggplot2::unit(.25, "cm"),
                   strip.text.x     = ggplot2::element_text(size=25+calibrate.text.size),
                   strip.text.y     = ggplot2::element_text(size=25+calibrate.text.size),
                   legend.title     = ggplot2::element_blank(),
                   legend.text      = ggplot2::element_text(size=18+calibrate.text.size),
                   legend.position = c(0.85,0.9),
                   legend.key = ggplot2::element_rect(fill = "white", colour = "black")) +
    ggplot2::theme_bw() + 
    ggplot2::ylab("Return Period") + 
    ggplot2::xlab("") +
    ggplot2::theme(legend.position = "none")+
    ggplot2::scale_x_discrete(labels=axis.x.names)
  rtpd.pi.plot
}
