\name{tilde.normalization}
\alias{tilde.normalization}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
데이터 정규화
}
\description{
Return Level 값을 이용해서 ref 데이터와 model 데이터의 정규화를 수행한다. 이 정규화는 PI-weight 를 계산할때 Sigma.D 값에 따른 Performance 가중값을 안정시킴.
}
\usage{
tilde.normalization(ref, model, target.rtlv)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ref}{
%%     ~~Describe \code{ref} here~~
}
  \item{model}{
%%     ~~Describe \code{model} here~~
}
  \item{target.rtlv}{
%%     ~~Describe \code{target.rtlv} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (ref, model, target.rtlv) 
{
    if (!all(names(ref) == names(model))) 
        stop("The variables order of \"ref\" and \"model\" should be the same.")
    cat("Calc \"Normalization Factors\"\n")
    norm.factors <- pbapply::pbmapply(function(ref_, model_) {
        total.df <- rbind(ref_, do.call(rbind, model_))
        normalize.factor <- data.frame(med = apply(total.df, 
            2, median, na.rm = T), max = apply(total.df, 2, max, 
            na.rm = T), min = apply(total.df, 2, min, na.rm = T))
        normalize.factor
    }, ref, model, SIMPLIFY = FALSE)
    cat("Normalize Target\n")
    if (!is.null(dim(target.rtlv[[1]]))) {
        mapply(function(target.rtlv_, norm.factors_) {
            tmp <- pbapply::pblapply(X = 1:nrow(target.rtlv_), 
                FUN = function(i.loc) {
                  rt.T <- target.rtlv_[i.loc, ]
                  med_ <- norm.factors_$med[i.loc]
                  min_ <- norm.factors_$min[i.loc]
                  max_ <- norm.factors_$max[i.loc]
                  norm.value <- rt.T
                  norm.value[rt.T >= med_] <- (norm.value[rt.T >= 
                    med_] - med_)/(max_ - med_)
                  norm.value[rt.T < med_] <- (norm.value[rt.T < 
                    med_] - med_)/(med_ - min_)
                  norm.value
                })
            tmp_ <- do.call(rbind, tmp)
            rownames(tmp_) <- paste0("G", 1:nrow(tmp_))
            tmp_
        }, target.rtlv, norm.factors, SIMPLIFY = FALSE)
    }
    else {
        mapply(function(multi.target.rtlv_, norm.factors_) {
            pbapply::pblapply(multi.target.rtlv_, function(target.rtlv_) {
                tmp <- lapply(X = 1:nrow(target.rtlv_), FUN = function(i.loc) {
                  rt.T <- target.rtlv_[i.loc, ]
                  med_ <- norm.factors_$med[i.loc]
                  min_ <- norm.factors_$min[i.loc]
                  max_ <- norm.factors_$max[i.loc]
                  norm.value <- rt.T
                  norm.value[rt.T >= med_] <- (norm.value[rt.T >= 
                    med_] - med_)/(max_ - med_)
                  norm.value[rt.T < med_] <- (norm.value[rt.T < 
                    med_] - med_)/(med_ - min_)
                  norm.value
                })
                tmp_ <- do.call(rbind, tmp)
                rownames(tmp_) <- paste0("G", 1:nrow(tmp_))
                tmp_
            })
        }, target.rtlv, norm.factors, SIMPLIFY = FALSE)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
