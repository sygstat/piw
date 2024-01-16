\name{calc.sigma.d}
\alias{calc.sigma.d}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
최적의 sigma.D 값 계산.
}
\description{
PI-weights 에서 P-weights 계산을 위한 적절한 \eqn{\sigma_D} 를 선택해야 하는데, 최적 \eqn{\sigma_D} 선택을 위한 함수이다. chi-square test 를 통하여 최적의 \eqn{\sigma_D}를 선택하며, 선택되는 \eqn{\sigma_D}는 주어진 \eqn{\sigma_D} 후보들 중에서 P-weights를 통계적으로 등가중치가 되지 않도록 하는 가장 큰 값이다. 
}
\usage{
calc.sigma.d(ref, model, simul.N = 10000, random.seed = 123456, significance.level = 0.05, candi.sigma.d = seq(1, 0.01, -0.01))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ref}{
model 자료의 reference 자료. \code{\link{readata.fromDirPath}} 또는 \code{\link{readata.fromFilePath}}를 통해 생성된 객체(\code{\link{tilde.normalization}} 도 가능). 예) 관측 데이터 
}
  \item{model}{
model 자료. \code{\link{readata.fromDirPath}} 또는 \code{\link{readata.fromFilePath}}를 통해 생성된 객체(\code{\link{tilde.normalization}} 도 가능).
}
  \item{simul.N}{
simulation의 횟수.
}
  \item{random.seed}{
Random seed.
}
  \item{significance.level}{
유의수준.
}
  \item{candi.sigma.d}{
\eqn{\sigma_D} 후보들. Vector.
}
}
\details{
To select an appropriate value of \eqn{\sigma_D} for the P-weights, we attempted to use the entropy criterian again, but we were not fortunate enough to obtain the optimal result, as with \eqn{\sigma_S}. Thus, a technique based on the \eqn{p}-value of the chi-square statistic is considered in this study.\cr
We denote \eqn{P_i(\sigma_D)} as a normalized P-weight for model \eqn{i} and for the given \eqn{\sigma_D}, which is defined as follows:\cr
\deqn{P_i(\sigma_D)=\frac{exp(-D_{i}/\sigma_D)}{\sum_{l=1}^{M} exp(-D_{i}/\sigma_D)}}
where \eqn{D_i} is the performance measure of the \eqn{i}-th model. \cr
For testing the hypothesis frame, the null hypothesis and the alternative hypothesis for \eqn{i=1, ..., M} are as follows:\cr
\deqn{H0 : all weights are equal}
\deqn{H1 : some weights are not equal}
For the given \eqn{P_i}, the chi-square statistic used to test the above hypothesis is as follows:
\deqn{\chi_{0}^{2}(\sigma_D)=\sum_{i=1}^{M}\frac{(1/m - P_i(\sigma_D))^2}{1/M}}.
Because we do not want to accept equal weights, \eqn{\sigma_D} should be selected to reject the null hypothesis. That is, the \eqn{p}-value obtained from the chi-square statistic should be less than a preassigned value \eqn{\alpha} (significance level), e.g., \eqn{\alpha=0.05}. In addition, because we also do not want aggressive weights, a \eqn{\sigma_D} can be selected as the maximum value of \eqn{\sigma_D} in which we still reject \eqn{H0} with \eqn{\alpha} level. That is, our select is 
\deqn{\sigma_D^*=max[\sigma_D:p-value(\sigma_D)<\alpha]},
where \eqn{p}-value \eqn{(\sigma_D)=Pr[\chi^2>\chi_0^2(\sigma_D)|H0]}.
}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{simul.res }{sigma.D 값에 따른 simulation 결과. 후보 sigma.D의 개수만큼 리스트 원소가 생성됨. list.}
  \item{var.agg }{각 variable 별, 최적 sigma.D 결과. data.frame.}
  \item{sigma.D }{최종 선택된 sigma.D 의 값. var.agg 의 sigma.D 값을 평균한 값임. numeric.}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
Atmosphere 2021, 12, 97. https://doi.org/10.3390/atmos12010097
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
\code{\link{readata.fromDirPath}}
\code{\link{readata.fromFilePath}}
\code{\link{tilde.normalization}}
}
\examples{
  if(FALSE){ ## NOT RUN
    obs.file.path = "/blabla/OBS.rds"
    obs = readata.fromFilePath(obs.file.path, 1973:2010, variable_names = c('AMP1','AMP5'), rfa9 = T)
    histDir = "...blabla/hisDir/"
    his = readata.fromDirPath(histDir   , 1973:2010, variable_names = c('AMP1','AMP5'), rfa9 = T)
    
    library(lmom)
    
    # Parameter estimation Functions by L-moments
    pfns <- list(pelgev, pelgev, pelnor, pelgev, pelgev)
    
    # Quantile Functions
    qfns <- list(quagev, quagev, quanor, quagev, quagev)
    
    obs.rtlv <- rtlv.multi.var(obs, pfns, qfns)
    his.rtlv <- rtlv.multi.var(his, pfns, qfns)
    his.rtlv.bc <- rtlv.multi.var(his, pfns, qfns)
    
    obs.tilde <- tilde.normalization(obs, his, target.rtlv = obs.rtlv)
    his.tilde <- tilde.normalization(obs, his, target.rtlv = his.rtlv)
  
    res.sigma.d <- calc.sigma.d(ref = obs.tilde, 
                                model = his.tilde, 
                                simul.N = 10000, 
                                random.seed = 197079,
                                significance.level = 0.05,
                                candi.sigma.d = seq(1.0, 0.01, -0.01))
  }
}