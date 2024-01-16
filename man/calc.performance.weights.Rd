\name{calc.performance.weights}
\alias{calc.performance.weights}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
P-weights 계산.
}
\description{
주어진 자료와 sigma.D 값으로부터 P-weights 를 계산한다.
}
\usage{
calc.performance.weights(ref, model, sigma.D)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ref}{
model 자료의 reference 자료. \code{\link{readata.fromDirPath}} 또는 \code{\link{readata.fromFilePath}}를 통해 생성된 객체(\code{\link{tilde.normalization}} 도 가능). 예) 관측 데이터 
}
  \item{model}{
model 자료. \code{\link{readata.fromDirPath}} 또는 \code{\link{readata.fromFilePath}}를 통해 생성된 객체(\code{\link{tilde.normalization}} 도 가능).
}
  \item{sigma.D}{
실수값.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{weight.df }{모델별 가중값 테이블. columns은 격자점이 된다. 단, 첫 번째 column은 model 이름.data.frame.}
  \item{weight.list }{모델별 가중값 리스트. list 원소가 격자점이 된다. list.}
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
}
\examples{
if(FALSE){##NOT RUN
  obs.file.path = "/blabla/OBS.rds"
  obs = readata.fromFilePath(obs.file.path, 1973:2010, variable_names = c('AMP1','AMP5'), rfa9 = T)
  
  histDir = "...blabla/hisDir/"
  his = readata.fromDirPath(histDir   , 1973:2010, variable_names = c('AMP1','AMP5'), rfa9 = T)
  
  calc.performance.weights(ref=obs, model = his, sigma.D=0.5)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
