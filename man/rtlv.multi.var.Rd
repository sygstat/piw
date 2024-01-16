\name{rtlv.multi.var}
\alias{rtlv.multi.var}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Return Level 계산하기
}
\description{
\code{readata.fromDirPath}, \code{readata.fromFilePath} 함수로 읽은 데이터에 대한 Return Level 을 계산하기 위한 함수이다. 
}
\usage{
rtlv.multi.var(data.obj, pfns, qfns, rt.lv = c(2, 5, 10, 20, 30, 50, 100))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data.obj}{
\code{readata.fromDirPath}, \code{readata.fromFilePath} 함수로 읽은 데이터 객체(list).
}
  \item{pfns}{
Parameter estimation for specific distributions by the method of L-moments. lmom 패키지에서 pelxxx 꼴로 구현된 함수를 list 로 묶어서 입력해야 한다. (반드시 data.obj가 가지고 있는 변수의 개수와 동일한 개수를 설정해야 한다.) 사용가능한 목록은 ?lmom::pelgev 를 콘솔에 실행해서 확인할 수 있다.
}
  \item{qfns}{
Quantile function for specific distributions by the method of L-moments. lmom 패키지에서 quaxxx 꼴로 구현된 함수를 list 로 묶어서 입력해야 한다. (반드시 data.obj가 가지고 있는 변수의 개수와 동일한 개수를 설정해야 한다.) 
}
  \item{rt.lv}{
계산할 return level. 예) c(20, 50) 을 인수로 주면 20년, 50년 Return Level이 계산됨.
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
계산된 Return Level이 변수 및 Return Level 개수에 따라 list 로 묶인 객체 반환된다. 
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
  \code{\link{lmom}}
  \code{\link{readata.fromFilePath}}
  \code{\link{readata.fromDirPath}}
}
\examples{
if(FALSE){ ## NOT RUN
  obs.file.path =  "/blabla/OBS.rds"
  obs = readata.fromFilePath(obs.file.path, 1973:2010, variable_names = c('AMP1','AMP5'), rfa9 = T)
  
  # Parameter estimation Functions by L-moments
  pfns = list(pelgev, pelnor)
  
  # Quantile Functions
  qfns = list(quagev, quanor)
  
  obs.rtlv = rtlv.multi.var(obs, pfns, qfns)
}
}