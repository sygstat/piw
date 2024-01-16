\name{readata.fromDirPath}
\alias{readata.fromDirPath}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
파일이 있는 폴더경로로부터 데이터를 불러오기.
}
\description{
\code{dir.path}으로 넘어온 경로의 파일을 불러온다. 반드시 폴더경로이어야 한다. 폴더내의 파일은 확장자는 모두 .rds 이어야 하고, JY-type의 데이터이어야 한다. 
JY-type 데이터는 ['lon', 'lat', 'time', 'info','변수명']을 갖는 list 형식의 객체이다. 
'lon'  : vector 형태의 경도(longitude) 정보,
'lat'  : vector 형태의 위도(latitude) 정보,
'time' : vector 형태의 time 정보, 
'info' : 데이터의 summary 정보(이 함수에서는 이 데이터를 사용하지 않으므로, 어떤 형태의 값으로 구성되든 상관 없다), 
'변수명' : 변수 데이터의 값, 예) AMP1, AMP5 등.(변수는 여러개가 있을수 있음)
}
\usage{
readata.fromDirPath(dir.path, filter.period, variable_names, filter.grid = NULL, rfa9 = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dir.path}{
JY-type 데이터로 구성된 rds 파일이 있는 "폴더"의 경로. 하나의 폴더 경로만 다룰수 있다.
}
  \item{filter.period}{
필터링 할 시간. vector 타입이어야 하며, 파일 내의 'time' vector에 존재하는 값만 필터링함.
}
  \item{variable_names}{
읽어올 변수명, 여러개를 설정할 수 있음. 예) c('AMP1', AMP5')
}
  \item{filter.grid}{
전체 격자점의 정보가 아닌 일부 격자점의 정보만 읽을 경우 사용. default 는 NULL 이고, 이는 전체 격자 정보를 읽는다. 격자점의 정보는 격자의 index 번호를 vector 타입으로 설정하면 된다.
}
  \item{rfa9}{
RFA 적용 여부, Default는 FALSE 이며, 이는 RFA를 적용하지 않는 것이다.
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
\code{variable_names}에 설정된 변수개수의 크기를 갖는 리스트 객체가 반환된다. 리스트의 각각의 원소는 폴더내의 파일들의 이름으로 구성된 리스트이며, 이것들의 원소는 data.frame 형태의 데이터이다. 
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
  \code{\link{readata.fromFilePath}}
}
\examples{
if(FALSE){ ## NOT RUN
  histDir = "...blabla/hisDir/"
  his = readata.fromDirPath(histDir   , 1973:2010, variable_names = c('AMP1','AMP5'), rfa9 = T)
}
}