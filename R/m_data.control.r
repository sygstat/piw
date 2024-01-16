##################################################
## DATA 읽어 오는 모듈들
## by Yonggwan, 2021.05.13
##################################################


#' 변수별 모델의 개수 및 순서가 동일한지 체크. 개수 및 순서가 동일하지 않으면 stop
#'
#' @param model.data 
#'
#' @return
#' @export
#'
#' @examples
check.valid <- function(model.data){
  cat("Validation Check...")
  if((length(model.data)>1) && ( is.null(dim(model.data[[1]])) )){
    number.of.models <- do.call(c, lapply(model.data, function(x) length(x)))
    size.check <- number.of.models == max(number.of.models)
    if(!all(size.check)){
      cat("# of models by variable\n")
      print(number.of.models)
      stop("The number of models by variable is different.")
    }   
    
    models.tmp <- data.frame(t(do.call(rbind, lapply(model.data, function(x) names(x)))))
    
    if(!all(apply(models.tmp, MARGIN = 2, function(x) x==models.tmp[[1]]))){
      print(models.tmp)
      stop("The order of models by variable is different.")
    }
  }
  cat("Ok\n")
}

#' Title
#'
#' @param gdata 
#' @param filter.period 
#' @param variable_name 
#' @param filter.gird 
#'
#' @return
#' @export
#'
#' @examples
read.basic <- function(gdata, filter.period, variable_name, filter.grid=NULL){
  
  data.check <- c("lon","lat","time") %in% names(gdata)
  if(!all(data.check)){
    stop(paste0("Check your file - Unvailable of attribute : ",c("lon","lat","time")[!data.check]))
  }
  var.check <- variable_name %in% names(gdata)
  if(!all(var.check)){
    stop(paste0("Check your file - Unvailable of variable : ",variable_name[!var.check]))
  }
  
  # filter.period 에 설정된 기간이 데이터에 존재하지 않을 경우
  if(sum(gdata$time %in% filter.period) != length(filter.period)){
    stop("The time period set in 'filter.period' does not exist in the data.")
  }
  
  target <- gdata[[variable_name]][,,gdata$time %in% filter.period]
  
  res <- matrix(target, nrow=dim(target)[3], ncol = prod(dim(target)[1:2]), byrow = T)
  
  res <- res[,complete.cases(t(res))] # remove NA
  colnames(res) <- paste0("G",1:ncol(res))
  
  # 특정 그리드만 뽑을 경우
  if(!is.null(filter.grid)){
    res <- res[,filter.grid]
  }
  res
}


#' Title
#'
#' @param gdata JY data
#' @param filter.period P
#' @param variable_name V
#'
#' @return
#' @export
#'
#' @examples
rfa9.basic <- function(gdata, filter.period, variable_name, filter.grid=NULL){
  
  data.check <- c("lon","lat","time") %in% names(gdata)
  if(!all(data.check)){
    stop(paste0("Check your file - Unvailable of attribute : ",c("lon","lat","time")[!data.check]))
  }
  var.check <- variable_name %in% names(gdata)
  if(!all(var.check)){
    stop(paste0("Check your file - Unvailable of variable : ",variable_name[!var.check]))
  }
  
  # filter.period 에 설정된 기간이 데이터에 존재하지 않을 경우
  if(sum(gdata$time %in% filter.period) != length(filter.period)){
    stop("The time period set in 'filter.period' does not exist in the data.")
  }
  
  target <- gdata[[variable_name]][,,gdata$time %in% filter.period]
  
  max.i <- nrow(target)
  max.j <- ncol(target)
  
  rfa.data <- matrix(nrow=dim(target)[3]*9)
  
  # i : matrix 에서 위에서 아래로 움직임
  # j : matrix 에서 왼쪽에서 오른쪽으로 움직임
  for(j in 1:max.j){
    for(i in 1:max.i){
      # 현재 셀(i,j, )의 값 중에 하나라도 NA 인 것은 제외
      if(all(!is.na(target[i,j, ]))){
        rfa.values <- 
          lapply(
            list(c(i,j),c(i+1, j),c(i-1, j),c(i , j+1),c(i , j-1),c(i+1, j+1),c(i+1, j-1),c(i-1, j+1),c(i-1, j-1)),
            function(cell){
              if(all(cell>0) & cell[1] <= max.i & cell[2] <= max.j){ # 범위 벗어나지 않는 셀만
                value <- target[cell[1], cell[2], ]
                value
              }else{ # 범위 벗어날경우 NA 반환
                rep(NA,dim(target)[3])
              }
            })
        rfa.data <- cbind(rfa.data, do.call(c,rfa.values))
      }
    }
  }
  rfa.data <- rfa.data[,-1]
  colnames(rfa.data) <- paste0("G",1:ncol(rfa.data))
  
  # 특정 그리드만 뽑을 경우
  if(!is.null(filter.grid)){
    rfa.data <- rfa.data[,filter.grid]
  }
  
  
  return(as.matrix(rfa.data))
}


#' Title
#'
#' @param file.path.list 
#' @param filter.period 
#' @param variable_name 
#'
#' @return
#' @export
#'
#' @examples
readata.fromFilePath <- function(file.path.list, filter.period, variable_names, filter.grid=NULL, rfa9=FALSE){
  
  # Model list
  model.list <- sort(gsub(x=file.path.list, pattern = ".*/|.rds", replacement = ""))
  
  model.data <-
    pbapply::pblapply(file.path.list, function(model.path){
      readr::read_rds(model.path)
    })
  names(model.data) <- model.list
  
  res <- 
    lapply(variable_names, function(variable_name){
      tmp.data <- 
        lapply(model.data, function(m.data){
          if(rfa9) rfa9.basic(m.data, filter.period, variable_name, filter.grid)
          else read.basic(m.data, filter.period, variable_name, filter.grid)
        })
      if(length(model.list)>1) names(tmp.data) <- model.list
      else tmp.data <- tmp.data[[1]]
      tmp.data
    })
  
  names(res) <- variable_names
  check.valid(res)
  res
}


#' Title
#'
#' @param dir.path 
#' @param filter.period 
#' @param variable_names 
#'
#' @return
#' @export
#'
#' @examples
readata.fromDirPath <- function(dir.path, filter.period, variable_names, filter.grid=NULL, rfa9=FALSE){
  
  model.file.list <- list.files(dir.path)
  
  model.list <- sort(gsub(x=model.file.list, pattern = ".rds", replacement = ""))
  
  model.data <-
    pbapply::pblapply(model.file.list, function(model.file){
      readr::read_rds(file.path(dir.path,model.file))
    })
  names(model.data) <- model.list
  
  res <- 
    pbapply::pblapply(variable_names, function(variable_name){
      tmp.data <- 
        lapply(model.data, function(m.data){
          if(rfa9) rfa9.basic(m.data, filter.period, variable_name, filter.grid)
          else read.basic(m.data, filter.period, variable_name, filter.grid)
        })
      if(length(model.list)>1) names(tmp.data) <- model.list
      else tmp.data <- tmp.data[[1]]
      tmp.data
    })
  # res <- 
  #   lapply(variable_names, function(variable_name){
  #     # Model list
  #     
  #     model.data <-
  #       pbapply::pblapply(model.file.list, function(model.file){
  #         tmp <- readr::read_rds(file.path(dir.path,model.file))
  #         if(rfa9) rfa9.basic(tmp, filter.period, variable_name, filter.grid)
  #         else read.basic(tmp, filter.period, variable_name, filter.grid)
  #       })
  #     if(length(model.list)>1) names(model.data) <- model.list
  #     else model.data <- model.data[[1]]
  #     model.data
  #   })
  
  names(res) <- variable_names
  check.valid(res)
  res
}

#' Title
#'
#' @param data.model 
#'
#' @return
#' @export
#'
#' @examples
get.model.length <- function(data.model){
  if(is.null(dim(data.model[[1]]))){
    length(data.model[[1]])
  }else{
    stop("No models")
  }
}

#' Title
#'
#' @param data.model 
#'
#' @return
#'
#' @examples
get.model.names <- function(data.model){
  if(is.null(dim(data.model[[1]]))){
    names(data.model[[1]])
  }else{
    stop("No models")
  }
}
