library(spatstat)
library(spatstat.linnet)
library(Rcpp)

## 도메인 트리에 대해 BFS를 수행한 후, 탐색 순서와 직전 객체, 간선의 탐색 방향을 기록한다.
treeparty.visit <- function(x, root = 1)
{
  if (!("lpp" %in% class(x)))
    stop("[x] should be an 'lpp' object.")
  res <- list()
  res$size.v <- x$domain$vertices$n
  res$size.e <- x$domain$lines$n
  if (res$size.v != res$size.e + 1)
    stop("Not a tree")
  adjlist <- vector("list", res$size.v) ## 각 정점에 인접한 간선 목록
  for (i in 1:res$size.v)
    adjlist[i] <- c(1)
  for (i in 1:res$size.e)
  {
    j <- x$domain$from[i]
    adjlist[[j]][1] <- adjlist[[j]][1] + 1
    adjlist[[j]][adjlist[[j]][1]] <- i
    j <- x$domain$to[i]
    adjlist[[j]][1] <- adjlist[[j]][1] + 1
    adjlist[[j]][adjlist[[j]][1]] <- i
  }
  res$queue <- rep(0L, res$size.v + res$size.e) ## 간선과 정점을 탐색한 순서. 양수일 경우 정점, 음수일 경우 간선.
  i <- 1 ## 탐색할 간선 또는 정점의 순서
  j <- 1 ## 지금까지 확인한 간선과 정점의 개수
  parent.v <- rep(0L, res$size.v)
  parent.e <- rep(0L, res$size.e)
  res$direction <- rep(FALSE, res$size.e) ## 간선의 탐색 방향. TRUE일 경우 from에서 to, FALSE일 경우 to에서 from
  res$queue[1] <- root ## <root>번 정점부터 시작
  parent.v[root] <- root
  while (i <= length(res$queue))
  {
    if (res$queue[i] > 0) ## 정점
    {
      if (adjlist[[res$queue[i]]][1] < 2) ## isolated vertex
        next
      for (k in 2:adjlist[[res$queue[i]]][1])
      {
        if (parent.e[adjlist[[res$queue[i]]][k]] != 0)
          next
        parent.e[adjlist[[res$queue[i]]][k]] <- res$queue[i]
        j <- j + 1
        res$queue[j] <- -adjlist[[res$queue[i]]][k]
        res$direction[-res$queue[j]] <- x$domain$from[-res$queue[j]] == res$queue[i]
      }
    }
    else ## 간선
    {
      if (parent.v[x$domain$from[-res$queue[i]]] == 0)
      {
        parent.v[x$domain$from[-res$queue[i]]] <- -res$queue[i]
        j <- j + 1
        res$queue[j] <- x$domain$from[-res$queue[i]]
      }
      if (parent.v[x$domain$to[-res$queue[i]]] == 0)
      {
        parent.v[x$domain$to[-res$queue[i]]] <- -res$queue[i]
        j <- j + 1
        res$queue[j] <- x$domain$to[-res$queue[i]]
      }
    }
    i <- i + 1
  }
  parent.v[root] <- 0L
  res$queue <- as.integer(res$queue)
  res$parent <- as.integer(parent.v[parent.e])
  class(res) <- "treeparty.visit"
  return(res)
}

## 관측치를 도메인 트리에 맞추어 트리의 형태로 배치한다.
## 데이터 트리에 대해 DFS를 수행하여 그 순서와 서브트리의 크기를 기록해둔다.
## 어떠한 두 관측치도 동일한 위치에 있지 않음을 전제로 한다.
## <x$domain>과 <visit>에서 사용한 <domain>은 동일함을 전제로 한다.
## 관측치를 변경할 경우, build부터 다시 수행해야 한다.
## <marks>는 <x$data$marks>를 벡터로 변환한 것이며,
## hyperframe의 개별 탐색 성능이 좋지 않은 문제를 해결했다.
## 
treeparty.build <- function(x, visit)
{
  ## 1 관측치를 트리 형태로 배치한다. 이때 관측치 트리에 junction이 존재하나 도메인 트리에서 간극이 확인될 경우, null vertex를 추가한다.
  ## 재귀함수를 사용하지 않는다.
  if (class(visit) != "treeparty.visit")
    stop("[visit] should be a 'treeparty.visit' object.")
  if (!("lpp" %in% class(x)))
    stop("[x] should be an 'lpp' object.")
  if (class(x$data$marks) != "factor")
    stop("[x$data$marks] should be a 'factor' object.")
  if (visit$size.v != x$domain$vertices$n || visit$size.e != x$domain$lines$n)
    stop("The size of [visit] and [x] are different.")
  #plot(x, legend = F, main = "")
  order.e <- 1:visit$size.e ## 큐에서 각 간선이 위치한 번호
  for (i in 1:(visit$size.v + visit$size.e))
    if (visit$queue[i] < 0)
      order.e[-visit$queue[i]] <- i
  ## > dendrite$data[1]
  ## Hyperframe:
  ##          x        y seg        tp marks
  ## 1 45.07634 220.7606 156 0.4087487  thin
  obs <- data.frame(index = 1L:nrow(x$data), seg = x$data$seg, tp = x$data$tp)
  res <- list()
  res$obs <- obs
  res$marks <- as.data.frame(x$data[,5])[[1]]
  ## BFS에서 방문한 간선의 역순으로 관측치 나열하기
  ## <root>에서 멀리 떨어진 관측치일수록 앞에 온다.
  for (i in 1:nrow(obs))
    if (visit$direction[obs$seg[i]])
      obs$tp[i] <- -obs$tp[i]
  obs <- obs[order(-order.e[obs$seg], obs$tp),]
  obs$tp <- abs(obs$tp)
  res$parent <- rep(0L, nrow(x$data))
  ## 각 간선에서 <root>에 가장 가까운 관측치를 임시로 저장하는 큐
  ## <qi>는 obs index, <qp>는 parent의 edge index를 저장하며,
  ## 구조체 큐를 구현하기 귀찮은 관계로 큐를 두 개 사용했다.
  ## <jf>과 <jb>는 큐의 front와 back의 위치를 의미한다.
  qi <- rep(0L, nrow(obs))
  qp <- rep(0L, nrow(obs))
  jf <- 1L
  jb <- 1L
  for (i in 1:(nrow(x$data) - 1))
  {
    idx <- obs$index[i]
    if (obs$seg[i] == obs$seg[i + 1])
    {
      res$parent[idx] <- obs$index[i + 1]
      next
    }
    qi[jb] <- idx
    qp[jb] <- visit$parent[obs$seg[i]]
    #plot(x$domain$lines[qp[jb]], add = TRUE, col = "pink")
    jb <- jb + 1L
    while (jf < jb)
    {
      if (order.e[qp[jf]] < order.e[obs$seg[i + 1]])
        break
      else if (qp[jf] == obs$seg[i + 1])
      {
        ## 큐에 관측치가 둘 이상이고,
        ## 두 관측치의 부모 간선이 동일하며,
        ## 해당 부모 간선의 끝점 중 <root>와 멀리 있는 끝점에 또 다른 관측치가 없어야 함
        ## 이 조건을 모두 만족할 경우,
        ## 가짜 관측치를 만들어 두 관측치의 부모로 둔다.
        if (jb == jf + 1)
          res$parent[qi[jf]] <- obs$index[i + 1]
        else if (qp[jf] != qp[jf + 1])
          res$parent[qi[jf]] <- obs$index[i + 1]
        else if (obs$seg[i + 1] != qp[jf] ||
                 obs$tp[i + 1] == visit$direction[obs$seg[i + 1]])
          res$parent[qi[jf]] <- obs$index[i + 1]
        else
        {
          res$parent[length(res$parent) + 1] <- qp[jf] # 임시 데이터
          res$obs[length(res$parent),] <- list(length(res$parent), obs$seg[i + 1], visit$direction[obs$seg[i + 1]] + 0)
          while (jf < jb)
            if (qp[jf] == res$parent[length(res$parent)])
            {
              res$parent[qi[jf]] <- length(res$parent)
              jf <- jf + 1L
            }
            else
              break
          jf <- jf - 1L
          res$parent[length(res$parent)] <- obs$index[i + 1]
        }
        jf <- jf + 1L
      }
      else
      {
        ## 큐의 앞에 있는 관측치가 아직 탐색하지 않은 관측치보다
        ## <root>에서 한참 멀리 있다.
        ## 큐의 앞에 있는 관측치를 맨 뒤로 보내거나,
        ## 여러 관측치의 공통 조상을 가짜 관측치로 생성한다.
        ## 생성한 가짜 관측치는 큐에 다시 넣는다.
        if (jb == jf + 1) ## 큐에 관측치가 하나뿐임
          qp[jf] <- visit$parent[qp[jf]]
        else if (qp[jf] != qp[jf + 1])
        {
          qi[jb] <- qi[jf]
          qp[jb] <- visit$parent[qp[jf]]
          ##plot(x$domain$lines[qp[jb]], add = TRUE, col = "pink")
          jf <- jf + 1L
          jb <- jb + 1L
        }
        else
        {
          res$parent[length(res$parent) + 1] <- qp[jf] # 임시 데이터
          res$obs[length(res$parent),] <- list(length(res$parent), qp[jf], visit$direction[qp[jf]] + 0)
          while (jf < jb)
            if (qp[jf] == res$parent[length(res$parent)])
            {
              res$parent[qi[jf]] <- length(res$parent)
              jf <- jf + 1L
            }
            else
              break
          qp[jb] <- visit$parent[res$parent[length(res$parent)]]
          qi[jb] <- length(res$parent)
          res$parent[length(res$parent)] <- 0L
          ##plot(x$domain$lines[qp[jb]], add = TRUE, col = "pink")
          jb <- jb + 1L
        }
      }
    }
  }
  ## 
  i <- nrow(x$data)
  idx <- obs$index[i]
  qi[jb] <- idx
  qp[jb] <- visit$parent[obs$seg[i]]
  ##plot(x$domain$lines[qp[jb]], add = TRUE, col = "pink")
  jb <- jb + 1L
  while (jf < jb + 1)
  {
    ##order.e[obs$seg[i + 1]] == -inf로 취급함
    if (jb == jf + 1) ## 큐에 관측치가 하나뿐이면 루트
    {
      res$root <- qi[jf]
      res$parent[res$root] <- 0L
      break
    }
    else if (qp[jf] != qp[jf + 1])
    {
      qi[jb] <- qi[jf]
      qp[jb] <- visit$parent[qp[jf]]
      ##plot(x$domain$lines[qp[jb]], add = TRUE, col = "pink")
      jf <- jf + 1L
      jb <- jb + 1L
    }
    else
    {
      res$parent[length(res$parent) + 1] <- qp[jf] # 임시 데이터
      res$obs[length(res$parent),] <- list(length(res$parent), qp[jf], visit$direction[qp[jf]] + 0)
      while (jf < jb)
        if (qp[jf] == res$parent[length(res$parent)])
        {
          res$parent[qi[jf]] <- length(res$parent)
          jf <- jf + 1
        }
        else
          break
      qi[jb] <- length(res$parent)
      qp[jb] <- visit$parent[res$parent[length(res$parent)]]
      res$parent[length(res$parent)] <- 0L
      ##plot(x$domain$lines[qp[jb]], add = TRUE, col = "pink")
      jb <- jb + 1
    }
  }
  res$order <- rep(0L, nrow(res$obs)) ## DFS를 통한 방문 순서
  res$invorder <- rep(0L, nrow(res$obs)) ## 의 역함수
  is.branch <- rep(FALSE, nrow(res$obs))
  adjlist <- vector("list", nrow(res$obs)) ## 각 정점에 인접한 간선 목록
  for (i in 1:nrow(res$obs))
    adjlist[i] <- c(1L)
  for (i in 1:nrow(res$obs))
  {
    j <- res$parent[i]
    if (j == 0)
      next
    adjlist[[j]][1] <- adjlist[[j]][1] + 1L
    adjlist[[j]][adjlist[[j]][1]] <- i
  }
  jt <- 1L
  jb <- nrow(res$obs)
  res$order[1] <- res$root
  
  
  while (jt > 0)
  {
    #cat('\n')
    if (adjlist[[res$order[jt]]][1] == 1)
    {
      while (jt > 0)
      {
        res$order[jb] <- res$order[jt]
        res$invorder[res$order[jb]] <- jb
        #cat(" >>", res$order[jb])
        jb <- jb - 1L
        jt <- jt - 1L
        if (is.branch[jt + 1])
          break
      }
    }
    else
    {
      jt.next <- jt
      for (j in 2:adjlist[[res$order[jt]]][1])
      {
        jt.next <- jt.next + 1L
        res$order[jt.next] <- adjlist[[res$order[jt]]][j]
        is.branch[jt.next] <- (j > 2)
        #cat(" <<", res$order[jt.next])
      }
      jt <- jt.next
      if (jt > jb)
        stop("!")
    }
  }
  res$subtree.size <- rep(1L, nrow(res$obs))
  for (i in res$order[nrow(res$obs):2])
    res$subtree.size[res$parent[i]] <- res$subtree.size[res$parent[i]] + res$subtree.size[i]
  class(res) <- "treeparty.build"
  return(res)
}

## 데이터 트리에서 각 관측치 정점을 루트로 하는 서브트리에 대해
## 종류별 관측치의 가중치 합을 구한다.
## 가중치를 변경할 경우, count부터 다시 수행해야 한다.
treeparty.count <- function(build, weight = rep(1, length(build$marks)))
{
  if (class(treeparty_count) != "function")
    stop("Source 'treeparty.cpp' is not compiled.")
  if (class(build) != "treeparty.build")
    stop("[build] should be a 'treeparty.build' object.")
  ## 데이터 트리에서 각 관측치 정점을 루트로 하는 서브트리에 대해
  ## 종류별 관측치의 가중치 합 카운트
  res <- list()
  res$labels <- factor(levels(build$marks), levels(build$marks))
  res$count <- data.frame(treeparty_count(as.integer(build$marks), length(res$labels), build$order, build$parent, weight)) # use Rcpp
  colnames(res$count) <- res$labels
  res$root <- build$root
  res$obs <- sum(res$count[res$root,])
  class(res) <- "treeparty.count"
  return(res)
}

## legacy version
#treeparty.count.R <- function(build, weight = rep(1, length(build$marks)))
#{
#  if (class(build) != "treeparty.build")
#    stop("[build] should be a 'treeparty.build' object.")
#  ## 데이터 트리에서 각 관측치 정점을 루트로 하는 서브트리에 대해
#  ## 종류별 관측치의 가중치 합 카운트
#  res <- list()
#  res$labels <- as.factor(levels(build$marks))
#  res$count <- data.frame(matrix(0, length(build$order), length(res$labels)))
#  colnames(res$count) <- res$labels
#  ## 현재 가장 느린 부분으로, R에서 for문 실행 속도가 느린 것으로 판단된다.
#  ## 만약 C++로 짤 수 있다면?
#  marks <- as.integer(build$marks)
#  for (i in 1:length(build$marks))
#  {
#    label <- marks[i]
#    res$count[i, label] <- weight[i]
#  }
#  for (i in length(build$order):2)
#    res$count[build$parent[build$order[i]],] <- res$count[build$parent[build$order[i]],] + res$count[build$order[i],]
#  res$root <- build$root
#  res$obs <- sum(res$count[res$root,])
#  class(res) <- "treeparty.count"
#  return(res)
#}

## 관측치 트리에 대해 결정 트리를 생성한다.
## 데이터셋을 분할하는 기준을 변경할 경우, split부터 다시 수행해야 한다.
## 개선 사항 : accuracy와 gini에 한해, apply 함수를 사용하여 계산 속도를 높였다.
treeparty.split <- function(count, minbucket = 30, eval = "accuracy", significance = 0.01)
{
  if (class(count) != "treeparty.count")
    stop("[count] should be a 'treeparty.count' object.")
  if (!(eval %in% c("accuracy", "chisq", "gini", "entropy")))
    stop("Invalid evaluation. [eval] should be 'accuracy', 'chisq', 'gini', or 'entropy.'")
  res <- list()
  res$root <- count$root
  res$eval <- eval
  res$subtree <- NA
  res$supertree <- NA
  res$pred <- NA
  res$pred.sub <- NA
  res$pred.sup <- NA
  class(res) <- "treeparty.split"
  if (count$obs <= minbucket)
  {
    res$pred <- count$labels[which.max(count$count[res$root,])] ## 가장 많은 관측치 선택, 동점 시 먼저 오는 거
    return(res)
  }
  else if (sum(count$count[res$root,] > 0) == 1)
  {
    if (eval == "accuracy")
      res$accuracy <- 1.0
    else if (eval == "chisq")
      res$p.value <- 1.0
    else
      res$impurity <- 0.0
    res$pred <- count$labels[count$count[res$root,] > 0]
    return(res)
  }
  if (eval == "chisq")
  {
    checked <- 0
    max.stat <- 0
    for (i in 1:nrow(count$count))
    {
      test.stat <- 0
      subtree <- sum(count$count[i,])
      if (subtree == 0 || subtree == count$obs)
        next
      for (j in 1:length(count$count))
      {
        expected <- subtree * count$count[count$root, j] / count$obs
        test.stat <- test.stat + (count$count[i, j] - expected) ^ 2 / expected
        expected <- count$count[count$root, j] - expected
        test.stat <- test.stat + (count$count[count$root, j] - count$count[i, j] - expected) ^ 2 / expected
      }
      if (test.stat > max.stat)
      {
        max.stat <- test.stat
        res$root <- i
      }
      checked <- checked + 1
    }
    res$p.value <- pchisq(max.stat, length(count$count) - 1, lower.tail = FALSE)
    if (res$p.value >= significance)
    {
      res$root <- count$root
      res$pred <- count$labels[which.max(count$count[res$root,])] ## 가장 많은 관측치 선택, 동점 시 먼저 오는 거
    }
    else
    {
      res$pred.sub <- count$labels[which.max(count$count[res$root,])]
      suptree <- count$count[count$root,] - count$count[res$root,]
      res$pred.sup <- count$labels[which.max(suptree)]
    }
  }
  else if (eval == "accuracy")
  {
    #remaining <- count$count
    #for (j in 1:length(count$count))
    #  remaining[,j] <- count$count[count$root, j] - remaining[,j]
    remaining <- rep(as.vector(unlist(count$count[count$root,])), each = nrow(count$count)) - count$count
    
    test.stats <- apply(count$count, 1, max) + apply(remaining, 1, max)
    max.stat <- max(test.stats)
    res$root <- which.max(test.stats)
    res$accuracy <- max.stat / count$obs
    res$pred.sub <- count$labels[which.max(count$count[res$root,])]
    suptree <- count$count[count$root,] - count$count[res$root,]
    res$pred.sup <- count$labels[which.max(suptree)]
  }
  else if (eval == "gini")
  {
    test1 <- apply(count$count ^ 2, 1, sum) / apply(count$count, 1, sum)
    #test2 <- count$count
    #for (j in 1:length(count$count))
    #  test2[,j] <- count$count[count$root, j] - test2[,j]
    test2 <- rep(as.vector(unlist(count$count[count$root,])), each = nrow(count$count)) - count$count
    test2 <- apply(test2 ^ 2, 1, sum) / apply(test2, 1, sum)
    test.stats <- test1 + test2
    
    max.stat <- max(test.stats, na.rm = T)
    res$root <- which.max(test.stats)
    res$impurity <- 1 - max.stat / count$obs
    
    res$pred.sub <- count$labels[which.max(count$count[res$root,])]
    suptree <- count$count[count$root,] - count$count[res$root,]
    res$pred.sup <- count$labels[which.max(suptree)]
  }
  else if (eval == "entropy")
  {
    checked <- 0
    min.stat <- count$obs
    sum1 <- apply(count$count, 1, sum)
    prop <- count$count / sum1
    test1 <- apply(prop * log(prop), 1, sum, na.rm = T)
    test2 <- rep(as.vector(unlist(count$count[count$root,])), each = nrow(count$count)) - count$count
    sum2 <- apply(test2, 1, sum)
    prop <- test2 / sum2
    test2 <- apply(prop * log(prop), 1, sum, na.rm = T)
    test.stats <- -(test1 * sum1 + test2 * sum2)
    
    min.stat <- min(test.stats, na.rm = T)
    res$root <- which.min(test.stats)
    res$impurity <- max.stat / count$obs
    
    res$pred.sub <- count$labels[which.max(count$count[res$root,])]
    suptree <- count$count[count$root,] - count$count[res$root,]
    res$pred.sup <- count$labels[which.max(suptree)]
  }
  return(res)
}

## predict
## DFS를 이용하면 서브트리 판정 문제를 간단하게 해결할 수 있다.
## 새로운 데이터를 고려하지 않을 때, build는 point 기반 tree이다.
## 새로운 데이터는 seg와 tp로 구성된 data.frame이며,
## 새로운 데이터를 고려할 때 build는 seg-tp 기반 tree이다.
## treeparty.predict.ada를 Rcpp 기반으로 재작성하기 위해
## treeparty.predict 또한 Rcpp 기반으로 작성해야 한다.
treeparty.predict <- function(build, split, newdata = NA, index = NA)
{
  if (!all(is.na(newdata)) && !all(is.na(index)))
    stop("[newdata] and [index] cannot exist together.")
  if (all(is.na(index)))
    index <- 1:length(build$marks)
  if (!is.na(split$pred))
    return(rep(split$pred, length(index)))
  if (!all(is.na(newdata)))
    stop("[newdata] is not supported yet.")
  
  #return(treeparty_predict(build, split, index))
  
  res <- factor(rep(NA, length(index)), levels(build$marks))
  on.subtree1 <- build$invorder[index] >= build$invorder[split$root]
  on.subtree2 <- build$invorder[index] < build$invorder[split$root] + build$subtree.size[split$root]
  on.subtree <- on.subtree1 & on.subtree2
  if (!is.na(split$pred.sub))
    res[on.subtree == TRUE] <- split$pred.sub
  else
    res[on.subtree == TRUE] <- treeparty.predict(build, split$subtree, NA, index[on.subtree == TRUE])
  if (!is.na(split$pred.sup))
    res[on.subtree == FALSE] <- split$pred.sup
  else
    res[on.subtree == FALSE] <- treeparty.predict(build, split$supertree, NA, index[on.subtree == FALSE])
  return(res)
}

treeparty.predict.ada <- function(build, split, limit = split$iter)
{
  if (class(split) != "treeparty.adaboost")
    stop("[split] should be 'treeparty.adaboost' object.")
  index <- 1:length(build$marks)
  labels <- levels(build$marks)
  preds <- data.frame(matrix(0, length(build$marks), length(labels)))
  colnames(preds) <- labels
  for (i in 1:limit)
  {
    pred <- treeparty.predict(build, split$stumps[[i]])
    for (j in 1:length(labels))
      preds[j] <- preds[j] + split$say[i] * (pred == labels[j])
  }
  res <- factor(levels(build$marks)[1], levels(build$marks))
  res <- rep(res, length(build$marks))
  for (j in 2:length(labels))
  {
    better <- preds[1] < preds[j]
    res[better] <- labels[j]
    preds[better, 1] <- preds[better, j]
  }
  return(res)
}

## 결정 트리
treeparty.decision <- function(build, weight = rep(1, length(build$marks)), minbucket = 30, depth = 3, eval = "gini", significance = 0.01, index = 1:length(build$marks))
{
  counted <- treeparty.count(build, weight)
  if (depth <= 0)
    minbucket <- sum(weight) + 1
  split <- treeparty.split(counted, minbucket, eval, significance)
  if (!is.na(split$pred))
    return(split)
  pred.sub <- split$pred.sub
  pred.sup <- split$pred.sup
  split$pred.sub <- TRUE
  split$pred.sup <- FALSE
  
  # prediction <- treeparty.predict(build, split)
  on.subtree1 <- build$invorder[index] >= build$invorder[split$root]
  on.subtree2 <- build$invorder[index] < build$invorder[split$root] + build$subtree.size[split$root]
  on.subtree <- on.subtree1 & on.subtree2
  on.suptree <- !on.subtree
  # end of prediction
  
  weight.subtree <- weight
  weight.subtree[index[on.suptree]] <- 0
  weight.suptree <- weight
  weight.suptree[index[on.subtree]] <- 0
  
  #split$pred.sub <- pred.sub
  #split$pred.sup <- pred.sup
  split$pred.sub <- NA
  split$pred.sup <- NA
  split$subtree <- treeparty.decision(build, weight.subtree, minbucket, depth - 1, eval, significance, index[on.subtree])
  split$supertree <- treeparty.decision(build, weight.suptree, minbucket, depth - 1, eval, significance, index[on.suptree])
  return(split)
}

print.treeparty.split <- function(split, depth = 0)
{
  if (depth == 0)
  {
    cat("Decision tree on linear network\n\n")
    cat("Roots of subtrees are determined by observations.\n\n")
    cat("root")
    if (!is.na(split$pred))
      cat(" -", levels(split$pred)[split$pred])
    else
    {
      if (split$eval == "chisq")
        cat(" (p = ", split$p.value, ")", sep = "")
      else if (split$eval == "accuracy")
        cat(" (acc = ", split$accuracy, ")", sep = "")
      else if (split$eval == "gini")
        cat(" (gini = ", split$impurity, ")", sep = "")
      len <- nchar(as.character(split$root))
      cat("\n  On subtree #", split$root, sep = "")
      if (len == 1)
        cat(' ')
      if (is.na(split$pred.sub))
        print.treeparty.split(split$subtree, depth + 1)
      else
        cat(" -", levels(split$pred.sub)[split$pred.sub])
      cat("\n  Out of subtree")
      if (len > 2)
        cat(rep(' ', len - 2))
      if (is.na(split$pred.sup))
        print.treeparty.split(split$supertree, depth + 1)
      else
        cat(" -", levels(split$pred.sup)[split$pred.sup])
    }
    cat('\n')
  }
  else
  {
    if (!is.na(split$pred))
      cat(" -", levels(split$pred)[split$pred])
    else
    {
      if (split$eval == "chisq")
        cat(" (p = ", split$p.value, ")", sep = "")
      else if (split$eval == "accuracy")
        cat(" (acc = ", split$accuracy, ")", sep = "")
      else if (split$eval == "gini")
        cat(" (gini = ", split$impurity, ")", sep = "")
      len <- nchar(as.character(split$root))
      cat("\n", rep("  ", depth + 1), "On subtree #", split$root, sep = "")
      if (len == 1)
        cat(' ')
      if (is.na(split$pred.sub))
        print.treeparty.split(split$subtree, depth + 1)
      else
        cat(" -", levels(split$pred.sub)[split$pred.sub])
      cat("\n", rep("  ", depth + 1), "Out of subtree", sep = "")
      if (len > 2)
        cat(rep(' ', len - 2))
      if (is.na(split$pred.sup))
        print.treeparty.split(split$supertree, depth + 1)
      else
        cat(" -", levels(split$pred.sup)[split$pred.sup])
    }
  }
}

## Hastie, T., Rosset, S., Zhu, J., & Zou, H. (2009). Multi-class adaboost. Statistics and its Interface, 2(3), 349-360.
treeparty.adaboost <- function(build, iter = 100, minbucket = 1, depth = 3, eval = "accuracy", weight = rep(1, length(build$marks)), verbose = FALSE)
{
  if (class(build) != "treeparty.build")
    stop("[build] should be a 'treeparty.build' object.")
  res <- list()
  res$iter <- iter
  res$say <- rep(0, iter)
  res$stumps <- vector("list", length = iter)
  res$classes <- length(levels(build$marks))
  for (i in 1:iter)
  {
    res$stumps[[i]] <- treeparty.decision(build, weight, minbucket, depth, eval = eval)
    pred <- treeparty.predict(build, res$stumps[[i]])
    error <- sum(weight[build$marks != pred]) / sum(weight)
    if (verbose == TRUE && i %% 10 == 0)
      cat("[", i, "] ", error, "\n", sep = '')
    res$say[i] <- log(1 - error) - log(error) + log(res$classes - 1)
    weight <- weight * exp(res$say[i] * (build$marks != pred))
    weight <- weight / sum(weight) * length(build$marks)
  }
  class(res) <- "treeparty.adaboost"
  res$eval <- eval
  return(res)
}


sourceCpp("treeparty.cpp")
