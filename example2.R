## sample 2 : maze
v <- ppp(x = rep(1:10, 10), y = rep(1:10, each = 10), c(0, 11), c(0, 11))
edge <- matrix(c(1:99,
                 2, 12, 2, 5, 15, 5, 6, 7, 19, 9, 
                 12, 13, 14, 4, 25, 17, 18, 19, 29, 10, 
                 11, 21, 24, 34, 24, 16, 37, 38, 28, 20, 
                 32, 42, 32, 44, 25, 26, 36, 48, 40, 50, 
                 31, 52, 33, 43, 35, 45, 37, 49, 39, 60, 
                 41, 62, 54, 44, 56, 57, 47, 57, 49, 70, 
                 62, 63, 64, 65, 55, 65, 77, 67, 59, 80, 
                 61, 82, 72, 73, 74, 66, 87, 68, 78, 79, 
                 82, 83, 93, 94, 75, 85, 86, 89, 90, 100, 
                 81, 91, 94, 95, 96, 97, 98, 88, 98),
               nrow = 99, ncol = 2)
maze <- linnet(v, edges = edge)
mainroute <- c(1, 2, 4, 5,
               12, 13, 14, 15, 16, 17, 18, 19,
               24, 25, 26, 28, 29,
               32, 33, 34, 36, 37, 38, 39, 40,
               42, 43, 44, 47, 48, 49, 50,
               52, 55, 56, 57, 60,
               62, 63, 64, 65, 67, 68, 70,
               72, 73, 74, 75, 77, 78, 79, 80,
               82, 83, 85, 86, 87, 88, 89, 90,
               93, 94, 95, 96, 97, 98)
junction <- c(2, 5, 12, 15, 25, 24, 44, 32, 62,
              65, 57, 37, 19, 49, 82, 94, 98)

set.seed(88)
X <- runiflpp(n = 396, L = maze, nsim = 1)
X$data <- cbind(X$data, marks = sapply(X$data$seg, function(x) ifelse(x %in% mainroute, "main", "sub")))
fromjunction <- sapply(X$data$seg, function(x) x %in% junction)
tojunction <- sapply(X$domain$to[X$data$seg], function(x) x %in% junction)
mutation <- (1 - X$data$tp) * fromjunction + X$data$tp * tojunction
X$data$marks[runif(396) < mutation] <- "junction"

X$data$marks <- factor(X$data$marks, c("main", "sub", "junction"))
X$ctype <- c(X$ctype, as.factor("mark"))

png("maze.png", 480, 400)
par(mar = c(0, 2, 0, 0))
par(cex = 2)
plot.lpp(X, main = "")
dev.off()

seed.min <- 1L
seed.max <- 100L
iter <- 100L
correct.treeparty <- matrix(0L, seed.max, iter)
correct.Euclidean <- matrix(0L, seed.max, iter)

starttime <- Sys.time()
for (seed in seed.min:seed.max)
{
  set.seed(seed)
  X <- runiflpp(n = 396, L = maze, nsim = 1)
  X$data <- cbind(X$data, marks = sapply(X$data$seg, function(x) ifelse(x %in% mainroute, "main", "sub")))
  fromjunction <- sapply(X$data$seg, function(x) x %in% junction)
  tojunction <- sapply(X$domain$to[X$data$seg], function(x) x %in% junction)
  mutation <- (1 - X$data$tp) * fromjunction + X$data$tp * tojunction
  X$data$marks[runif(396) < mutation] <- "junction"
  X$data$marks <- as.factor(X$data$marks)
  X$ctype <- c(X$ctype, as.factor("mark"))
  plot(X, main = paste("Seed", as.character(seed)))
  
  ## Leave-one-out cross validation for treeparty.adaboost
  
  visited <- treeparty.visit(X)
  built <- treeparty.build(X, visited)
  labels <- levels(built$marks)
  for (i in 1L:length(built$marks))
  {
    weight <- rep(1, length(built$marks))
    weight[i] <- 0
    adapted <- treeparty.adaboost(built, depth = 3, iter = iter, weight = weight, eval = "accuracy")
    
    ## iteration에 따른 예측값의 변화를 모두 확인하기 위해
    ## treeparty.predict.ada의 코드를 긁어왔다.
    preds <- rep(0, length(labels))
    for (it in 1L:iter)
    {
      pred <- treeparty.predict(built, adapted$stumps[[it]], index = i)
      pred <- as.integer(pred)
      preds[pred] <- preds[pred] + adapted$say[it]
      pred <- which.max(preds)
      if (pred == as.integer(built$marks[i]))
        correct.treeparty[seed, it] <- correct.treeparty[seed, it] + 1L
    }
    #cat(i, "/", 396, '\n')
    if (i < 396 || seed < seed.max)
      print((Sys.time() - starttime) * (396 - i + (seed.max - seed) * 396) / (i + (seed - seed.min) * 396))
  }
  write.csv(correct.treeparty, "treeparty.accuracy.csv")
  
  y <- X$data$marks
  x <- data.frame(x = X$data$x, y = X$data$y)
  
  labels <- levels(y)
  for (i in 1:length(y))
  {
    adapted2 <- ada(x[-i,], y[-i], n_rounds = iter, verbose = FALSE, progress = FALSE, split = "accuracy")
    preds <- rep(0, length(labels))
    for (it in 1L:iter)
    {
      pred <- stats::predict(adapted2$trees[[it]], x[i,], type = "class")
      pred <- as.integer(pred)
      preds[pred] <- preds[pred] + adapted2$alphas[it]
      pred <- which.max(preds)
      if (pred == as.integer(y[i]))
        correct.Euclidean[seed, it] <- correct.Euclidean[seed, it] + 1L
    }
    if (i < 396 || seed < seed.max)
      print((Sys.time() - starttime) * (396 - i + (seed.max - seed) * 396) / (i + (seed - seed.min) * 396))
  }
  write.csv(correct.Euclidean, "Euclidean.accuracy.csv")
}

seed.min <- 1L
seed.max <- 100L
iter <- 20L
correct.knn <- matrix(0, seed.max, iter)
starttime <- Sys.time()
for (seed in seed.min:seed.max)
{
  set.seed(seed)
  X <- runiflpp(n = 396, L = maze, nsim = 1)
  X$data <- cbind(X$data, marks = sapply(X$data$seg, function(x) ifelse(x %in% mainroute, "main", "sub")))
  fromjunction <- sapply(X$data$seg, function(x) x %in% junction)
  tojunction <- sapply(X$domain$to[X$data$seg], function(x) x %in% junction)
  mutation <- (1 - X$data$tp) * fromjunction + X$data$tp * tojunction
  X$data$marks[runif(396) < mutation] <- "junction"
  X$data$marks <- as.factor(X$data$marks)
  X$ctype <- c(X$ctype, as.factor("mark"))
  plot(X, main = paste("Seed", as.character(seed)))
  
  ## Leave-one-out cross validation for treeparty.adaboost
  
  
  y <- X$data$marks
  x <- data.frame(x = X$data$x, y = X$data$y)
  
  labels <- levels(y)
  for (i in 1:length(y))
  {
    for (it in 1L:iter)
      if (knn(x[-i,], x[i,], y[-i], it) == y[i])
        correct.knn[seed, it] <- correct.knn[seed, it] + 1
    if ((i < 396 || seed < seed.max) && i %% 99 == 0)
      print((Sys.time() - starttime) * (396 - i + (seed.max - seed) * 396) / (i + (seed - seed.min) * 396))
  }
  write.csv(correct.knn, "knn.csv")
}