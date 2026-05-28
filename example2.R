MUTATEMORE = F
JUNCTION1.5 = F
ROTATE45 = F


## sample 2 : maze
set.seed(18)
vx <- rep(0:9, 10)
vy <- rep(0:9, each = 10)
if (MUTATEMORE)
{
  vx <- vx + runif(100)
  vy <- vy + runif(100)
}
if (ROTATE45)
  v <- ppp(x = vx - vy, y = vx + vy, c(-10, 10), c(0, 20))
if (!ROTATE45)
  v <- ppp(vx, vy, c(0, 10), c(0, 10))
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

X <- runiflpp(n = 396, L = maze, nsim = 1)
X$data <- cbind(X$data, marks = sapply(X$data$seg, function(x) ifelse(x %in% mainroute, "main", "sub")))

mutation <- rep(0, 396)

for (j in junction)
{
  jx <- v$x[j]
  jy <- v$y[j]
  if (ROTATE45)
    for (i in 1:396)
      mutation[i] <- mutation[i] + max(0, 1 - sqrt(((jx - X$data$x[i]) ^ 2 + (jy - X$data$y[i]) ^ 2) / 2))
  else
    for (i in 1:396)
      mutation[i] <- mutation[i] + max(0, 1 - sqrt((jx - X$data$x[i]) ^ 2 + (jy - X$data$y[i]) ^ 2))
}

if (JUNCTION1.5)
  X$data$marks[runif(396) < mutation * 1.5] <- "junction"
if (!JUNCTION1.5)
  X$data$marks[runif(396) < mutation] <- "junction"

X$data$marks <- factor(X$data$marks, c("main", "sub", "junction"))
if (MUTATEMORE)
{
  noise <- sample(1:99, 5)
  X.noise <- X$data$seg %in% noise
  X$data$marks[X.noise] <- levels(X$data$marks)[as.integer(runif(sum(X.noise), 1, 4))]
}
X$ctype <- c(X$ctype, as.factor("mark"))
table(X$data$marks)

png("maze.png", 480, 400)
par(mar = c(0, 2, 0, 0))
par(cex = 2)
plot(X, main = "")
dev.off()

seed.min <- 1L
seed.max <- 100L
iter <- 100L

evals <- c("accuracy", "gini", "entropy")
## correct.*  : for proposed
## correct2.* : for comparison

for (eval in evals)
{
correct.treeparty <- matrix(0L, seed.max, iter)
correct.Euclidean <- matrix(0L, seed.max, iter)

starttime <- Sys.time()
for (seed in seed.min:seed.max)
{
  set.seed(seed)
  vx <- rep(0:9, 10)
  vy <- rep(0:9, each = 10)
  if (MUTATEMORE)
  {
    vx <- vx + runif(100)
    vy <- vy + runif(100)
  }
  if (ROTATE45)
    v <- ppp(x = vx - vy, y = vx + vy, c(-10, 10), c(0, 20))
  if (!ROTATE45)
    v <- ppp(vx, vy, c(0, 10), c(0, 10))
  maze <- linnet(v, edges = edge)
  X <- runiflpp(n = 396, L = maze, nsim = 1)
  X$data <- cbind(X$data, marks = sapply(X$data$seg, function(x) ifelse(x %in% mainroute, "main", "sub")))
  mutation <- rep(0, 396)
  for (j in junction)
  {
    jx <- v$x[j]
    jy <- v$y[j]
    for (i in 1:396)
      mutation[i] <- mutation[i] + max(0, 1 - sqrt((jx - X$data$x[i]) ^ 2 + (jy - X$data$y[i]) ^ 2) / 1)
  }
  X$data$marks[runif(396) < mutation] <- "junction"
  X$data$marks <- factor(X$data$marks, c("main", "sub", "junction"))
  if (MUTATEMORE)
  {
    noise <- sample(1:99, 5)
    X.noise <- X$data$seg %in% noise
      X$data$marks[X.noise] <- levels(X$data$marks)[as.integer(runif(sum(X.noise), 1, 4))]
  }
  X$ctype <- c(X$ctype, as.factor("mark"))
  tb <- table(X$data$marks)
  plot(X, main = paste("Seed", as.character(seed)), cex = 2)
  
  ## Leave-one-out cross validation for treeparty.adaboost
  
  visited <- treeparty.visit(X)
  built <- treeparty.build(X, visited)
  labels <- levels(built$marks)
  for (i in 1L:length(built$marks))
  {
    weight <- rep(1, length(built$marks))
    weight[i] <- 0
    adapted <- treeparty.adaboost(built, depth = 3, iter = iter, weight = weight, eval = eval)
    
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
  write.csv(correct.treeparty, paste("treeparty.", eval, ".csv", sep = ""))
  next
  
  y <- X$data$marks
  x <- data.frame(x = X$data$x, y = X$data$y)
  
  labels <- levels(y)
  for (i in 1:length(y))
  {
    adapted2 <- ada(x[-i,], y[-i], n_rounds = iter, verbose = F, progress = F, split = eval)
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
  write.csv(correct.Euclidean, paste("Euclidean.", eval, ".csv", sep = ""))
}
}

seed.min <- 1L
seed.max <- 100L
iter <- 20L
correct.knn <- matrix(0, seed.max, iter)
starttime <- Sys.time()
for (seed in seed.min:seed.max)
{
  set.seed(seed)
  vx <- rep(0:9, 10)
  vy <- rep(0:9, each = 10)
  if (MUTATEMORE)
  {
    vx <- vx + runif(100)
    vy <- vy + runif(100)
  }
  if (ROTATE45)
    v <- ppp(x = vx - vy, y = vx + vy, c(-10, 10), c(0, 20))
  if (!ROTATE45)
    v <- ppp(vx, vy, c(0, 10), c(0, 10))
  maze <- linnet(v, edges = edge)
  X <- runiflpp(n = 396, L = maze, nsim = 1)
  X$data <- cbind(X$data, marks = sapply(X$data$seg, function(x) ifelse(x %in% mainroute, "main", "sub")))
  mutation <- rep(0, 396)
  for (j in junction)
  {
    jx <- v$x[j]
    jy <- v$y[j]
    for (i in 1:396)
      mutation[i] <- mutation[i] + max(0, 1 - sqrt((jx - X$data$x[i]) ^ 2 + (jy - X$data$y[i]) ^ 2) / 1)
  }
  X$data$marks[runif(396) < mutation] <- "junction"
  X$data$marks <- factor(X$data$marks, c("main", "sub", "junction"))
  if (MUTATEMORE)
  {
    noise <- sample(1:99, 5)
    X.noise <- X$data$seg %in% noise
    X$data$marks[X.noise] <- levels(X$data$marks)[as.integer(runif(sum(X.noise), 1, 4))]
  }
  X$ctype <- c(X$ctype, as.factor("mark"))
  tb <- table(X$data$marks)
  folds <- rep(1, 396)
  for (i in 1:3)
    folds[X$data$marks == names(tb[i])] <- as.integer((order(runif(tb[i])) - 1) / tb[i] * 5) + 1
  plot(X, main = paste("Seed", as.character(seed)), cex = 2)
  
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

## 5-fold cross validation

seed.min <- 1L
seed.max <- 100L
iter <- 100L

for (eval in evals)
{
correct.treeparty <- matrix(0L, seed.max, iter)
correct.Euclidean <- matrix(0L, seed.max, iter)

for (seed in seed.min:seed.max)
{
  set.seed(seed)
  vx <- rep(0:9, 10)
  vy <- rep(0:9, each = 10)
  if (MUTATEMORE)
  {
    vx <- vx + runif(100)
    vy <- vy + runif(100)
  }
  if (ROTATE45)
    v <- ppp(x = vx - vy, y = vx + vy, c(-10, 10), c(0, 20))
  if (!ROTATE45)
    v <- ppp(vx, vy, c(0, 10), c(0, 10))
  maze <- linnet(v, edges = edge)
  X <- runiflpp(n = 396, L = maze, nsim = 1)
  X$data <- cbind(X$data, marks = sapply(X$data$seg, function(x) ifelse(x %in% mainroute, "main", "sub")))
  mutation <- rep(0, 396)
  for (j in junction)
  {
    jx <- v$x[j]
    jy <- v$y[j]
    for (i in 1:396)
      mutation[i] <- mutation[i] + max(0, 1 - sqrt((jx - X$data$x[i]) ^ 2 + (jy - X$data$y[i]) ^ 2) / 1)
  }
  X$data$marks[runif(396) < mutation] <- "junction"
  X$data$marks <- factor(X$data$marks, c("main", "sub", "junction"))
  if (MUTATEMORE)
  {
    noise <- sample(1:99, 5)
    X.noise <- X$data$seg %in% noise
    X$data$marks[X.noise] <- levels(X$data$marks)[as.integer(runif(sum(X.noise), 1, 4))]
  }
  X$ctype <- c(X$ctype, as.factor("mark"))
  tb <- table(X$data$marks)
  folds <- rep(1, 396)
  for (i in 1:3)
    folds[X$data$marks == names(tb[i])] <- as.integer((order(runif(tb[i])) - 1) / tb[i] * 5) + 1
  plot(X, main = paste("Seed", as.character(seed)), cex = 2)
  
  visited <- treeparty.visit(X)
  built <- treeparty.build(X, visited)
  labels <- levels(built$marks)
  for (i in 1L:5L)
  {
    weight <- ifelse(folds == i, 0, 1)
    adapted <- treeparty.adaboost(built, depth = 3, iter = iter, weight = weight, eval = eval)
    
    ## iteration에 따른 예측값의 변화를 모두 확인하기 위해
    ## treeparty.predict.ada의 코드를 긁어왔다.
    where <- which(folds == i)
    preds <- matrix(0, length(where), length(labels))
    for (it in 1L:iter)
    {
      pred <- treeparty.predict(built, adapted$stumps[[it]], index = where)
      pred <- as.integer(pred)
      for (f in 1:nrow(preds))
        preds[f, pred[f]] <- preds[f, pred[f]] + adapted$say[it]
      pred <- apply(preds, 1, which.max)
      correct.treeparty[seed, it] <- correct.treeparty[seed, it] + sum(pred == as.integer(built$marks[where]))
    }
    cat(i, "/", 5, '\n')
  }
  write.csv(correct.treeparty, paste("treeparty5.", eval, ".csv", sep = ""))
  y <- X$data$marks
  x <- data.frame(x = X$data$x, y = X$data$y)
  
  labels <- levels(y)
  for (i in 1L:5L)
  {
    where <- which(folds == i)
    adapted2 <- ada(x[-where,], y[-where], n_rounds = iter, verbose = F, progress = F, split = eval)
    preds <- matrix(0, length(where), length(labels))
    for (it in 1L:iter)
    {
      pred <- stats::predict(adapted2$trees[[it]], x[where,], type = "class")
      pred <- as.integer(pred)
      for (f in 1:nrow(preds))
        preds[f, pred[f]] <- preds[f, pred[f]] + adapted2$alphas[it]
      pred <- apply(preds, 1, which.max)
      correct.Euclidean[seed, it] <- correct.Euclidean[seed, it] + sum(pred == as.integer(y[where]))
    }
    cat(i, "/", 5, '\n')
  }
  write.csv(correct.Euclidean, paste("Euclidean5.", eval, ".csv", sep = ""))
}
}
seed.min <- 1L
seed.max <- 100L
iter <- 20L
correct.knn <- matrix(0, seed.max, iter)
starttime <- Sys.time()
for (seed in seed.min:seed.max)
{
  set.seed(seed)
  vx <- rep(0:9, 10)
  vy <- rep(0:9, each = 10)
  if (MUTATEMORE)
  {
    vx <- vx + runif(100)
    vy <- vy + runif(100)
  }
  if (ROTATE45)
    v <- ppp(x = vx - vy, y = vx + vy, c(-10, 10), c(0, 20))
  if (!ROTATE45)
    v <- ppp(vx, vy, c(0, 10), c(0, 10))
  maze <- linnet(v, edges = edge)
  X <- runiflpp(n = 396, L = maze, nsim = 1)
  X$data <- cbind(X$data, marks = sapply(X$data$seg, function(x) ifelse(x %in% mainroute, "main", "sub")))
  mutation <- rep(0, 396)
  for (j in junction)
  {
    jx <- v$x[j]
    jy <- v$y[j]
    for (i in 1:396)
      mutation[i] <- mutation[i] + max(0, 1 - sqrt((jx - X$data$x[i]) ^ 2 + (jy - X$data$y[i]) ^ 2) / 1)
  }
  X$data$marks[runif(396) < mutation] <- "junction"
  X$data$marks <- factor(X$data$marks, c("main", "sub", "junction"))
  if (MUTATEMORE)
  {
    noise <- sample(1:99, 5)
    X.noise <- X$data$seg %in% noise
    X$data$marks[X.noise] <- levels(X$data$marks)[as.integer(runif(sum(X.noise), 1, 4))]
  }
  X$ctype <- c(X$ctype, as.factor("mark"))
  tb <- table(X$data$marks)
  folds <- rep(1, 396)
  for (i in 1:3)
    folds[X$data$marks == names(tb[i])] <- as.integer((order(runif(tb[i])) - 1) / tb[i] * 5) + 1
  plot(X, main = paste("Seed", as.character(seed)), cex = 2)
  
  ## Leave-one-out cross validation for treeparty.adaboost
  
  y <- X$data$marks
  x <- data.frame(x = X$data$x, y = X$data$y)
  
  labels <- levels(y)
  for (i in 1:5)
  {
    for (it in 1L:iter)
      correct.knn[seed, it] <- correct.knn[seed, it] + sum(knn(x[folds != i,], x[folds == i,], y[folds != i], it) == y[folds == i])
    cat(i, "/", 5, '\n')
  }
  write.csv(correct.knn, "knn5.csv")
}


## 10-fold cross validation

seed.min <- 1L
seed.max <- 100L
iter <- 100L

for (eval in evals)
{
correct.treeparty <- matrix(0L, seed.max, iter)
correct.Euclidean <- matrix(0L, seed.max, iter)

for (seed in seed.min:seed.max)
{
  set.seed(seed)
  vx <- rep(0:9, 10)
  vy <- rep(0:9, each = 10)
  if (MUTATEMORE)
    
    set.seed(seed)
  vx <- rep(0:9, 10)
  vy <- rep(0:9, each = 10)
  if (MUTATEMORE)
  {
    vx <- vx + runif(100)
    vy <- vy + runif(100)
  }
  if (ROTATE45)
    v <- ppp(x = vx - vy, y = vx + vy, c(-10, 10), c(0, 20))
  if (!ROTATE45)
    v <- ppp(vx, vy, c(0, 10), c(0, 10))
  maze <- linnet(v, edges = edge)
  X <- runiflpp(n = 396, L = maze, nsim = 1)
  X$data <- cbind(X$data, marks = sapply(X$data$seg, function(x) ifelse(x %in% mainroute, "main", "sub")))
  mutation <- rep(0, 396)
  for (j in junction)
  {
    jx <- v$x[j]
    jy <- v$y[j]
    for (i in 1:396)
      mutation[i] <- mutation[i] + max(0, 1 - sqrt((jx - X$data$x[i]) ^ 2 + (jy - X$data$y[i]) ^ 2) / 1)
  }
  X$data$marks[runif(396) < mutation] <- "junction"
  X$data$marks <- factor(X$data$marks, c("main", "sub", "junction"))
  if (MUTATEMORE)
  {
    noise <- sample(1:99, 5)
    X.noise <- X$data$seg %in% noise
    X$data$marks[X.noise] <- levels(X$data$marks)[as.integer(runif(sum(X.noise), 1, 4))]
  }
  X$ctype <- c(X$ctype, as.factor("mark"))
  tb <- table(X$data$marks)
  folds <- rep(1, 396)
  for (i in 1:3)
    folds[X$data$marks == names(tb[i])] <- as.integer((order(runif(tb[i])) - 1) / tb[i] * 10) + 1
  plot(X, main = paste("Seed", as.character(seed)), cex = 2)
  
  visited <- treeparty.visit(X)
  built <- treeparty.build(X, visited)
  labels <- levels(built$marks)
  for (i in 1L:10L)
  {
    weight <- ifelse(folds == i, 0, 1)
    adapted <- treeparty.adaboost(built, depth = 3, iter = iter, weight = weight, eval = eval)
    
    ## iteration에 따른 예측값의 변화를 모두 확인하기 위해
    ## treeparty.predict.ada의 코드를 긁어왔다.
    where <- which(folds == i)
    preds <- matrix(0, length(where), length(labels))
    for (it in 1L:iter)
    {
      pred <- treeparty.predict(built, adapted$stumps[[it]], index = where)
      pred <- as.integer(pred)
      for (f in 1:nrow(preds))
        preds[f, pred[f]] <- preds[f, pred[f]] + adapted$say[it]
      pred <- apply(preds, 1, which.max)
      correct.treeparty[seed, it] <- correct.treeparty[seed, it] + sum(pred == as.integer(built$marks[where]))
    }
    cat(i, "/", 10, '\n')
  }
  next
  write.csv(correct.treeparty, paste("treeparty10.5.", eval, ".csv", sep = ""))
  y <- X$data$marks
  x <- data.frame(x = X$data$x, y = X$data$y)
  
  labels <- levels(y)
  for (i in 1L:10L)
  {
    where <- which(folds == i)
    adapted2 <- ada(x[-where,], y[-where], n_rounds = iter, verbose = F, progress = F, split = eval)
    preds <- matrix(0, length(where), length(labels))
    for (it in 1L:iter)
    {
      pred <- stats::predict(adapted2$trees[[it]], x[where,], type = "class")
      pred <- as.integer(pred)
      for (f in 1:nrow(preds))
        preds[f, pred[f]] <- preds[f, pred[f]] + adapted2$alphas[it]
      pred <- apply(preds, 1, which.max)
      correct.Euclidean[seed, it] <- correct.Euclidean[seed, it] + sum(pred == as.integer(y[where]))
    }
    cat(i, "/", 10, '\n')
  }
  write.csv(correct.Euclidean, paste("Euclidean10.", eval, ".csv", sep = ""))
}
}


## 18번 시드 재현

seed <- 18
iter <- 100
set.seed(seed)
vx <- rep(0:9, 10)
vy <- rep(0:9, each = 10)
if (MUTATEMORE)
{
  vx <- vx + runif(100)
  vy <- vy + runif(100)
}
if (ROTATE45)
  v <- ppp(x = vx - vy, y = vx + vy, c(-10, 10), c(0, 20))
if (!ROTATE45)
  v <- ppp(vx, vy, c(0, 10), c(0, 10))
maze <- linnet(v, edges = edge)
X <- runiflpp(n = 396, L = maze, nsim = 1)
X$data <- cbind(X$data, marks = sapply(X$data$seg, function(x) ifelse(x %in% mainroute, "main", "sub")))
mutation <- rep(0, 396)
for (j in junction)
{
  jx <- v$x[j]
  jy <- v$y[j]
  for (i in 1:396)
    mutation[i] <- mutation[i] + max(0, 1 - sqrt((jx - X$data$x[i]) ^ 2 + (jy - X$data$y[i]) ^ 2) / 1)
}
X$data$marks[runif(396) < mutation] <- "junction"
X$data$marks <- factor(X$data$marks, c("main", "sub", "junction"))
if (MUTATEMORE)
{
  noise <- sample(1:99, 5)
  X.noise <- X$data$seg %in% noise
  X$data$marks[X.noise] <- levels(X$data$marks)[as.integer(runif(sum(X.noise), 1, 4))]
}
X$ctype <- c(X$ctype, as.factor("mark"))
plot(X, main = paste("Seed", as.character(seed)))

visited <- treeparty.visit(X)
built <- treeparty.build(X, visited)
labels <- levels(built$marks)
pred.treeparty <- rep(0, length(built$marks))
starttime <- Sys.time()
for (i in 1L:length(built$marks))
{
  weight <- rep(1, length(built$marks))
  weight[i] <- 0
  adapted <- treeparty.adaboost(built, depth = 3, iter = iter, weight = weight, eval = "accuracy")
  
  preds <- rep(0, 3)
  ## iteration에 따른 예측값의 변화를 모두 확인하기 위해
  ## treeparty.predict.ada의 코드를 긁어왔다.
  for (it in 1L:iter)
  {
    pred <- treeparty.predict(built, adapted$stumps[[it]], index = i)
    pred <- as.integer(pred)
    preds[pred] <- preds[pred] + adapted$say[it]
  }
  pred.treeparty[i] <- which.max(preds)
  
  #cat(i, "/", 396, '\n')
  if (i < 396)
    print((Sys.time() - starttime) * (396 - i) / i)
}

pred.treeparty <- factor(pred.treeparty, levels = 1:3)
levels(pred.treeparty) <- levels(built$marks)

Xpred <- X
Xpred$data$marks <- pred.treeparty
X.wrong <- Xpred
X.wrong$data <- Xpred$data[built$marks != pred.treeparty]
png("maze.treeparty.accuracy.png", 480, 400)
par(mar = c(0, 0, 5, 0))
par(cex = 1)
plot(Xpred, cex = 3, cols = "pink", main = "(a) Proposed: MR")
plot(X.wrong, add = T, cols = "black", cex = 3)
dev.off()

y <- X$data$marks
x <- data.frame(x = X$data$x, y = X$data$y)

labels <- levels(y)
pred.Euclidean <- rep(0, length(built$marks))
starttime <- Sys.time()
for (i in 1:length(y))
{
  adapted2 <- ada(x[-i,], y[-i], n_rounds = iter, verbose = FALSE, progress = FALSE)
  preds <- rep(0, 3)
  for (it in 1L:iter)
  {
    pred <- stats::predict(adapted2$trees[[it]], x[i,], type = "class")
    pred <- as.integer(pred)
    preds[pred] <- preds[pred] + adapted2$alphas[it]
  }
  pred.Euclidean[i] <- which.max(preds)
  if (i < 396)
    print((Sys.time() - starttime) * (396 - i) / i)
}

pred.Euclidean <- factor(pred.Euclidean, levels = 1:3)
levels(pred.Euclidean) <- levels(built$marks)

Xpred <- X
Xpred$data$marks <- pred.Euclidean
X.wrong <- Xpred
X.wrong$data <- Xpred$data[built$marks != pred.Euclidean]
png("maze.Euclidean.accuracy.png", 480, 400)
par(mar = c(0, 0, 5, 0))
par(cex = 1)
plot(Xpred, cex = 3, cols = "pink", main = "(b) Comparison: MR")
plot(X.wrong, add = T, cols = "black", cex = 3)
dev.off()

sum.tb <- c(0L, 0L, 0L)

## counting the sample observations
for (seed in seed.min:seed.max)
{
  set.seed(seed)
  vx <- rep(0:9, 10)
  vy <- rep(0:9, each = 10)
  if (MUTATEMORE)
  {
    vx <- vx + runif(100)
    vy <- vy + runif(100)
  }
  if (ROTATE45)
    v <- ppp(x = vx - vy, y = vx + vy, c(-10, 10), c(0, 20))
  if (!ROTATE45)
    v <- ppp(vx, vy, c(0, 10), c(0, 10))
  maze <- linnet(v, edges = edge)
  X <- runiflpp(n = 396, L = maze, nsim = 1)
  X$data <- cbind(X$data, marks = sapply(X$data$seg, function(x) ifelse(x %in% mainroute, "main", "sub")))
  mutation <- rep(0, 396)
  for (j in junction)
  {
    jx <- v$x[j]
    jy <- v$y[j]
    for (i in 1:396)
      mutation[i] <- mutation[i] + max(0, 1 - sqrt((jx - X$data$x[i]) ^ 2 + (jy - X$data$y[i]) ^ 2) / 1)
  }
  X$data$marks[runif(396) < mutation] <- "junction"
  X$data$marks <- factor(X$data$marks, c("main", "sub", "junction"))
  if (MUTATEMORE)
  {
    noise <- sample(1:99, 5)
    X.noise <- X$data$seg %in% noise
    X$data$marks[X.noise] <- levels(X$data$marks)[as.integer(runif(sum(X.noise), 1, 4))]
  }
  X$ctype <- c(X$ctype, as.factor("mark"))
  tb <- table(X$data$marks)
  sum.tb <- sum.tb + tb
}
print(sum.tb)


set.seed(seed <- 18)
vx <- rep(0:9, 10)
vy <- rep(0:9, each = 10)
if (MUTATEMORE)
{
  vx <- vx + runif(100)
  vy <- vy + runif(100)
}
if (ROTATE45)
  v <- ppp(x = vx - vy, y = vx + vy, c(-10, 10), c(0, 20))
if (!ROTATE45)
  v <- ppp(vx, vy, c(0, 10), c(0, 10))
maze <- linnet(v, edges = edge)
X <- runiflpp(n = 396, L = maze, nsim = 1)
X$data <- cbind(X$data, marks = sapply(X$data$seg, function(x) ifelse(x %in% mainroute, "main", "sub")))
mutation <- rep(0, 396)
for (j in junction)
{
  jx <- v$x[j]
  jy <- v$y[j]
  for (i in 1:396)
    mutation[i] <- mutation[i] + max(0, 1 - sqrt((jx - X$data$x[i]) ^ 2 + (jy - X$data$y[i]) ^ 2) / 1)
}
X$data$marks[runif(396) < mutation] <- "junction"
X$data$marks <- factor(X$data$marks, c("main", "sub", "junction"))
if (MUTATEMORE)
{
  noise <- sample(1:99, 5)
  X.noise <- X$data$seg %in% noise
  X$data$marks[X.noise] <- levels(X$data$marks)[as.integer(runif(sum(X.noise), 1, 4))]
}
X$ctype <- c(X$ctype, as.factor("mark"))
tb <- table(X$data$marks)
plot(X, main = paste("Seed", as.character(seed)), cex = 2)
