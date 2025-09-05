setwd("C:\\Users\\donts\\OneDrive\\Desktop\\강의실\\1S\\연구과제")
#setwd("C:\\Users\\A\\Desktop\\강의실\\1S\\연구과제")
rm(list = ls())

source("treeparty.R")
#load(".RData")

## Comparison methods

ada <- function(X, y, tree_depth = 3, n_rounds = 100, verbose = FALSE,
                control = NULL, progress = FALSE, gini = TRUE)
{
  labels <- length(levels(y))
  if (!is.data.frame(X))
    stop("[X] must be a data frame.")
  if (length(X$marks) > 0)
    stop("[X] should not contain a column [marks].")
  df <- cbind(X, data.frame(marks = y))
  if (is.null(control))
    control = rpart::rpart.control(minsplit = 0, minbucket = 1, 
                                   cp = -1, maxcompete = 0, maxsurrogate = 0, usesurrogate = 0, 
                                   xval = 0, maxdepth = tree_depth)
  else if (control$maxdepth != tree_depth)
    warning(paste("tree_depth set to: ", control$maxdepth))
  n = nrow(X)
  w = rep(1 / n, n)
  trees = list()
  alphas = list()
  count <- matrix(0, n, labels)
  if (progress)
    acc <- rep(0, n_rounds)
  for (i in seq(n_rounds))
  {
    if (gini)
      tree = rpart::rpart(marks ~ ., data = df, weights = w, 
                          method = "class", control = control,
                          x = FALSE, y = FALSE, model = FALSE)
    else
      tree = rpart::rpart(marks ~ ., data = df, weights = w, 
                          method = alist, control = control,
                          x = FALSE, y = FALSE, model = FALSE)
    tree$where = NULL
    tree$call = NULL
    tree$cptable = NULL
    tree$functions = NULL
    tree$control = NULL
    tree$variable.importance = NULL
    tree$parms = NULL
    pred = stats::predict(tree, X, type = "class")
    #pred = xpred.rpart(tree2)
    if (progress)
    {
      yhat <- levels(y)[apply(count, 1, which.max)]
      acc[i] <- sum(yhat == y) / n
    }
    e = sum(w * (pred != y))
    if (abs(e) < 1e-08)
    {
      if (i == 1)
      {
        trees[[i]] = tree
        alphas[[i]] = 1
        terms = tree$terms
      }
      n_rounds = i
      break
    }
    alpha = log((1 - e)/e) + log(labels - 1)
    for (j in seq(labels))
      count[as.integer(pred) == j, j] <- count[as.integer(pred) == j, j] + alpha
    w = w * exp(alpha * (pred != y))
    w = w / sum(w)
    #if (i == 1)
    terms = tree$terms
    #else
    #  tree$terms = NULL
    trees[[i]] = tree
    alphas[[i]] = alpha / 2
    if (verbose & (i %% 10 == 0))
      cat("Iteration: ", i, "\n")
  }
  if (progress)
    out = list(alphas = unlist(alphas), trees = trees,
               tree_depth = tree_depth, terms = terms, acc = acc)
  else
    out = list(alphas = unlist(alphas), trees = trees,
               tree_depth = tree_depth, terms = terms)
  class(out) = "adaboost"
  
  #yhat = stats::predict(out, as.matrix(X))
  out$yhat <- levels(y)[apply(count, 1, which.max)]
  out$confusion_matrix = table(y, out$yhat)
  out
}

# for knn
library(class)

data("dendrite")
plot(dendrite)

X <- dendrite

## Predicting itself using treeparty.decision and treeparty.adaboost

visited <- treeparty.visit(X)
built <- treeparty.build(X, visited)
counted <- treeparty.count(built)
#counted.R <- treeparty.count.R(built)

split <- treeparty.split(counted)
treeparty.predict(built, split)

{
  y <- X$data$marks
  x <- data.frame(x = X$data$x, y = X$data$y)
  decided <- treeparty.decision(built, depth = 3, eval = "gini", minbucket = 10)
  decided2 <- ada(x, y, 3, 1, verbose = TRUE, progress = TRUE)
  treeparty.predict(built, decided)
  factor(decided2$yhat, levels(X$data$marks))
  decided2$trees[[1]]
}

decided <- treeparty.decision(built, depth = 3, eval = "gini", minbucket = 10)
decided
pred <- treeparty.predict(built, decided)
sum(pred == built$marks) / length(built$marks)

visited <- treeparty.visit(X)
built <- treeparty.build(X, visited)
adapted <- treeparty.adaboost(built, depth = 3, iter = 100, eval = "gini", verbose = TRUE)
pred <- treeparty.predict.ada(built, adapted)
sum(pred == built$marks) / length(built$marks)

adapted2 <- ada(x, y, 3, 100, verbose = TRUE, progress = TRUE)
pred2 <- factor(adapted2$yhat, levels(X$data$marks))
sum(pred2 == built$marks) / length(built$marks)

correct = pred == built$marks

acc <- rep(0, adapted$iter)
for (i in 1:adapted$iter)
{
  acc[i] <- sum(treeparty.predict.ada(built, adapted, i) == built$marks) / length(built$marks)
  cat(max(acc), '/', i, '\n')
}
plot(acc, xlab = "iteration", ylab = "accuracy",
     ylim = c(0, 1),
     main = "Accuracy by adaboost on linear network")

acc2 <- rep(0, adapted$iter)
for (i in 1:adapted$iter)
{
  adapted2 <- ada(x, y, 3, i, verbose = F, progress = TRUE)
  pred2 <- factor(adapted2$yhat, levels(X$data$marks))
  acc2[i] <- sum(pred2 == built$marks) / length(built$marks)
  cat(max(acc), '/', i, '\n')
}
plot(acc2, xlab = "iteration", ylab = "accuracy",
     ylim = c(0, 1),
     main = "Accuracy by Euclidean adaboost")



## sample 1 : t
v <- ppp(x = c(-1, 0, 1, 0), y = c(0, 0, 0, 1), c(-1.01, 1.01), c(-0.01, 1.01))
edge <- matrix(c(2, 2, 2, 1, 3, 4), 3, 2)
tshape <- linnet(v, edges = edge)

set.seed(1557)
X <- runiflpp(n = 300, L = tshape, nsim = 1)

marks <- matrix(0, 3, 4)
marks[,1] <- 0:2
for (i in 2:4)
  marks[,i] <- marks[,i - 1] + as.integer(runif(3, 1, 3))
marks <- marks %% 3
marks <- as.vector(marks)

mark <- X$data$seg * 4 + 1
X$data <- cbind(X$data, marks = sapply(X$data$tp, function(x) ifelse(sin(x * pi * 12) * 5 > runif(1, -1, 1), "A", "B")))
X$data$marks[X$data$tp > 3/24 & X$data$tp < 7/24 & X$data$marks == "A"] <- "C"
X$data$marks[X$data$tp > 9/24 & X$data$tp < 13/24 & X$data$marks == "B"] <- "C"
X$data$marks[X$data$tp > 15/24 & X$data$tp < 19/24 & X$data$marks == "A"] <- "C"
X$data$marks[X$data$tp > 21/24 & X$data$marks == "B"] <- "C"
#10% 확률로 돌연변이 생성
noise <- runif(nrow(X$data), 0, 20)
X$data$marks[noise < 3] <- "C"
X$data$marks[noise < 2] <- "B"
X$data$marks[noise < 1] <- "A"

X$data$marks <- as.factor(X$data$marks)
X$ctype <- c(X$ctype, as.factor("mark"))


png("tshape.png", 960, 480)
par(mar = c(0, 0, 0, 0))
par(cex = 2)
plot(X, main = "", cex = 1.5)
dev.off()

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

# sample 3 : 금강
lst <- readRDS("Geum_linnet.RDS") # df$Geum$data에는 toc가 들어 있다
Geum <- lst$Geum
rm(lst)

# Good: 2급수 이상
# Bad : 3급수 이하
X <- Geum
Geumdf <- as.data.frame(Geum$data)
X$data <- Geum$data[, c(1:4, 10)] # 2022년 6월
X$data <- Geum$data[, c(1:4, 14)] # 2022년 10월
#X$data <- Geum$data[, c(1:4, 22)] # 2023년 6월

X$data[,5] <- factor(ifelse(Geumdf[, 10] <= 4, "Good", "Bad"), c("Good", "Bad"))
X$ctype <- Geum$ctype[1:5]
X$data <- as.data.frame(X$data) ## We cannot use hyperframe since
colnames(X$data)[5] <- "marks"  ## this removes a column of X$data
png("Geum.png", 640, 480)
par(mar = c(0, 2, 0, 0))
par(cex = 3)
plot(X, cols = c("blue", "red"), main = "")
dev.off()
plot(X, cols = c("blue", "red"), main = "") # 2022년 6월 금강 총유기탄소


# sample 4 : dendrite
X <- dendrite

## self-prediction for identity test
iter <- 100L
visited <- treeparty.visit(X)
built <- treeparty.build(X, visited)
adapted <- treeparty.adaboost(built, depth = 3, iter = iter, verbose = TRUE, eval = "gini")
pred <- treeparty.predict.ada(built, adapted)
sum(pred == built$marks) / nrow(X$data)

y <- X$data$marks
x <- data.frame(x = X$data$x, y = X$data$y)
adapted2 <- ada(x, y, n_rounds = iter, verbose = FALSE, progress = FALSE, tree_depth = 3)
pred2 <- factor(adapted2$yhat, levels(built$marks))

sum(pred2 == built$marks) / nrow(X$data)

sum(pred != pred2)

correct <- pred == built$marks

Xpred <- X
Xpred$data$marks <- pred
X.wrong <- Xpred
X.wrong$data <- X$data[!correct]
png("linear.T.png", 960, 192)
par(mar = c(0, 0, 0, 0))
par(cex = 2)
plot(Xpred, cex = 1.5, cols = "pink", main = "")
plot(X.wrong, add = T, cols = "black", cex = 1.5)
dev.off()

acc <- rep(0, adapted$iter)
for (i in 1:adapted$iter)
  acc[i] <- sum(treeparty.predict.ada(built, adapted, i) == built$marks) / nrow(X$data)
plot(acc, xlab = "iteration", ylab = "accuracy",
     ylim = c(0, 1),
     main = "Accuracy by adaboost on linear network")

acc2 <- rep(0, adapted$iter)
for (i in 1:adapted$iter)
{
  adapted2 <- ada(x, y, n_rounds = i, verbose = FALSE, progress = FALSE)
  pred2 <- factor(adapted2$yhat, levels(built$marks))
  acc2[i] <- sum(pred2 == built$marks) / length(built$marks)
}
plot(acc2, xlab = "iteration", ylab = "accuracy",
     ylim = c(0, 1),
     main = "Accuracy by Euclidean adaboost")


## Leave-one-out cross validation for treeparty.adaboost

iter <- 20L

visited <- treeparty.visit(X)
built <- treeparty.build(X, visited)
correct <- rep(0L, iter)
correct.last <- rep(FALSE, length(built$marks))
labels <- levels(built$marks)
starttime <- Sys.time()
for (i in 1L:length(built$marks))
{
  weight <- rep(1, length(built$marks))
  weight[i] <- 0
  adapted <- treeparty.adaboost(built, depth = 3, iter = iter, weight = weight, eval = "entropy")
  
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
    {
      correct[it] <- correct[it] + 1L
      if (it == iter)
        correct.last[i] <- TRUE
    }
  }
  #cat(max(correct), "/", i, '\n')
  if (i < length(built$marks))
    print((Sys.time() - starttime) * (length(built$marks) - i) / i)
}

X.wrong <- X
X.wrong$data <- X$data[!correct.last]
png("geum.T2.png", 640, 480)
par(mar = c(0, 2, 0, 0))
par(cex = 3)
plot(X, cols = "pink", main = "")
plot(X.wrong, add = T, cols = "red")
dev.off()


y <- X$data$marks
x <- data.frame(x = X$data$x, y = X$data$y)

correct2 <- rep(0L, iter)
correct2.last <- rep(FALSE, length(y))
labels <- levels(y)
starttime <- Sys.time()
for (i in 1:length(y))
{
  adapted2 <- ada(x[-i,], y[-i], n_rounds = iter, verbose = FALSE, progress = FALSE, tree_depth = 3, gini = TRUE)
  preds <- rep(0, length(labels))
  for (it in 1L:iter)
  {
    pred <- stats::predict(adapted2$trees[[it]], x[i,], type = "class")
    pred <- as.integer(pred)
    preds[pred] <- preds[pred] + adapted2$alphas[it]
    pred <- which.max(preds)
    if (pred == as.integer(y[i]))
    {
      #correct[it] <- correct[it] - 1L
      correct2[it] <- correct2[it] + 1L
      if (it == iter)
        correct2.last[i] <- TRUE
    }
  }
  #cat(max(correct2), "/", i, '\n')
  if (i < length(y))
    print((Sys.time() - starttime) * (length(y) - i) / i)
}

correct3 <- rep(0L, iter)
y <- X$data$marks
x <- data.frame(x = X$data$x, y = X$data$y)

labels <- levels(y)
for (i in 1:length(y))
{
  for (it in 1L:iter)
    if (knn(x[-i,], x[i,], y[-i], it) == y[i])
      correct3[it] <- correct3[it] + 1
}



X.wrong <- X
X.wrong$data <- X$data[!correct2.last]
png("geum.E.png", 640, 480)
par(mar = c(0, 2, 0, 0))
par(cex = 3)
plot(X, cols = "pink", main = "")
plot(X.wrong, add = T, cols = "red")
dev.off()

png("correct.E.png", width = 640, height = 480)
X.correct2 <- X
X.correct2$data <- X$data[correct2.last]
plot(X, main = "Euclidean adaboost")
plot(X.correct2, add = TRUE, cols = "pink")
dev.off()

png("accuracy.geum2.png", width = 640, height = 363)
par(mfrow = c(1, 2))
plot(correct / length(built$marks) * 100, ylim = c(70, 90), xlab = "iteration", ylab = "accuracy (%)", main = "(a) Proposed method")
abline(h = 80, col = "lightgray")
plot(correct2 / length(y) * 100, ylim = c(70, 90), xlab = "iteration", ylab = "accuracy", main = "(b) Contrast method")
abline(h = 80, col = "lightgray")
dev.off()

plot(correct / length(built$marks) * 100, ylim = c(70, 90), xlab = "iteration", ylab = "accuracy (%)", main = "(a) Proposed algorithm")
abline(h = 80, col = "lightgray")
plot(correct2 / length(y) * 100, ylim = c(70, 90), xlab = "iteration", ylab = "accuracy", main = "(b) Euclidean adaboost")
abline(h = 80, col = "lightgray")


# 
seed.max <- 100L
iter <- 100L
correct.treeparty2 <- matrix(0L, seed.max, iter)
correct.Euclidean <- matrix(0L, seed.max, iter)

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
    adapted <- treeparty.adaboost(built, depth = 3, iter = iter, weight = weight, eval = "gini")
    
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
  write.csv(correct.treeparty, "treeparty.entropy.csv")
  next
  
  y <- X$data$marks
  x <- data.frame(x = X$data$x, y = X$data$y)
  
  labels <- levels(y)
  for (i in 1:length(y))
  {
    adapted2 <- ada(x[-i,], y[-i], n_rounds = iter, verbose = FALSE, progress = FALSE)
    preds <- rep(0, length(labels))
    for (it in 1L:iter)
    {
      pred <- stats::predict(adapted2$trees[[it]], x[i,], type = "class")
      pred <- as.integer(pred)
      preds[pred] <- preds[pred] + adapted2$alphas[it]
      pred <- which.max(preds)
      if (pred == as.integer(y[i]))
        correct.Euclidean2[seed, it] <- correct.Euclidean2[seed, it] + 1L
    }
    cat(i + 396, "/", 792, '\n')
  }
  write.csv(correct.Euclidean, "Euclidean.csv")
}

difference <- correct.treeparty[1:seed,] - correct.Euclidean[1:seed,]

png("boxplot.gini.png", 629, 480)
boxplot(difference, xlab = "iteration", ylab = "difference", axes = F,
        main = "")
axis(1, c(-10, 1:11 * 10))
axis(2, -5:5 * 50)
axis(3, c(-10, 110))
axis(4, c(-200, 200))
abline(h = 0)
dev.off()
better <- rep(0, 100)
for (i in 1:seed)
  better <- better + (difference[i,] > 0) + (difference[i,] == 0) * 0.5
png("better.gini.png", 629, 480)
plot(better / seed, xlab = "iteration", ylab = "probability",
     axes = F, main = "")
axis(1, c(-10, 1:11 * 10))
axis(2, 0:10 * 0.1)
axis(3, c(-10, 110))
axis(4, c(-100, 200))
abline(h = 0.5)
abline(h = 0.58, lty = 2)
abline(h = 0.1)
dev.off()

## Remove a seed
correct.treeparty2[seed,]
correct.treeparty2[seed,] <- 0
correct.Euclidean2[seed,] <- 0
seed.min <- seed

seed.min <- seed + 1

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


# seed, iter
correct.Euclidean <- as.matrix(read.csv("Euclidean.csv"))[,-1]
correct.treeparty <- as.matrix(read.csv("treeparty.accuracy.csv"))[,-1]
avg.Euclidean <- rep(0, 100)
avg.treeparty <- rep(0, 100)
avg.treeparty2 <- rep(0, 100)

for (seed in 1:100)
  avg.Euclidean <- avg.Euclidean + as.vector(correct.Euclidean[seed,])
avg.treeparty <- rep(0, 100)
for (seed in 1:100)
  avg.treeparty <- avg.treeparty + as.vector(correct.treeparty[seed,])
for (seed in 1:100)
  avg.treeparty2 <- avg.treeparty2 + as.vector(correct.treeparty2[seed,])

avg.Euclidean <- avg.Euclidean / 100
avg.treeparty <- avg.treeparty / 100
avg.treeparty2 <- avg.treeparty2 / 100
plot(avg.Euclidean, type = 'l', ylim = c(0, 396))
points(avg.treeparty)
abline(h = 396 - 44, lty = 2)

for (iter in 1:100)
  print(mean(correct.treeparty[,iter] - correct.Euclidean[,iter]) / sd(correct.treeparty[,iter] - correct.Euclidean[,iter]))

for (iter in c(5, 10, 20, 30, 50, 75, 100))
  print(sum(correct.treeparty[,iter] > correct.Euclidean[,iter]) + sum(correct.treeparty[,iter] == correct.Euclidean[,iter]) / 2)

png("accuracy.maze.png", width = 750, height = 300)
par(mfrow = c(1, 3))
plot(avg.treeparty2 / 396 * 100, type = 'l', ylim = c(65, 85), xlab = "iteration", ylab = "accuracy (%)", main = "(a) Proposed: MR")
abline(h = 80, col = "lightgray")
plot(avg.treeparty / 396 * 100, type = 'l', ylim = c(65, 85), xlab = "iteration", ylab = "accuracy (%)", main = "(b) Proposed: gini")
abline(h = 80, col = "lightgray")
plot(avg.Euclidean / 396 * 100, type = 'l', ylim = c(65, 85), xlab = "iteration", ylab = "accuracy (%)", main = "(c) Comparison: gini")
abline(h = 80, col = "lightgray")
dev.off()

pred2 <- stats::predict(adapted2$trees[[it]], x[i,], type = "class")


## 88번 시드 재현

seed <- 88
iter <- 67
set.seed(88)
X <- runiflpp(n = 396, L = maze, nsim = 1)
X$data <- cbind(X$data, marks = sapply(X$data$seg, function(x) ifelse(x %in% mainroute, "main", "sub")))
fromjunction <- sapply(X$data$seg, function(x) x %in% junction)
tojunction <- sapply(X$domain$to[X$data$seg], function(x) x %in% junction)
mutation <- (1 - X$data$tp) * fromjunction + X$data$tp * tojunction
X$data$marks[runif(396) < mutation] <- "junction"
X$data$marks <- as.factor(X$data$marks)
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
png("maze.Euclidean.png", 480, 400)
par(mar = c(0, 0, 5, 0))
par(cex = 1)
plot(Xpred, cex = 3, cols = "pink", main = "(c) Comparison: gini")
plot(X.wrong, add = T, cols = "black", cex = 3)
dev.off()

