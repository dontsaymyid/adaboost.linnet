# sample 3 : 금강
lst <- readRDS("Geum_linnet.RDS") # df$Geum$data에는 toc가 들어 있다
Geum <- lst$Geum
rm(lst)

# Good: 2급수 이상
# Bad : 3급수 이하
X <- Geum
Geumdf <- as.data.frame(Geum$data)
month <- 1 # 2022년 1월
good <- 2
X$data <- Geum$data[, c(1:4, month + 4)]

X$data[,5] <- factor(ifelse(Geumdf[, month + 4] <= good, "Good", "Bad"), c("Good", "Bad"))
X$ctype <- Geum$ctype[1:5]
X$data <- as.data.frame(X$data) ## We cannot use hyperframe since
colnames(X$data)[5] <- "marks"  ## this removes a column of X$data
rm(list = c("month", "good"))
table(X$data$marks)
png("Geum.png", 640, 480)
par(mar = c(0, 2, 0, 0))
par(cex = 3)
plot(X, cols = c("blue", "red"), main = "")
dev.off()
plot(X, cols = c("blue", "red"), main = "") # 2022년 6월 금강 총유기탄소

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
correct.gini <- correct

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
  adapted2 <- ada(x[-i,], y[-i], n_rounds = iter, verbose = FALSE, progress = FALSE, tree_depth = 3, split = "gini")
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
correct2.gini <- correct2

correct3 <- rep(0L, iter)
y <- X$data$marks
x <- data.frame(x = X$data$x, y = X$data$y)

labels <- levels(y)
for (i in 1:length(y))
{
  for (it in seq(1L, iter, 2))
    if (knn(x[-i,], x[i,], y[-i], it) == y[i])
      correct3[it] <- correct3[it] + 1
}



X.wrong <- X
X.wrong2 <- X
X.wrong$data <- X$data[!correct.last,]
X.wrong2$data <- X$data[!correct2.last,]
png("geum.pred.png", 1280, 480)
par(mfrow = c(1, 2))
par(mar = c(0, 2, 0, 0))
par(cex = 3)
plot(X, cols = "pink", main = "")
plot(X.wrong, add = T, cols = "red")
plot(X, cols = "pink", main = "")
plot(X.wrong2, add = T, cols = "red")
dev.off()

png("correct.E.png", width = 640, height = 480)
X.correct2 <- X
X.correct2$data <- X$data[correct2.last]
plot(X, main = "Euclidean adaboost")
plot(X.correct2, add = TRUE, cols = "pink")
dev.off()

png("accuracy.geum.png", 1000, 576)
par(mfrow = c(2, 4))
plot(correct.accuracy / length(built$marks) * 100, ylim = c(70, 90), xlab = "iteration", ylab = "accuracy (%)", main = "(a) Proposed: MR")
abline(h = 80, col = "lightgray")
plot(correct.gini / length(built$marks) * 100, ylim = c(70, 90), xlab = "iteration", ylab = "accuracy (%)", main = "(b) Proposed: gini")
abline(h = 80, col = "lightgray")
plot(correct.entropy / length(built$marks) * 100, ylim = c(70, 90), xlab = "iteration", ylab = "accuracy (%)", main = "(c) Proposed: entropy")
abline(h = 80, col = "lightgray")
plot.new()
plot(correct2.accuracy / length(built$marks) * 100, ylim = c(70, 90), xlab = "iteration", ylab = "accuracy (%)", main = "(d) Comparison: MR")
abline(h = 80, col = "lightgray")
plot(correct2.gini / length(built$marks) * 100, ylim = c(70, 90), xlab = "iteration", ylab = "accuracy (%)", main = "(e) Comparison: gini")
abline(h = 80, col = "lightgray")
plot(correct2.entropy / length(built$marks) * 100, ylim = c(70, 90), xlab = "iteration", ylab = "accuracy (%)", main = "(f) Comparison: entropy")
abline(h = 80, col = "lightgray")
plot(correct3 / length(built$marks) * 100, ylim = c(70, 90), xlab = "neighbors", ylab = "accuracy (%)", main = "(g) Comparison: knn")
abline(h = 80, col = "lightgray")
dev.off()

## 5-fold cross validation

set.seed(0)

tb <- table(X$data$marks)
folds <- rep(1, 129)
for (i in 1:2)
  folds[X$data$marks == names(tb[i])] <- as.integer((order(runif(tb[i])) - 1) / tb[i] * 5) + 1

iter <- 20L

visited <- treeparty.visit(X)
built <- treeparty.build(X, visited)
labels <- levels(built$marks)

correct <- rep(0L, iter)
for (i in 1L:5L)
{
  weight <- ifelse(folds == i, 0, 1)
  adapted <- treeparty.adaboost(built, depth = 3, iter = iter, weight = weight, eval = "accuracy")
  
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
    correct[it] <- correct[it] + sum(pred == as.integer(built$marks[where]))
  }
  cat(i, "/", 5, '\n')
}
correct.accuracy <- correct

correct2 <- rep(0L, iter)
labels <- levels(y)
for (i in 1:5)
{
  where <- which(folds == i)
  adapted2 <- ada(x[-where,], y[-where], n_rounds = iter, verbose = F, progress = FALSE, tree_depth = 3, split = "accuracy")
  preds <- matrix(0, length(where), length(labels))
  for (it in 1L:iter)
  {
    pred <- stats::predict(adapted2$trees[[it]], x[where,], type = "class")
    pred <- as.integer(pred)
    for (f in 1:nrow(preds))
      preds[f, pred[f]] <- preds[f, pred[f]] + adapted2$alphas[it]
    pred <- apply(preds, 1, which.max)
    correct2[it] <- correct2[it] + sum(pred == as.integer(y[where]))
  }
  cat(i, "/", 5, '\n')
}
correct2.accuracy <- correct2

correct3 <- rep(0L, iter)
y <- X$data$marks
x <- data.frame(x = X$data$x, y = X$data$y)

labels <- levels(y)
for (i in 1:5)
{
  where <- which(folds == i)
  for (it in seq(1L, iter, 2))
    correct3[it] <- correct3[it] + sum(knn(x[-where,], x[where,], y[-where], it) == y[where])
}

png("accuracy.geum5.png", 1000, 576)
par(mfrow = c(2, 4))
plot(correct.accuracy / length(built$marks) * 100, ylim = c(70, 90), xlab = "iteration", ylab = "accuracy (%)", main = "(a) Proposed: MR")
abline(h = 80, col = "lightgray")
plot(correct.gini / length(built$marks) * 100, ylim = c(70, 90), xlab = "iteration", ylab = "accuracy (%)", main = "(b) Proposed: gini")
abline(h = 80, col = "lightgray")
plot(correct.entropy / length(built$marks) * 100, ylim = c(70, 90), xlab = "iteration", ylab = "accuracy (%)", main = "(c) Proposed: entropy")
abline(h = 80, col = "lightgray")
plot.new()
plot(correct2.accuracy / length(built$marks) * 100, ylim = c(70, 90), xlab = "iteration", ylab = "accuracy (%)", main = "(d) Comparison: MR")
abline(h = 80, col = "lightgray")
plot(correct2.gini / length(built$marks) * 100, ylim = c(70, 90), xlab = "iteration", ylab = "accuracy (%)", main = "(e) Comparison: gini")
abline(h = 80, col = "lightgray")
plot(correct2.entropy / length(built$marks) * 100, ylim = c(70, 90), xlab = "iteration", ylab = "accuracy (%)", main = "(f) Comparison: entropy")
abline(h = 80, col = "lightgray")
plot(correct3 / length(built$marks) * 100, ylim = c(70, 90), xlab = "neighbors", ylab = "accuracy (%)", main = "(g) Comparison: knn")
abline(h = 80, col = "lightgray")
dev.off()



## 10-fold cross validation

set.seed(0)

tb <- table(X$data$marks)
folds <- rep(1, 129)
for (i in 1:2)
  folds[X$data$marks == names(tb[i])] <- as.integer((order(runif(tb[i])) - 1) / tb[i] * 10) + 1

iter <- 20L

visited <- treeparty.visit(X)
built <- treeparty.build(X, visited)
labels <- levels(built$marks)

correct <- rep(0L, iter)
for (i in 1L:10L)
{
  weight <- ifelse(folds == i, 0, 1)
  adapted <- treeparty.adaboost(built, depth = 3, iter = iter, weight = weight, eval = "entropy")
  
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
    correct[it] <- correct[it] + sum(pred == as.integer(built$marks[where]))
  }
  cat(i, "/", 10, '\n')
}
correct.entropy <- correct

labels <- levels(y)
correct2 <- rep(0L, iter)
for (i in 1:10)
{
  where <- which(folds == i)
  adapted2 <- ada(x[-where,], y[-where], n_rounds = iter, verbose = F, progress = FALSE, tree_depth = 3, split = "entropy")
  preds <- matrix(0, length(where), length(labels))
  for (it in 1L:iter)
  {
    pred <- stats::predict(adapted2$trees[[it]], x[where,], type = "class")
    pred <- as.integer(pred)
    for (f in 1:nrow(preds))
      preds[f, pred[f]] <- preds[f, pred[f]] + adapted2$alphas[it]
    pred <- apply(preds, 1, which.max)
    correct2[it] <- correct2[it] + sum(pred == as.integer(y[where]))
  }
  cat(i, "/", 10, '\n')
}
correct2.entropy <- correct2

correct3 <- rep(0L, iter)
y <- X$data$marks
x <- data.frame(x = X$data$x, y = X$data$y)

labels <- levels(y)
for (i in 1:10)
{
  where <- which(folds == i)
  for (it in seq(1L, iter, 2))
    correct3[it] <- correct3[it] + sum(knn(x[-where,], x[where,], y[-where], it) == y[where])
}

png("accuracy.geum10.png", 1000, 576)
par(mfrow = c(2, 4))
plot(correct.accuracy / length(built$marks) * 100, ylim = c(60, 90), xlab = "iteration", ylab = "accuracy (%)", main = "(a) Proposed: MR")
abline(h = 80, col = "lightgray")
plot(correct.gini / length(built$marks) * 100, ylim = c(60, 90), xlab = "iteration", ylab = "accuracy (%)", main = "(b) Proposed: gini")
abline(h = 80, col = "lightgray")
plot(correct.entropy / length(built$marks) * 100, ylim = c(60, 90), xlab = "iteration", ylab = "accuracy (%)", main = "(c) Proposed: entropy")
abline(h = 80, col = "lightgray")
plot.new()
plot(correct2.accuracy / length(built$marks) * 100, ylim = c(60, 90), xlab = "iteration", ylab = "accuracy (%)", main = "(d) Comparison: MR")
abline(h = 80, col = "lightgray")
plot(correct2.gini / length(built$marks) * 100, ylim = c(60, 90), xlab = "iteration", ylab = "accuracy (%)", main = "(e) Comparison: gini")
abline(h = 80, col = "lightgray")
plot(correct2.entropy / length(built$marks) * 100, ylim = c(60, 90), xlab = "iteration", ylab = "accuracy (%)", main = "(f) Comparison: entropy")
abline(h = 80, col = "lightgray")
plot(correct3 / length(built$marks) * 100, ylim = c(60, 90), xlab = "neighbors", ylab = "accuracy (%)", main = "(g) Comparison: knn")
abline(h = 80, col = "lightgray")
dev.off()

# Displaying the performances
res <- rbind(correct.accuracy[1:4 * 5], correct.gini[1:4 * 5], correct.entropy[1:4 * 5],
             correct2.accuracy[1:4 * 5], correct2.gini[1:4 * 5], correct2.entropy[1:4 * 5],
             correct3[c(5, 9, 15, 19)])
res <- cbind(res,
             c(max(correct.accuracy), max(correct.gini), max(correct.entropy),
               max(correct2.accuracy), max(correct2.gini), max(correct2.entropy),
               max(correct3))) / 1.29
round(res, 1)