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