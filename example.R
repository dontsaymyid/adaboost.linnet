source("treeparty.R")
#load(".RData")

## Comparison methods

ada <- function(X, y, tree_depth = 3, n_rounds = 100, verbose = FALSE,
                control = NULL, progress = FALSE, split = "gini")
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
    if (split == "gini")
      tree = rpart::rpart(marks ~ ., data = df, weights = w, 
                          method = "class", control = control,
                          x = FALSE, y = FALSE, model = FALSE)
    else if (split == "accuracy")
    {
      tree = rpart::rpart(marks ~ ., data = df, weights = w, 
                          method = alist, control = control,
                          x = FALSE, y = FALSE, model = FALSE)
      attr(tree, "ylevels") <- levels(y)    
    }
    else if (split == "entropy")
    {
      tree = rpart::rpart(marks ~ ., data = df, weights = w, 
                          method = elist, control = control,
                          x = FALSE, y = FALSE, model = FALSE)
      attr(tree, "ylevels") <- levels(y)    
    }
    else
      stop("invalid split")
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
    if (verbose)
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

adapted2 <- ada(x[-i,], y[-i], n_rounds = iter, verbose = T, progress = F, split = "gini")


# for knn
library(class)