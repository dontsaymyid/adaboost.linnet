library(spatstat)
library(spatstat.linnet)
source('treeparty_SC.R', chdir = TRUE)
source('comparison.split.R', chdir = TRUE)
source('treeparty.loop1_SC.R', chdir = TRUE)
library(class)

library(pROC) 


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

get_metrics <- function(actual, predicted, scores = NULL) {
  target_levels <- c("main", "sub", "junction")
  
  # 팩터 변환 (레벨 누락 시에도 안전하게 테이블 생성)
  actual <- factor(actual, levels = target_levels)
  predicted <- factor(predicted, levels = target_levels)
  
  cm <- table(Actual = actual, Predicted = predicted)
  cm_mat <- as.matrix(cm)
  
  # 각 클래스별 정밀도/재현율 계산 (에러 방지를 위해 + 1e-10 추가)
  precision <- diag(cm_mat) / (colSums(cm_mat) + 1e-10)
  recall <- diag(cm_mat) / (rowSums(cm_mat) + 1e-10)
  
  # F1 (Macro-averaged) 및 Balanced Accuracy
  f1_macro <- mean(2 * (precision * recall) / (precision + recall + 1e-10), na.rm = TRUE)
  balanced_acc <- mean(recall, na.rm = TRUE)
  
  return(list(F1 = f1_macro, B_Acc = balanced_acc, AUC = NA))
}

calc_impurity <- function(counts, eval_type) {
  total <- sum(counts)
  if(total == 0) return(0)
  p <- counts / total
  p <- p[p > 0]
  if (eval_type == "gini") {
    return(1 - sum(p^2))
  } else if (eval_type == "entropy") {
    return(-sum(p * log(p)))
  } else { # accuracy (MR)
    return(1 - max(p))
  }
}


# 사이클 정점들의 리스트(cycle_nodes)와 C++ 결과(count_matrix)가 주어졌을 때
condense_cycle_weights <- function(count_matrix, cycle_nodes, parent_array) {
  
  n_cycle <- length(cycle_nodes)
  condensed_weights <- matrix(0, nrow = n_cycle, ncol = ncol(count_matrix))
  
  for (i in 1:n_cycle) {
    v <- cycle_nodes[i]
    total_w <- count_matrix[v, ]
    
    # v의 자식들 중, '사이클 위에 있는 자식' 찾기
    # parent_array를 뒤져서 부모가 v인 노드들 중 cycle_nodes에 포함된 노드
    cycle_children <- which(parent_array == v)
    cycle_children <- cycle_children[cycle_children %in% cycle_nodes]
    
    for (child in cycle_children) {
      total_w <- total_w - count_matrix[child, ]
    }
    
    condensed_weights[i, ] <- total_w
  }
  
  return(condensed_weights) 
}

# ---------------------------------------------------------
# Data generation
# ---------------------------------------------------------
v <- ppp(x = c(-3, 3, 3, -3, 1, -1), 
         y = c(-3, -3, 3, 3, -4, 4), 
         window = owin(c(-4, 4), c(-4, 4)))


edge <- matrix(c(1, 2,
                 2, 3,
                 3, 4,
                 4, 1,
                 1, 5,
                 3, 6), nrow = 6, ncol = 2, byrow = TRUE)

loop_net <- linnet(v, edges = edge)

res_final_list <- vector("list", 100)

for (seed in 1:100) {
  cat("\n========================================\n")
  cat("Running Simulation for Seed:", seed, "\n")
  cat("========================================\n")
  
  # -------------------------------------------------------------------------
  # Data preparation
  # -------------------------------------------------------------------------
  set.seed(seed) 
  
  v <- ppp(x = c(-0.5, 0.5, 0.5, -0.5, 1, -1), 
           y = c(-0.5, -0.5, 0.5, 0.5, -4, 4), 
           window = owin(c(-4, 4), c(-4, 4)))
  edge <- matrix(c(1, 2, 2, 3, 3, 4, 4, 1, 1, 5, 3, 6), nrow = 6, ncol = 2, byrow = TRUE)
  loop_net <- linnet(v, edges = edge)
  mainroute <- c(6, 1, 2)
  
  # Generate n
  n <- 150
  X <- runiflpp(n = n, L = loop_net)
  
  # Generate network
  X$data$marks[X$data$seg %in% c(1, 2, 5)] <- "main"
  X$data$marks[X$data$seg %in% c(3, 4, 6)] <- "sub"
  
  # 소수 접합점 생성 
  # 노드 1: (-3, -3), 노드 3: (3, 3)
  dist_to_node1 <- sqrt((X$data$x - (-3))^2 + (X$data$y - (-3))^2)
  dist_to_node3 <- sqrt((X$data$x - 3)^2 + (X$data$y - 3)^2)
  
  # 반경 1.0 이내의 점들만 junction으로 덮어씌움 (전체 데이터의 약 10~15% 정도)
  junction_radius <- 1.0
  X$data$marks[dist_to_node1 < junction_radius | dist_to_node3 < junction_radius] <- "junction"
  
  # add 15% noise
  noise_idx <- sample(1:n, round(0.15 * n))
  X$data$marks[noise_idx] <- sample(c("main", "sub", "junction"), length(noise_idx), replace = TRUE)
  X$data$marks <- factor(X$data$marks, levels = c("main", "sub", "junction"))
  
  # visualization
  plot(X, cols = c("blue", "red", "green"), main = paste("Micro-Loop Network - Seed", seed), cex = 1.2)
  legend("topleft", legend=c("main", "sub", "junction"), col=c("blue", "red", "green"), pch=1)
  
  # indexing
  df_X <- as.data.frame(X)
  df_X$orig_id <- 1:nrow(df_X)
  is_cycle <- (X$domain$from %in% 1:4) & (X$domain$to %in% 1:4)
  ordered_seg_ids <- c(which(X$domain$from == 1 & X$domain$to == 2),
                       which(X$domain$from == 2 & X$domain$to == 3),
                       which(X$domain$from == 3 & X$domain$to == 4),
                       which(X$domain$from == 4 & X$domain$to == 1))
  
  cycle_data <- do.call(rbind, lapply(ordered_seg_ids, function(s) {
    df <- df_X[df_X$seg == s, ]; return(df[order(df$tp), ])
  }))
  n_cycle <- nrow(cycle_data)
  cycle_data$cycle_k <- 1:n_cycle
  df_X$cycle_k <- NA
  df_X$cycle_k[match(cycle_data$orig_id, df_X$orig_id)] <- cycle_data$cycle_k
  df_X$cycle_k[df_X$seg %in% which(!is_cycle & (X$domain$from == 1 | X$domain$to == 1))] <- 1
  df_X$cycle_k[df_X$seg %in% which(!is_cycle & (X$domain$from == 3 | X$domain$to == 3))] <- n_cycle
  
  # 데이터를 cycle_k 기준으로 정렬한 후, 다시 df_X에 매핑해 준다
  df_X <- df_X[order(df_X$cycle_k), ] 
  
  n_obs <- nrow(df_X); y <- df_X$marks; x <- data.frame(x = df_X$x, y = df_X$y)
  labels <- levels(y); pos_col_idx <- which(labels == "main")
  folds <- sample(rep(1:10, length.out = n_obs))
  
  # -------------------------------------------------------------------------
  # 10-fold CV (Proposed vs Euclidean Comparison)
  # -------------------------------------------------------------------------
  x_center <- mean(v$x[is_cycle])
  y_center <- mean(v$y[is_cycle])
  v.condensed <- ppp(x = ifelse(is_cycle, x_center, v$x), 
                     y = ifelse(is_cycle, y_center, v$y), 
                     window = v$window)
  
  loop_net.condensed <- linnet(v.condensed, edges = edge)
  
  v.core <- which.max(is_cycle)
  X.condensed <- lpp(L = loop_net.condensed)
  X.condensed$data <- X$data
  
  # contract the loop
  X.condensed$domain$from[is_cycle[X.condensed$domain$from] & !is_cycle] <- v.core
  X.condensed$domain$to[is_cycle[X.condensed$domain$to] & !is_cycle] <- v.core
  X.condensed$domain$lines$ends <- X.condensed$domain$lines$ends[-v.core,]
  X.condensed$domain$lines$n <- X.condensed$domain$lines$n - 1
  X.condensed$domain$from <- X.condensed$domain$from[-v.core]
  X.condensed$domain$to <- X.condensed$domain$to[-v.core]
  
  # find seg-tp of the specific point
  {
    if (any(X.condensed$domain$from == v.core))
    {
      seg.core <- which(X.condensed$domain$from == v.core)[1]
      tp.core <- 0
    }
    else
    {
      seg.core <- which(X.condensed$domain$to == v.core)[1]
      tp.core <- 1
    }
  }
  
  # move observations to the specific point
  X.condensed$data$seg[is_cycle[X.condensed$data$seg]] <- 0
  X.condensed$data$seg <- X.condensed$data$seg - (X.condensed$data$seg > v.core)
  X.condensed$data$tp[X.condensed$data$seg == 0] <- tp.core
  X.condensed$data$seg[X.condensed$data$seg == 0] <- seg.core
  
  X.condensed$data$x <- X.condensed$domain$lines$ends$x0[X.condensed$data$seg] * (1 - X.condensed$data$tp) + X.condensed$domain$lines$ends$x1[X.condensed$data$seg] * X.condensed$data$tp
  X.condensed$data$y <- X.condensed$domain$lines$ends$y0[X.condensed$data$seg] * (1 - X.condensed$data$tp) + X.condensed$domain$lines$ends$y1[X.condensed$data$seg] * X.condensed$data$tp
  
  visited <- treeparty.visit(X.condensed)
  built <- treeparty.build(X.condensed, visited)
  
  
  iter <- 20L
  evals <- c("accuracy", "gini", "entropy")
  metrics_list_prop <- list(); metrics_list_comp <- list()
  correct.accuracy <- rep(0L, iter); correct.gini <- rep(0L, iter); correct.entropy <- rep(0L, iter)
  correct2.accuracy <- rep(0L, iter); correct2.gini <- rep(0L, iter); correct2.entropy <- rep(0L, iter)
  
  
  for (eval in evals)
  {
    # ---  Proposed Method (Cycle Condensation AdaBoost) ---
    all_preds_prop <- factor(rep(NA, n_obs), levels = labels)
    all_scores_prop <- rep(0, n_obs)
    correct_prop <- rep(0L, iter)
    
    for (i in 1:10) {
      where_test <- which(folds == i)
      where_train <- which(folds != i)
      
      n_train <- length(where_train)
      w_train <- rep(1/n_train, n_train)
      
      all_preds_matrix_prop <- matrix(0, nrow=n_obs, ncol=length(labels))
      colnames(all_preds_matrix_prop) <- labels
      
      for (it in 1:iter) {
        # Condensation
        count_matrix <- matrix(0, nrow=n_cycle, ncol=length(labels))
        colnames(count_matrix) <- labels
        for (idx in 1:n_train) {
          orig_i <- where_train[idx]
          k <- df_X$cycle_k[orig_i]
          count_matrix[k, y[orig_i]] <- count_matrix[k, y[orig_i]] + w_train[idx]
        }
        
        # Optimization
        best_p1 <- 1; best_p2 <- 2; best_impurity <- Inf
        for (opt_iter in 1:3) {
          for (cand2 in 1:n_cycle) {
            if (best_p1 == cand2) next
            p1 <- min(best_p1, cand2); p2 <- max(best_p1, cand2)
            groupA <- p1:(p2-1); groupB <- setdiff(1:n_cycle, groupA)
            cA <- colSums(count_matrix[groupA, , drop=FALSE])
            cB <- colSums(count_matrix[groupB, , drop=FALSE])
            imp <- (sum(cA)*calc_impurity(cA, eval) + sum(cB)*calc_impurity(cB, eval))
            if (imp < best_impurity) { best_impurity <- imp; best_p2 <- cand2 }
          }
          for (cand1 in 1:n_cycle) {
            if (cand1 == best_p2) next
            p1 <- min(cand1, best_p2); p2 <- max(cand1, best_p2)
            groupA <- p1:(p2-1); groupB <- setdiff(1:n_cycle, groupA)
            cA <- colSums(count_matrix[groupA, , drop=FALSE])
            cB <- colSums(count_matrix[groupB, , drop=FALSE])
            imp <- (sum(cA)*calc_impurity(cA, eval) + sum(cB)*calc_impurity(cB, eval))
            if (imp < best_impurity) { best_impurity <- imp; best_p1 <- cand1 }
          }
        }
        
        # condensation
        w_train2 <- rep(0, n_obs)
        w_train2[where_train] <- w_train
        
        # tree-based classification
        counted <- treeparty.count(built, w_train2)
        split <- treeparty.split(counted, minbucket = 1, eval = eval)
        
        # prediction
        p1 <- min(best_p1, best_p2); p2 <- max(best_p1, best_p2)
        groupA <- p1:(p2-1); groupB <- setdiff(1:n_cycle, groupA)
        predA <- labels[which.max(colSums(count_matrix[groupA, , drop=FALSE]))]
        predB <- labels[which.max(colSums(count_matrix[groupB, , drop=FALSE]))]
        
        pred_train <- factor(rep(NA, n_train), levels = labels)
        for (idx in 1:n_train) {
          k <- df_X$cycle_k[where_train[idx]]
          pred_train[idx] <- if (k %in% groupA) predA else predB
        }
        error <- sum(w_train[pred_train != y[where_train]])
        
        # Prediction
        pred_train2 <- treeparty.predict(built, split, index = where_train)
        error2 <- sum(w_train[pred_train2 != y[where_train]])
        
        err_clamped <- max(1e-10, min(1 - 1e-10, error, error2))
        
        # update adaboost weights
        alpha <- log((1 - err_clamped) / err_clamped) + log(length(labels) - 1)
        
        if (error < error2)
          w_train <- w_train * exp(alpha * (pred_train != y[where_train]))
        else
          w_train <- w_train * exp(alpha * (pred_train2 != y[where_train]))
        w_train <- w_train / sum(w_train)
        
        # Recorld all voting
        if (error < error2)
          for (orig_i in 1:n_obs) {
            k <- df_X$cycle_k[orig_i]
            pred_k <- if (k %in% groupA) predA else predB
            all_preds_matrix_prop[orig_i, pred_k] <- all_preds_matrix_prop[orig_i, pred_k] + alpha
          }
        else
        {
          pred_all2 <- treeparty.predict(built, split)
          all_preds_matrix_prop[df_X$orig_id, pred_all2] <- all_preds_matrix_prop[df_X$orig_id, pred_all2] + alpha
        }
        
        # calculate current accuracy
        pred_at_it <- apply(all_preds_matrix_prop[where_test, , drop=FALSE], 1, which.max)
        correct_prop[it] <- correct_prop[it] + sum(pred_at_it == as.integer(y[where_test]))
      }
      
      all_preds_prop[where_test] <- labels[apply(all_preds_matrix_prop[where_test, , drop=FALSE], 1, which.max)]
      all_scores_prop[where_test] <- all_preds_matrix_prop[where_test, pos_col_idx]
    }
    
    metrics_list_prop[[eval]] <- get_metrics(y, factor(all_preds_prop, levels = labels), all_scores_prop)
    if (eval == "accuracy") correct.accuracy <- correct_prop
    if (eval == "gini")     correct.gini <- correct_prop
    if (eval == "entropy")  correct.entropy <- correct_prop
    
    
    # --- \Comparison Method (Euclidean AdaBoost) ---
    correct_comp <- rep(0L, iter)
    all_preds_comp <- rep(NA, n_obs)
    all_scores_comp <- rep(0, n_obs)
    
    for (i in 1:10) {
      where <- which(folds == i)
      adapted2 <- ada(x[-where, ], y[-where], n_rounds = iter, verbose = FALSE, progress = FALSE, tree_depth = 1, split = eval)
      
      preds_c <- matrix(0, length(where), length(labels))
      for (it in 1L:iter) {
        pred <- stats::predict(adapted2$trees[[it]], x[where, ], type = "class")
        pred <- as.integer(pred)
        for (f in 1:nrow(preds_c)) preds_c[f, pred[f]] <- preds_c[f, pred[f]] + adapted2$alphas[it]
        pred_at_it <- apply(preds_c, 1, which.max)
        correct_comp[it] <- correct_comp[it] + sum(pred_at_it == as.integer(y[where]))
      }
      
      all_preds_comp[where] <- labels[apply(preds_c, 1, which.max)]
      all_scores_comp[where] <- preds_c[, pos_col_idx]
    }
    
    metrics_list_comp[[eval]] <- get_metrics(y, factor(all_preds_comp, levels = labels), all_scores_comp)
    if (eval == "accuracy") correct2.accuracy <- correct_comp
    if (eval == "gini")     correct2.gini <- correct_comp
    if (eval == "entropy")  correct2.entropy <- correct_comp
  }
  
  for (eval in evals)
  {
    # --- Proposed Method (Cycle Condensation AdaBoost) ---
    all_preds_prop <- factor(rep(NA, n_obs), levels = labels)
    all_scores_prop <- rep(0, n_obs)
    correct_prop <- rep(0L, iter)
    
    for (i in 1:10)
    {
      where_test <- which(folds == i)
      where_train <- which(folds != i)
      
      n_train <- length(where_train)
      w_train <- ifelse(folds == i, 0, 1/n_train)
      
      all_preds_matrix_prop <- matrix(0, nrow=n_obs, ncol=length(labels))
      colnames(all_preds_matrix_prop) <- labels
      
      for (it in 1:iter)
      {
        split <- treeparty1.split(weight = w_train, cycle.cut = 0, n_obs = n_obs, 
                                  n_cycle = n_cycle, df_X = df_X, y = y, 
                                  labels = labels, eval = eval, built = built)
        if (!is.null(split$p1))
        {
          wA <- w_train
          wB <- w_train
          groupA <- split$p1:(split$p2-1); groupB <- setdiff(1:n_cycle, groupA)
          wA[df_X$cycle_k %in% groupB] <- 0
          wB[df_X$cycle_k %in% groupA] <- 0
          split$subtree <- treeparty1.split(weight = wA, cycle.cut =groupB[1], n_cycle=n_cycle, labels=labels, n_obs=n_obs, df_X=df_X, y=y, built=built)
          split$supertree <- treeparty1.split(weight = wB, cycle.cut =groupA[1], n_cycle=n_cycle, labels=labels, n_obs=n_obs, df_X=df_X, y=y, built=built)
        }
        else
        {
          on.subtree1 <- built$invorder >= built$invorder[split$root]
          on.subtree2 <- built$invorder < built$invorder[split$root] + built$subtree.size[split$root]
          on.subtree <- on.subtree1 & on.subtree2
          on.suptree <- !on.subtree
          
          weight.subtree <- w_train
          weight.subtree[on.suptree] <- 0
          weight.suptree <- w_train
          weight.suptree[on.subtree] <- 0
          split$subtree <- treeparty1.split(weight.subtree, split$cycle.cut, n_cycle=n_cycle, labels=labels, n_obs=n_obs, df_X=df_X, y=y, built=built)
          split$supertree <- treeparty1.split(weight.suptree, split$cycle.cut, n_cycle=n_cycle, labels=labels, n_obs=n_obs, df_X=df_X, y=y, built=built)
        }
        
        pred <- factor(rep(NA, n_obs), levels = labels)
        if (!is.null(split$p1))
        {
          on.subtree1 <- df_X$cycle_k >= split$p1
          on.subtree2 <- df_X$cycle_k < split$p2
          groupA <- on.subtree1 & on.subtree2
          groupB <- !groupA
        }
        else
        {
          on.subtree1 <- built$invorder >= built$invorder[split$root]
          on.subtree2 <- built$invorder < built$invorder[split$root] + built$subtree.size[split$root]
          groupA <- on.subtree1 & on.subtree2
          groupB <- !groupA
        }
        if (!is.null(split$subtree$p1))
        {
          on.subtree1 <- df_X$cycle_k >= split$subtree$p1
          on.subtree2 <- df_X$cycle_k < split$subtree$p2
          groupAA <- on.subtree1 & on.subtree2
          groupAB <- !groupAA
        }
        else
        {
          on.subtree1 <- built$invorder >= built$invorder[split$subtree$root]
          on.subtree2 <- built$invorder < built$invorder[split$subtree$root] + built$subtree.size[split$subtree$root]
          groupAA <- on.subtree1 & on.subtree2
          groupAB <- !groupAA
          
          if (is.null(split$subtree$pred) || is.na(split$subtree$pred))
          {
            split$subtree$predA <- split$subtree$pred.sub
            split$subtree$predB <- split$subtree$pred.sup
          }
          else
          {
            split$subtree$predA <- split$subtree$pred
            split$subtree$predB <- split$subtree$pred
          }
        }
        if (!is.null(split$supertree$p1))
        {
          on.subtree1 <- df_X$cycle_k >= split$supertree$p1
          on.subtree2 <- df_X$cycle_k < split$supertree$p2
          groupBA <- on.subtree1 & on.subtree2
          groupBB <- !groupBA
        }
        else
        {
          on.subtree1 <- built$invorder >= built$invorder[split$supertree$root]
          on.subtree2 <- built$invorder < built$invorder[split$supertree$root] + built$subtree.size[split$supertree$root]
          groupBA <- on.subtree1 & on.subtree2
          groupBB <- !groupBA
          
          if (is.null(split$supertree$pred) || is.na(split$supertree$pred))
          {
            split$supertree$predA <- split$supertree$pred.sub
            split$supertree$predB <- split$supertree$pred.sup
          }
          else
          {
            split$supertree$predA <- split$supertree$pred
            split$supertree$predB <- split$supertree$pred
          }
        }
        pred[groupA & groupAA] <- split$subtree$predA
        pred[groupA & groupAB] <- split$subtree$predB
        pred[groupB & groupBA] <- split$supertree$predA
        pred[groupB & groupBB] <- split$supertree$predB
        
        wrong <- pred != y
        error <- sum(w_train[wrong])
        
        err_clamped <- max(1e-10, min(1 - 1e-10, error))
        
        # update Adaboost weights
        alpha <- log((1 - err_clamped) / err_clamped) + log(length(labels) - 1)
        w_train <- w_train * exp(alpha * (pred != y))
        w_train <- w_train / sum(w_train)
        
        # record all votes
        for (idx in 1:n_obs)
          all_preds_matrix_prop[idx, pred[idx]] <- all_preds_matrix_prop[idx, pred[idx]] + alpha
        
        # Calculate accuracy
        pred_at_it <- apply(all_preds_matrix_prop[where_test, , drop=FALSE], 1, which.max)
        correct_prop[it] <- correct_prop[it] + sum(pred_at_it == as.integer(y[where_test]))
      }
      
      all_preds_prop[where_test] <- labels[apply(all_preds_matrix_prop[where_test, , drop=FALSE], 1, which.max)]
      all_scores_prop[where_test] <- all_preds_matrix_prop[where_test, pos_col_idx]
    }
    
    metrics_list_prop[[eval]] <- get_metrics(y, factor(all_preds_prop, levels = labels), all_scores_prop)
    if (eval == "accuracy") correct.accuracy <- correct_prop
    if (eval == "gini")     correct.gini <- correct_prop
    if (eval == "entropy")  correct.entropy <- correct_prop
    
    
    # --- Comparison Method (Euclidean AdaBoost) ---
    correct_comp <- rep(0L, iter)
    all_preds_comp <- rep(NA, n_obs)
    all_scores_comp <- rep(0, n_obs)
    
    for (i in 1:10) {
      where <- which(folds == i)
      adapted2 <- ada(x[-where, ], y[-where], n_rounds = iter, verbose = FALSE, progress = FALSE, tree_depth = 2, split = eval)
      
      preds_c <- matrix(0, length(where), length(labels))
      for (it in 1L:iter) {
        pred <- stats::predict(adapted2$trees[[it]], x[where, ], type = "class")
        pred <- as.integer(pred)
        for (f in 1:nrow(preds_c)) preds_c[f, pred[f]] <- preds_c[f, pred[f]] + adapted2$alphas[it]
        pred_at_it <- apply(preds_c, 1, which.max)
        correct_comp[it] <- correct_comp[it] + sum(pred_at_it == as.integer(y[where]))
      }
      
      all_preds_comp[where] <- labels[apply(preds_c, 1, which.max)]
      all_scores_comp[where] <- preds_c[, pos_col_idx]
    }
    
    metrics_list_comp[[eval]] <- get_metrics(y, factor(all_preds_comp, levels = labels), all_scores_comp)
    if (eval == "accuracy") correct2.accuracy <- correct_comp
    if (eval == "gini")     correct2.gini <- correct_comp
    if (eval == "entropy")  correct2.entropy <- correct_comp
  }
  
  # -------------------------------------------------------------------------
  # KNN method
  # -------------------------------------------------------------------------
  correct3 <- rep(0L, iter)
  all_preds_knn <- rep(NA, n_obs)
  all_scores_knn <- rep(0, n_obs)
  
  for (i in 1:10) {
    where <- which(folds == i)
    for (it in seq(1L, iter, 2)) {
      knn_res <- knn(x[-where, ], x[where, ], y[-where], k = it, prob = TRUE)
      correct3[it] <- correct3[it] + sum(knn_res == y[where])
      
      all_preds_knn[where] <- as.character(knn_res)
      probs <- attr(knn_res, "prob")
      all_scores_knn[where] <- ifelse(knn_res == "main", probs, 1 - probs)
    }
  }
  metrics_knn <- get_metrics(y, factor(all_preds_knn, levels = labels), all_scores_knn)
  
  
  # check correct/incorrect
  correct_prop_idx <- (all_preds_prop == y)
  correct_comp_idx <- (all_preds_comp == y)
  
  # -------------------------------------------------------------------------
  # Generate res_final
  # -------------------------------------------------------------------------
  denom <- n_obs / 100
  
  res <- rbind(
    correct.accuracy[1:4 * 5], correct.gini[1:4 * 5], correct.entropy[1:4 * 5], 
    correct2.accuracy[1:4 * 5], correct2.gini[1:4 * 5], correct2.entropy[1:4 * 5],
    correct3[c(5, 9, 15, 19)]
  )
  
  res <- cbind(res, c(
    max(correct.accuracy), max(correct.gini), max(correct.entropy), 
    max(correct2.accuracy), max(correct2.gini), max(correct2.entropy),
    max(correct3)
  )) / denom
  
  extra_metrics <- matrix(NA, nrow = 7, ncol = 3)
  colnames(extra_metrics) <- c("F1", "B_Acc", "AUC")
  
  extra_metrics[1, ] <- unlist(metrics_list_prop[["accuracy"]])
  extra_metrics[2, ] <- unlist(metrics_list_prop[["gini"]])
  extra_metrics[3, ] <- unlist(metrics_list_prop[["entropy"]])
  extra_metrics[4, ] <- unlist(metrics_list_comp[["accuracy"]])
  extra_metrics[5, ] <- unlist(metrics_list_comp[["gini"]])
  extra_metrics[6, ] <- unlist(metrics_list_comp[["entropy"]])
  extra_metrics[7, ] <- unlist(metrics_knn)
  
  res_final <- cbind(res, extra_metrics)
  rownames(res_final) <- c("MR", "Gini", "Entropy", "MR_C", "Gini_C", "Entropy_C", "KNN") 
  colnames(res_final)[1:5] <- c("5", "10", "15", "20", "MAX_Acc")
  
  res_final_list[[seed]] <- res_final
}

cat("\n[ (Average over 100 seeds)]\n")
res_final_mean <- Reduce("+", res_final_list) / length(res_final_list)

saveRDS(res_final_list, "res_final_list.RDS")
round(res_final_mean, 3)
write.csv(res_final_mean, "res_final_mean.csv")
