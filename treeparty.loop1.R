n_train <- length(where_train)
w_train <- rep(1/n_train, n_train)
all_preds_matrix_prop <- matrix(0, nrow=n_obs, ncol=length(labels))
colnames(all_preds_matrix_prop) <- labels

# 임시 함수
treeparty1.split <- function(weight, cycle.cut = 0)
{
  # 1B. 가지 응축 (Condensation)
  count_matrix <- matrix(0, nrow=n_cycle, ncol=length(labels))
  colnames(count_matrix) <- labels
  for (idx in 1:n_obs)
  {
    k <- df_X$cycle_k[idx]
    count_matrix[k, y[idx]] <- count_matrix[k, y[idx]] + weight[idx]
  }
  
  # 2B. 2-Point Split 교대 최적화
  # cycle.cut은 이전에 사이클을 잘라서 현재 그래프에 없는 부분이다.
  if (cycle.cut)
  {
    best_p <- 1; best_impurity <- Inf
    for (cand in 1:n_cycle) {
      if (cycle.cut == cand) next
      p1 <- min(cycle.cut, cand); p2 <- max(cycle.cut, cand)
      groupA <- p1:(p2-1); groupB <- setdiff(1:n_cycle, groupA)
      cA <- colSums(count_matrix[groupA, , drop=FALSE])
      cB <- colSums(count_matrix[groupB, , drop=FALSE])
      imp <- (sum(cA)*calc_impurity(cA, eval) + sum(cB)*calc_impurity(cB, eval))
      if (imp < best_impurity) { best_impurity <- imp; best_p <- cand }
    }
  }
  else
  {
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
  }
  best_impurity <- best_impurity / sum(weight)
  
  # 2C. 트리 기반 분류
  counted <- treeparty.count(built, weight)
  split <- treeparty.split(counted, minbucket = 1, eval = eval)
  if (counted$obs < 2)
    split$impurity = Inf
  else if (eval == "accuracy")
    split$impurity <- 1 - split$accuracy
  
  # 3. 불순도 비교
  if (best_impurity < split$impurity)
  {
    # 사이클 절단
    res <- list()
    if (cycle.cut)
    {
      res$p1 <- min(cycle.cut, best_p)
      res$p2 <- max(cycle.cut, best_p)
    }
    else
    {
      res$p1 <- min(best_p1, best_p2)
      res$p2 <- max(best_p1, best_p2)
    }
    groupA <- p1:(p2-1); groupB <- setdiff(1:n_cycle, groupA)
    res$predA <- labels[which.max(colSums(count_matrix[groupA, , drop=FALSE]))]
    res$predB <- labels[which.max(colSums(count_matrix[groupB, , drop=FALSE]))]
    res$impurity <- best_impurity
    res$cycle.cut <- cycle.cut
  }
  else
  {
    # 트리 절단
    res <- split
    res$cycle.cut <- cycle.cut
  }
  return(res)
}



for (it in 1:iter) {
  
  
  # 3B. 그룹별 예측
  if (cycle.cut)
  {
    p1 <- min(cycle.cut, best_p)
    p2 <- max(cycle.cut, best_p)
  }
  else
  {
    p1 <- min(best_p1, best_p2)
    p2 <- max(best_p1, best_p2)
  }
  groupA <- p1:(p2-1); groupB <- setdiff(1:n_cycle, groupA)
  predA <- labels[which.max(colSums(count_matrix[groupA, , drop=FALSE]))]
  predB <- labels[which.max(colSums(count_matrix[groupB, , drop=FALSE]))]
  
  pred_train <- factor(rep(NA, n_train), levels = labels)
  for (idx in 1:n_train) {
    k <- df_X$cycle_k[where_train[idx]]
    pred_train[idx] <- if (k %in% groupA) predA else predB
  }
  error <- sum(w_train[pred_train != y[where_train]])
  
  # 3C. 그룹별 예측
  pred_train2 <- treeparty.predict(built, split, index = where_train)
  error2 <- sum(w_train[pred_train2 != y[where_train]])
  
  err_clamped <- max(1e-10, min(1 - 1e-10, error, error2))
  
  # 4. AdaBoost 가중치 업데이트
  
  alpha <- log((1 - err_clamped) / err_clamped) + log(length(labels) - 1)
  
  if (error < error2)
  {
    cat("Cycle cut ", exp(alpha), "\n")
    w_train <- w_train * exp(alpha * (pred_train != y[where_train]))
  }
  else
  {
    cat("Tree cut ", exp(alpha), "\n")
    w_train <- w_train * exp(alpha * (pred_train2 != y[where_train]))
  }
  w_train <- w_train / sum(w_train)
  
  # 테스트 셋 포함 전체 투표 기록
  if (error < error2)
    for (orig_i in 1:n_obs) {
      k <- df_X$cycle_k[orig_i]
      pred_k <- if (k %in% groupA) predA else predB
      all_preds_matrix_prop[orig_i, pred_k] <- all_preds_matrix_prop[orig_i, pred_k] + alpha
    }
  else
  {
    pred_all2 <- treeparty.predict(built, split)
    for (idx in 1:n_obs) {
      all_preds_matrix_prop[idx, pred_all2[idx]] <- all_preds_matrix_prop[idx, pred_all2[idx]] + alpha
    }
  }  
  # Fold별 현재 iteration 테스트 적중 수 집계
  pred_at_it <- apply(all_preds_matrix_prop[where_test, , drop=FALSE], 1, which.max)
  correct_prop[it] <- correct_prop[it] + sum(pred_at_it == as.integer(y[where_test]))
}