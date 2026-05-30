library(spatstat)
library(spatstat.linnet)
source('C:/Users/donts/Desktop/adaboost.linnet-main/treeparty.R', chdir = TRUE)
source('C:/Users/donts/Desktop/adaboost.linnet-main/comparison.split.R', chdir = TRUE)
library(class)

library(pROC) # AUC 계산용


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
  # actual, predicted는 factor ("Good", "Bad")
  cm <- table(Actual = actual, Predicted = predicted)
  
  # 클래스가 2개인 경우 (Good, Bad)
  # Good을 Positive(1)로 가정할 경우의 인덱스 설정
  tp <- cm["Good", "Good"]
  tn <- cm["Bad", "Bad"]
  fp <- cm["Bad", "Good"]
  fn <- cm["Good", "Bad"]
  
  precision <- tp / (tp + fp)
  recall <- tp / (tp + fn)
  f1 <- 2 * (precision * recall) / (precision + recall)
  
  # Balanced Accuracy: (Sensitivity + Specificity) / 2
  sensitivity <- recall
  specificity <- tn / (tn + fp)
  balanced_acc <- (sensitivity + specificity) / 2
  
  # AUC 계산
  auc_val <- NA
  if (!is.null(scores)) {
    # scores는 Good일 확률 혹은 가중치 합
    auc_val <- as.numeric(auc(actual, scores, levels = c("Bad", "Good"), direction = "<"))
  }
  
  return(list(F1 = f1, B_Acc = balanced_acc, AUC = auc_val))
}


# 사이클 정점들의 리스트(cycle_nodes)와 C++ 결과(count_matrix)가 주어졌을 때
condense_cycle_weights <- function(count_matrix, cycle_nodes, parent_array) {
  
  n_cycle <- length(cycle_nodes)
  # labels(분류 클래스 수) 개수만큼 열을 가진 빈 행렬 생성
  condensed_weights <- matrix(0, nrow = n_cycle, ncol = ncol(count_matrix))
  
  for (i in 1:n_cycle) {
    v <- cycle_nodes[i]
    total_w <- count_matrix[v, ]
    
    # v의 자식들 중, '사이클 위에 있는 자식' 찾기
    # parent_array를 뒤져서 부모가 v인 노드들 중 cycle_nodes에 포함된 노드
    cycle_children <- which(parent_array == v)
    cycle_children <- cycle_children[cycle_children %in% cycle_nodes]
    
    # 사이클 위로 이어지는 가중치만 뺀다. (바깥쪽 가지들의 가중치만 남김)
    for (child in cycle_children) {
      total_w <- total_w - count_matrix[child, ]
    }
    
    condensed_weights[i, ] <- total_w
  }
  
  return(condensed_weights) # 사이클 위 정점들로 완벽히 응축된 데이터!
}

# 1. 정점(vertices) 정의 (사각형 모양의 루프)
v <- ppp(x = c(-3, 3, 3, -3, 1, -1), 
         y = c(-3, -3, 3, 3, -4, 4), 
         window = owin(c(-4, 4), c(-4, 4)))

# 2. 에지(edges) 정의 (1-2-3-4-1 로 이어지는 루프 + 위아래로 뻗은 가지)
# 1-2, 2-3, 3-4, 4-1 (루프 형성)
# 1-5, 3-6 (외부 연결 가지)
edge <- matrix(c(1, 2,
                 2, 3,
                 3, 4,
                 4, 1,
                 1, 5,
                 3, 6), nrow = 6, ncol = 2, byrow = TRUE)

loop_net <- linnet(v, edges = edge)

# 3. 경로 설정
# 메인 경로는 아래쪽 가지에서 루프 절반을 타고 위로 올라가는 경로로 설정
mainroute <- c(6, 1, 2) # 에지 인덱스 (1번 에지: 1-2, 5번 에지: 1-5 등)
# 실제 에지 번호 확인은 plot(loop_net, show.edge.numbers=T)로 가능합니다.

# 4. 데이터 생성 (기존 maze와 동일한 방식)
set.seed(0)
X <- runiflpp(n = 300, L = loop_net)
# 에지 번호에 따라 main/sub 마킹
X$data$marks <- ifelse(X$data$seg %in% mainroute, "main", "sub")

# 5. Junction(교차점) 근처에서 변이 발생 (정점 1, 3이 교차점)
junction_nodes <- c(1, 3)
from_j <- X$domain$from[X$data$seg] %in% junction_nodes
to_j <- X$domain$to[X$data$seg] %in% junction_nodes
mutation_prob <- (1 - X$data$tp) * from_j + X$data$tp * to_j
X$data$marks[runif(300) < mutation_prob] <- "junction"

X$data$marks <- factor(X$data$marks, levels = c("main", "sub", "junction"))

# 시각화 확인
plot(X, cols = c("blue", "red", "green"), main = "Loop Network Simulation",
     cex = 2)

table(X$data$marks)

# 최대 정확도 83.33%

# ---------------------------------------------------------
# 1. 동적 에지 판별 및 사이클/가지 데이터 분리 (Robust 버전)
# ---------------------------------------------------------

labels <- levels(X$data$marks)

# [수정] spatstat이 에지 번호를 섞어버려도 찾을 수 있도록 동적 탐색
# 정점 1, 2, 3, 4 사이를 연결하는 에지가 사이클입니다.
is_cycle <- (X$domain$from %in% 1:4) & (X$domain$to %in% 1:4)
cycle_seg_ids <- which(is_cycle)
branch_seg_ids <- which(!is_cycle)

# 동적으로 찾은 ID를 바탕으로 데이터 분리
cycle_data <- X$data[X$data$seg %in% cycle_seg_ids, ] 
branch_data <- X$data[X$data$seg %in% branch_seg_ids, ]      

# 사이클 물리적 정렬
cycle_data <- cycle_data[order(cycle_data$seg, cycle_data$tp), ]
n_cycle <- nrow(cycle_data)

# 가중치 행렬 초기화
count_matrix <- matrix(0, nrow = n_cycle, ncol = length(labels))
colnames(count_matrix) <- labels

for(i in seq_len(n_cycle)) {
  lbl <- as.character(cycle_data$marks[i])
  if(!is.na(lbl)) count_matrix[i, lbl] <- 1
}

# ---------------------------------------------------------
# 2. 가지(Branch) 가중치 응축 (Condensation)
# ---------------------------------------------------------

# 정점 1에 연결된 외부 가지 찾기
branch1_id <- which(!is_cycle & (X$domain$from == 1 | X$domain$to == 1))
if(length(branch1_id) > 0) {
  for(lbl in labels) {
    b_count <- sum(branch_data$seg == branch1_id & branch_data$marks == lbl, na.rm = TRUE)
    count_matrix[1, lbl] <- count_matrix[1, lbl] + b_count
  }
}

# 정점 3에 연결된 외부 가지 찾기 및 가장 가까운 사이클 데이터 노드 찾기
branch3_id <- which(!is_cycle & (X$domain$from == 3 | X$domain$to == 3))
idx_node3 <- which(cycle_data$seg %in% which((X$domain$from == 3 | X$domain$to == 3) & is_cycle))[1]
if(is.na(idx_node3)) idx_node3 <- max(1, n_cycle) # 안전장치

if(length(branch3_id) > 0) {
  for(lbl in labels) {
    b_count <- sum(branch_data$seg == branch3_id & branch_data$marks == lbl, na.rm = TRUE)
    count_matrix[idx_node3, lbl] <- count_matrix[idx_node3, lbl] + b_count
  }
}

cat("사이클 응축 완료. 사이클 상의 노드 수:", n_cycle, "\n")

# ---------------------------------------------------------
# 3. 2-Point Split 교대 최적화 (Alternating Optimization)
# ---------------------------------------------------------

calc_gini <- function(counts) {
  total <- sum(counts)
  if(total == 0) return(0)
  return(1 - sum((counts / total)^2))
}

evaluate_2point_split <- function(p1, p2, count_matrix) {
  n <- nrow(count_matrix)
  if (n < 2 || p1 == p2) return(Inf)
  
  # [수정] 인덱스 오류(out of bounds) 원천 차단
  # p1을 무조건 작은 쪽으로 맞춤 (어차피 원형이므로 무관)
  if(p1 > p2) { temp <- p1; p1 <- p2; p2 <- temp }
  
  # 사이클 분할 그룹핑 (setdiff 활용으로 인덱스 오류 절대 발생 안함)
  groupA_idx <- p1:(p2 - 1)
  groupB_idx <- setdiff(1:n, groupA_idx)
  
  countsA <- colSums(count_matrix[groupA_idx, , drop = FALSE])
  countsB <- colSums(count_matrix[groupB_idx, , drop = FALSE])
  
  total <- sum(countsA) + sum(countsB)
  if(total == 0) return(Inf)
  
  gini <- (sum(countsA) / total) * calc_gini(countsA) + 
    (sum(countsB) / total) * calc_gini(countsB)
  return(gini)
}

best_p1 <- 1
best_p2 <- 2
best_impurity <- Inf

# 교대 최적화 진행
if (n_cycle >= 2) {
  for (iter in 1:5) { 
    for (cand2 in 1:n_cycle) {
      imp <- evaluate_2point_split(best_p1, cand2, count_matrix)
      if (imp < best_impurity) { best_impurity <- imp; best_p2 <- cand2 }
    }
    for (cand1 in 1:n_cycle) {
      imp <- evaluate_2point_split(cand1, best_p2, count_matrix)
      if (imp < best_impurity) { best_impurity <- imp; best_p1 <- cand1 }
    }
  }
}

cat(sprintf("최적 분할점 도출: Node %d, Node %d (Gini Impurity: %.4f)\n", best_p1, best_p2, best_impurity))

# ---------------------------------------------------------
# 분할 결과 시각화 (수정된 안전한 좌표 추출 방식)
# ---------------------------------------------------------

# spatstat의 lpp 객체에서 x, y 좌표를 포함한 전체 데이터프레임 추출
df_X <- as.data.frame(X)

# 사이클 데이터만 추출하여 동일한 조건(seg, tp)으로 정렬
cycle_df <- df_X[df_X$seg %in% cycle_seg_ids, ]
cycle_df <- cycle_df[order(cycle_df$seg, cycle_df$tp), ]

plot(X, cols = c("lightgray", "lightgray", "lightgray"), main = "2-Point Split Result on Cycle")

if (n_cycle >= 2) {
  # 정렬된 cycle_df에서 best_p1, best_p2 번째 행의 x, y 좌표를 직접 가져옵니다.
  points(cycle_df$x[best_p1], cycle_df$y[best_p1], col="blue", pch=19, cex=3)
  points(cycle_df$x[best_p2], cycle_df$y[best_p2], col="red", pch=19, cex=3)
  
  legend("topright", legend=c(paste("Split 1 (Node", best_p1, ")"), 
                              paste("Split 2 (Node", best_p2, ")")), 
         col=c("blue", "red"), pch=19)
}

# -------------------------------------------------------------------------
# 1. 평가 지표 및 불순도 계산 함수
# -------------------------------------------------------------------------
library(pROC)
library(class)

get_metrics <- function(actual, predicted, scores = NULL) {
  cm <- table(Actual = actual, Predicted = predicted)
  classes <- rownames(cm)
  
  pos_class <- "main"
  neg_class <- "sub"
  
  tp <- if(pos_class %in% classes && pos_class %in% colnames(cm)) cm[pos_class, pos_class] else 0
  tn <- if(neg_class %in% classes && neg_class %in% colnames(cm)) cm[neg_class, neg_class] else 0
  fp <- if(neg_class %in% classes && pos_class %in% colnames(cm)) cm[neg_class, pos_class] else 0
  fn <- if(pos_class %in% classes && neg_class %in% colnames(cm)) cm[pos_class, neg_class] else 0
  
  precision <- ifelse((tp + fp) == 0, 0, tp / (tp + fp))
  recall <- ifelse((tp + fn) == 0, 0, tp / (tp + fn))
  f1 <- ifelse((precision + recall) == 0, 0, 2 * (precision * recall) / (precision + recall))
  
  sensitivity <- recall
  specificity <- ifelse((tn + fp) == 0, 0, tn / (tn + fp))
  balanced_acc <- (sensitivity + specificity) / 2
  
  auc_val <- NA
  if (!is.null(scores)) {
    try({
      auc_val <- as.numeric(auc(actual, scores, levels = c(neg_class, pos_class), direction = "<"))
    }, silent = TRUE)
  }
  return(list(F1 = f1, B_Acc = balanced_acc, AUC = auc_val))
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

# -------------------------------------------------------------------------
# 2. 루프 네트워크 데이터 준비 및 사이클 투영(Projection)
# ---------------------------------------------------------
set.seed(0)
# loop_net은 앞서 생성된 객체를 사용합니다.
X <- runiflpp(n = 300, L = loop_net)
X$data$marks <- ifelse(X$data$seg %in% mainroute, "main", "sub")

junction_nodes <- c(1, 3)
from_j <- X$domain$from[X$data$seg] %in% junction_nodes
to_j <- X$domain$to[X$data$seg] %in% junction_nodes
mutation_prob <- (1 - X$data$tp) * from_j + X$data$tp * to_j
X$data$marks[runif(300) < mutation_prob] <- "junction"
X$data$marks <- factor(X$data$marks, levels = c("main", "sub", "junction"))

# ---- 학생 논리: 모든 관측치를 사이클 위로 응축하기 위한 인덱싱 ----
df_X <- as.data.frame(X)
df_X$orig_id <- 1:nrow(df_X)

is_cycle <- (X$domain$from %in% 1:4) & (X$domain$to %in% 1:4)
cycle_seg_ids <- which(is_cycle)
branch1_id <- which(!is_cycle & (X$domain$from == 1 | X$domain$to == 1))
branch3_id <- which(!is_cycle & (X$domain$from == 3 | X$domain$to == 3))

cycle_data <- df_X[df_X$seg %in% cycle_seg_ids, ]
cycle_data <- cycle_data[order(cycle_data$seg, cycle_data$tp), ]
n_cycle <- nrow(cycle_data)
cycle_data$cycle_k <- 1:n_cycle

df_X$cycle_k <- NA
df_X$cycle_k[match(cycle_data$orig_id, df_X$orig_id)] <- cycle_data$cycle_k
df_X$cycle_k[df_X$seg %in% branch1_id] <- 1
node3_k <- which(cycle_data$seg %in% which((X$domain$from == 3 | X$domain$to == 3) & is_cycle))[1]
if(is.na(node3_k)) node3_k <- max(1, n_cycle)
df_X$cycle_k[df_X$seg %in% branch3_id] <- node3_k
# ---------------------------------------------------------------------

n_obs <- nrow(df_X)
tb <- table(df_X$marks)
folds <- rep(1, n_obs)
for (i in 1:length(tb)) {
  if(tb[i] > 0) folds[df_X$marks == names(tb[i])] <- as.integer((order(runif(tb[i])) - 1) / tb[i] * 10) + 1
}

iter <- 20L
evals <- c("accuracy", "gini", "entropy")
metrics_list_prop <- list(); metrics_list_comp <- list()
correct.accuracy <- rep(0L, iter); correct.gini <- rep(0L, iter); correct.entropy <- rep(0L, iter)
correct2.accuracy <- rep(0L, iter); correct2.gini <- rep(0L, iter); correct2.entropy <- rep(0L, iter)

y <- df_X$marks
x <- data.frame(x = df_X$x, y = df_X$y)
labels <- factor(levels(y), levels = levels(y))
pos_col_idx <- which(labels == "main")

# -------------------------------------------------------------------------
# 3. 10-fold CV 메인 루프 (Proposed vs Euclidean Comparison)
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
# This can be done by deleting an edge
# and moving the endpoints of some edges to a specific point

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

plot(X.condensed)

visited <- treeparty.visit(X.condensed)
built <- treeparty.build(X.condensed, visited)


for (eval in evals)
{
  cat("\n=== Evaluating with split criterion:", eval, "===\n")
  
  # --- [A] Proposed Method (Cycle Condensation AdaBoost) ---
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
      # 1B. 가지 응축 (Condensation)
      
      count_matrix <- matrix(0, nrow=n_cycle, ncol=length(labels))
      colnames(count_matrix) <- labels
      for (idx in 1:n_train) {
        orig_i <- where_train[idx]
        k <- df_X$cycle_k[orig_i]
        count_matrix[k, y[orig_i]] <- count_matrix[k, y[orig_i]] + w_train[idx]
      }
      
      # 2B. 2-Point Split 교대 최적화
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
      
      # 1C. 사이클 응축
      w_train2 <- rep(0, n_obs)
      w_train2[where_train] <- w_train
      
      # 2C. 트리 기반 분류
      counted <- treeparty.count(built, w_train2)
      split <- treeparty.split(counted, minbucket = 1, eval = eval)
      
      # 3B. 그룹별 예측
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
      
      # 3C. 그룹별 예측
      pred_train2 <- treeparty.predict(built, split, index = where_train)
      error2 <- sum(w_train[pred_train2 != y[where_train]])
      
      err_clamped <- max(1e-10, min(1 - 1e-10, error, error2))
      
      # 4. AdaBoost 가중치 업데이트
      
      alpha <- log((1 - err_clamped) / err_clamped) + log(length(labels) - 1)
      
      if (error < error2)
        w_train <- w_train * exp(alpha * (pred_train != y[where_train]))
      else
        w_train <- w_train * exp(alpha * (pred_train2 != y[where_train]))
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
    
    all_preds_prop[where_test] <- labels[apply(all_preds_matrix_prop[where_test, , drop=FALSE], 1, which.max)]
    all_scores_prop[where_test] <- all_preds_matrix_prop[where_test, pos_col_idx]
    cat("  Proposed Fold", i, "/ 10 done.\n")
  }
  
  #metrics_list_prop[[eval]] <- get_metrics(y, factor(all_preds_prop, levels = labels), all_scores_prop)
  if (eval == "accuracy") correct.accuracy <- correct_prop
  if (eval == "gini")     correct.gini <- correct_prop
  if (eval == "entropy")  correct.entropy <- correct_prop
  
  
  # --- [B] Comparison Method (Euclidean AdaBoost) ---
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
    cat("  Comparison Fold", i, "/ 10 done.\n")
  }
  
  #metrics_list_comp[[eval]] <- get_metrics(y, factor(all_preds_comp, levels = labels), all_scores_comp)
  if (eval == "accuracy") correct2.accuracy <- correct_comp
  if (eval == "gini")     correct2.gini <- correct_comp
  if (eval == "entropy")  correct2.entropy <- correct_comp
}

for (eval in evals)
{
  cat("\n=== Evaluating with split criterion:", eval, "===\n")
  
  # --- [A] Proposed Method (Cycle Condensation AdaBoost) ---
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
      split <- treeparty1.split(w_train)
      if (is.null(split$cycle.cut))
      {
        wA <- w_train
        wB <- w_train
        groupA <- split$p1:(split$p2-1); groupB <- setdiff(1:n_cycle, groupA)
        wA[df_X$cycle_k %in% groupB] <- 0
        wB[df_X$cycle_k %in% groupA] <- 0
        split$subtree <- treeparty1.split(wA, groupB[1])
        split$supertree <- treeparty1.split(wB, groupA[1])
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
        split$subtree <- treeparty1.split(weight.subtree, split$cycle.cut)
        split$supertree <- treeparty1.split(weight.suptree, split$cycle.cut)
      }
      
      pred <- factor(rep(NA, n_obs), levels = labels)
      # 3. 그룹별 예측
      if (is.null(split$cycle.cut))
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
      if (is.null(split$subtree$cycle.cut))
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
        if (is.na(split$subtree$pred))
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
      if (is.null(split$supertree$cycle.cut))
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
        if (is.na(split$supertree$pred))
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
      
      # 4. AdaBoost 가중치 업데이트
      
      alpha <- log((1 - err_clamped) / err_clamped) + log(length(labels) - 1)
      w_train <- w_train * exp(alpha * (pred != y))
      w_train <- w_train / sum(w_train)
      
      # 테스트 셋 포함 전체 투표 기록
      for (idx in 1:n_obs)
        all_preds_matrix_prop[idx, pred[idx]] <- all_preds_matrix_prop[idx, pred[idx]] + alpha

      # Fold별 현재 iteration 테스트 적중 수 집계
      pred_at_it <- apply(all_preds_matrix_prop[where_test, , drop=FALSE], 1, which.max)
      correct_prop[it] <- correct_prop[it] + sum(pred_at_it == as.integer(y[where_test]))
    }
    
    all_preds_prop[where_test] <- labels[apply(all_preds_matrix_prop[where_test, , drop=FALSE], 1, which.max)]
    all_scores_prop[where_test] <- all_preds_matrix_prop[where_test, pos_col_idx]
    cat("  Proposed Fold", i, "/ 10 done.\n")
  }
  
  #metrics_list_prop[[eval]] <- get_metrics(y, factor(all_preds_prop, levels = labels), all_scores_prop)
  if (eval == "accuracy") correct.accuracy <- correct_prop
  if (eval == "gini")     correct.gini <- correct_prop
  if (eval == "entropy")  correct.entropy <- correct_prop
  
  
  # --- [B] Comparison Method (Euclidean AdaBoost) ---
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
    cat("  Comparison Fold", i, "/ 10 done.\n")
  }
  
  #metrics_list_comp[[eval]] <- get_metrics(y, factor(all_preds_comp, levels = labels), all_scores_comp)
  if (eval == "accuracy") correct2.accuracy <- correct_comp
  if (eval == "gini")     correct2.gini <- correct_comp
  if (eval == "entropy")  correct2.entropy <- correct_comp
}

# -------------------------------------------------------------------------
# 4. KNN 비교
# -------------------------------------------------------------------------
cat("\n=== Evaluating KNN ===\n")
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
#metrics_knn <- get_metrics(y, factor(all_preds_knn, levels = labels), all_scores_knn)

# -------------------------------------------------------------------------
# 5. 결과 테이블 (res_final) 결합 및 출력
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

round(res, 3)

extra_metrics <- matrix(NA, nrow = 7, ncol = 3)
colnames(extra_metrics) <- c("F1", "AUC", "B_Acc")

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

cat("\n[최종 평가 지표 결과]\n")
round(res_final, 3)