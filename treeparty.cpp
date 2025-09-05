#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]

NumericMatrix treeparty_count(IntegerVector marks, int labels, IntegerVector order, IntegerVector parent, NumericVector weight)
{
  NumericMatrix res(order.length(), labels);
  for (int i = 0; i < marks.length(); i++)
    res(i, marks[i] - 1) = weight[i];
  for (int i = order.length() - 1; i; i--)
    res.row(parent[order[i] - 1] - 1) = res.row(parent[order[i] - 1] - 1) + res.row(order[i] - 1);
  return res;
}
/*
RObject treeparty_predict(List build, List split, IntegerVector index)
{
  IntegerVector marks = build[1];
  CharacterVector labels = marks.attr("levels");
  IntegerVector invorder = build[5];
  IntegerVector subtreeSize = build[6];
  
  IntegerVector res(index.length());
  
  List* branch;
  
  res.attr("class") = "factor";
  res.attr("levels") = labels;
  return res;
}
*/