#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector ponderacao_temporal_cpp(NumericVector serie_temporal, NumericVector kt, int kt_max, int kt_min) {
  
  int N = kt_min + kt_max + 1;
  double soma_kt = sum(kt);
  NumericVector serie_temporal_ponderada(serie_temporal.size() - N + 1);
  
  for (int i = 0; i < serie_temporal_ponderada.size(); ++i) {
    NumericVector janela = serie_temporal[Range(i, i + N - 1)];
    
    // Inverte a ordem dos pesos kt
    NumericVector pesos = rev(kt[Range(2 - kt_max, 2 + kt_min)]);
    
    double ponderado = sum(janela * pesos) / sum(pesos);
    serie_temporal_ponderada[i] = ponderado * soma_kt;
  }
  
  return serie_temporal_ponderada;
}