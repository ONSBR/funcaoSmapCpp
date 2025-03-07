#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix ponderacao_temporal_cenario_cpp(NumericMatrix serie_temporal, NumericVector kt, int kt_max, int kt_min) {
  
  int N = kt_min + kt_max + 1;
  double soma_kt = sum(kt);
  int num_cenarios = serie_temporal.ncol(); // Linhas representam cenários
  int num_dias = serie_temporal.nrow(); // Colunas representam dias
  
  // A nova matriz terá menos colunas após aplicar a janela de tamanho N
  NumericMatrix serie_temporal_ponderada(num_cenarios, num_dias - N + 1);
  
  // Inverte a ordem dos pesos kt
  NumericVector pesos = rev(kt[Range(2 - kt_max, 2 + kt_min)]);
  
  // Loop sobre cada cenário (linha)
  for (int i = 0; i < num_cenarios; ++i) {
    NumericVector linha_atual = serie_temporal(_, i); // Série temporal para um cenário específico
    
    // Loop para calcular a ponderação temporal ao longo dos dias
    for (int j = 0; j < num_dias - N + 1; ++j) {
      NumericVector janela = linha_atual[Range(j, j + N - 1)];
      double ponderado = sum(janela * pesos) / sum(pesos);
      serie_temporal_ponderada(i, j) = ponderado * soma_kt;
    }
  }
  
  return serie_temporal_ponderada;
}
