#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

// Função propaga_tv
// [[Rcpp::export]]
NumericVector propaga_tv_cpp(NumericVector vazao_montante, NumericVector vazao_jusante, double tempo_viagem) {
  int Nper = vazao_jusante.size();
  NumericVector vazao_propagada(Nper, 0.0);
  int lag_dias = std::ceil(tempo_viagem / 24.0);
  double fator_segundo_dia = ((24 * lag_dias) - tempo_viagem) / 24.0;
  double fator_primeiro_dia = (tempo_viagem - (24 * (lag_dias - 1))) / 24.0;
  
  // Loop: em R o loop inicia em (1+lag_dias) até Nper; 
  // aqui usamos índices 0-based: i = lag_dias até Nper-1.
  for (int i = lag_dias; i < Nper; i++) {
    if (lag_dias > 0) {
      // Conversão dos índices:
      //   iper (R) = i + 1 (C++), assim:
      //   vazao_montante[iper - lag_dias] -> vazao_montante[i - lag_dias]
      //   vazao_montante[iper - lag_dias + 1] -> vazao_montante[i - lag_dias + 1]
      vazao_propagada[i] = (vazao_montante[i - lag_dias + 1] * fator_segundo_dia +
                            vazao_montante[i - lag_dias] * fator_primeiro_dia) +
                            vazao_jusante[i];
    } else {
      // Caso lag_dias == 0 (tempo_viagem < 24)
      vazao_propagada[i] = (vazao_montante[i] * fator_segundo_dia +
                            vazao_montante[i] * fator_primeiro_dia) +
                            vazao_jusante[i];
    }
  }
  return vazao_propagada;
}


// Função propaga_muskingum
// [[Rcpp::export]]
NumericVector propaga_muskingum_cpp(NumericVector vazao_montante, NumericVector vazao_jusante, int n, NumericVector coeficientes) {
  int Nper = vazao_jusante.size();
  NumericMatrix vazao_passo_n(Nper, n + 1);
  
  // Inicializa a primeira linha: todos os elementos recebem o primeiro valor de vazao_montante
  for (int j = 0; j < (n + 1); j++) {
    vazao_passo_n(0, j) = vazao_montante[0];
  }
  
  // Define a primeira coluna com a vazao_montante (coluna 0 em C++)
  for (int i = 0; i < Nper; i++) {
    vazao_passo_n(i, 0) = vazao_montante[i];
  }
  
  NumericVector vazao_propagada(Nper, 0.0);
  
  // Loop: em R iper varia de 2 a Nper; aqui i varia de 1 a Nper-1
  for (int i = 1; i < Nper; i++) {
    // Para cada linha, as colunas de 2 até (n+1) em R correspondem a j de 1 até n (0-indexado)
    for (int j = 1; j < (n + 1); j++) {
      // Note que em R os coeficientes são indexados a partir de 1; 
      // em C++ usamos coeficientes[0], coeficientes[1] e coeficientes[2]
      vazao_passo_n(i, j) = vazao_passo_n(i, j - 1) * coeficientes[0] +
                            vazao_passo_n(i - 1, j - 1) * coeficientes[1] +
                            vazao_passo_n(i - 1, j) * coeficientes[2];
    }
    // A última coluna (n+1 em R, ou n em C++) compõe a vazão propagada somada à vazao_jusante
    vazao_propagada[i] = vazao_passo_n(i, n) + vazao_jusante[i];
  }
  
  return vazao_propagada;
}