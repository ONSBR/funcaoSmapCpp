# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

ponderacao_temporal_cenario_cpp <- function(serie_temporal, kt, kt_max, kt_min) {
    .Call(`_funcaoSmapCpp_ponderacao_temporal_cenario_cpp`, serie_temporal, kt, kt_max, kt_min)
}

ponderacao_temporal_cpp <- function(serie_temporal, kt, kt_max, kt_min) {
    .Call(`_funcaoSmapCpp_ponderacao_temporal_cpp`, serie_temporal, kt, kt_max, kt_min)
}

propaga_tv_cpp <- function(vazao_montante, vazao_jusante, tempo_viagem) {
    .Call(`_funcaoSmapCpp_propaga_tv_cpp`, vazao_montante, vazao_jusante, tempo_viagem)
}

propaga_muskingum_cpp <- function(vazao_montante, vazao_jusante, n, coeficientes) {
    .Call(`_funcaoSmapCpp_propaga_muskingum_cpp`, vazao_montante, vazao_jusante, n, coeficientes)
}

rodada_varios_dias_cpp2 <- function(modelo, inicializacao, area, precipitacao, evapotranspiracao, Emarg, numero_dias) {
    .Call(`_funcaoSmapCpp_rodada_varios_dias_cpp2`, modelo, inicializacao, area, precipitacao, evapotranspiracao, Emarg, numero_dias)
}

rodada_cenarios_dias_cpp2 <- function(modelo, inicializacao, area, precipitacao, evapotranspiracao, Emarg, numero_dias, numero_cenarios) {
    .Call(`_funcaoSmapCpp_rodada_cenarios_dias_cpp2`, modelo, inicializacao, area, precipitacao, evapotranspiracao, Emarg, numero_dias, numero_cenarios)
}

rodada_pmur_cpp <- function(modelo, inicializacao, area, precipitacao, evapotranspiracao, Emarg, numero_dias) {
    .Call(`_funcaoSmapCpp_rodada_pmur_cpp`, modelo, inicializacao, area, precipitacao, evapotranspiracao, Emarg, numero_dias)
}

rodada_pmur_cpp_cenario <- function(modelo, inicializacao, area, precipitacao, evapotranspiracao, Emarg, numero_dias, numero_cenarios) {
    .Call(`_funcaoSmapCpp_rodada_pmur_cpp_cenario`, modelo, inicializacao, area, precipitacao, evapotranspiracao, Emarg, numero_dias, numero_cenarios)
}

