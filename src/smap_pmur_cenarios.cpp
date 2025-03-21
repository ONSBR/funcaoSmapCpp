#include <Rcpp.h>
using namespace Rcpp;

// C++ implementation of rodada_varios_dias function
// [[Rcpp::export]]
List rodada_pmur_cpp_cenario(NumericVector modelo, NumericMatrix inicializacao, double area,
                                 NumericMatrix precipitacao, NumericMatrix evapotranspiracao,
                                 NumericMatrix Emarg, int numero_dias, int numero_cenarios) {

	int ncols = 19; // Number of columns in the output matrix

	// Create the output matrix
	List listaSaida(numero_cenarios);
	NumericMatrix matrizSaida(numero_dias, ncols);
	colnames(matrizSaida) = CharacterVector::create("Qcalc", "Rsolo", "Rsup", "Rsup2", "Rsub",
												  "Es", "Er", "Rec", "Marg", "Ed", "Ed3", "Ed2",
												  "Eb", "Tu", "Qsup1", "Qplan", "Qsup2", "Qbase",
												  "Ai");

	double Capc_tmp = (modelo(3) / 100) * modelo(0);
	double Rsup_tmp = 0;
	double Pmur_tmp = (modelo(14) / 100) * Capc_tmp;
	double capacidade_armazenamento = 0;
	// Coefficients Ks
	double K_kts = pow(0.5, 1.0 / modelo(4));
	double K_1ts = pow(0.5, 1.0 / modelo(9));
	double K_2ts = pow(0.5, 1.0 / modelo(1));
	double K_2t2s = pow(0.5, 1.0 / modelo(6));
	double K_3ts = pow(0.5, 1.0 / modelo(10));
	
	for (int icenario = 0; icenario < numero_cenarios; icenario++){
		NumericMatrix matrizSaida(numero_dias, ncols);
		colnames(matrizSaida) = CharacterVector::create("Qcalc", "Rsolo", "Rsup", "Rsup2", "Rsub",
												  "Es", "Er", "Rec", "Marg", "Ed", "Ed3", "Ed2",
												  "Eb", "Tu", "Qsup1", "Qplan", "Qsup2", "Qbase",
												  "Ai");

		for (int idia = 0; idia < numero_dias; idia++) {

			capacidade_armazenamento = (modelo(0) - inicializacao(icenario, 4));
			matrizSaida(idia, 18) = capacidade_armazenamento * modelo(7) / 100;

			// Calculation of transfer functions
			matrizSaida(idia, 13) = inicializacao(icenario, 4) / modelo(0); // Eq.19 Manual

			if (precipitacao(icenario, idia) > matrizSaida(idia, 18)) {
			  matrizSaida(idia, 5) = pow(precipitacao(icenario, idia) - matrizSaida(idia, 18), 2) /
									 (precipitacao(icenario, idia) - matrizSaida(idia, 18) + capacidade_armazenamento); // Eq.11
			} else {
			  matrizSaida(idia, 5) = 0;
			}

			if (inicializacao(icenario, 4) > Pmur_tmp){
				if ((precipitacao(icenario, idia) - matrizSaida(idia, 5)) > evapotranspiracao(icenario, idia)) {
				matrizSaida(idia, 6) = evapotranspiracao(icenario, idia);
				} else {
				matrizSaida(idia, 6) = (precipitacao(icenario, idia) - matrizSaida(idia, 5)) +
										(evapotranspiracao(icenario, idia) - (precipitacao(icenario, idia) - matrizSaida(idia, 5))) *
											matrizSaida(idia, 13); // Eq.12
				}
			} else{
				matrizSaida(idia, 6) = 0;
			}

			if (inicializacao(icenario, 4) > Capc_tmp) {
			  matrizSaida(idia, 7) = (modelo(2) / 100) * matrizSaida(idia, 13) * (inicializacao(icenario, 4) - Capc_tmp); // Eq.13
			} else {
			  matrizSaida(idia, 7) = 0;
			}

			// 4th Reservoir
			if (inicializacao(icenario, 5) > modelo(8)) {
			  matrizSaida(idia, 8) = (inicializacao(icenario, 5) - modelo(8)) * (1 - K_1ts); // Eq.14
			} else {
			  matrizSaida(idia, 8) = 0;
			}

			matrizSaida(idia, 9) = std::min(inicializacao(icenario, 5) - matrizSaida(idia, 8), modelo(5)) * (1 - K_2ts); // Eq.15

			matrizSaida(idia, 10) = inicializacao(icenario, 3) * (1 - K_3ts); // Eq.17

			matrizSaida(idia, 11) = std::max(inicializacao(icenario, 5) - modelo(5) - matrizSaida(idia, 8), 0.0) * (1 - K_2t2s); // Eq.16

			matrizSaida(idia, 12) = inicializacao(icenario, 6) * (1 - K_kts); // Eq.18

			// Calculation of state variables
			

			matrizSaida(idia, 1) = std::min(inicializacao(icenario, 4) + precipitacao(icenario, idia) - 
								matrizSaida(idia, 5) - matrizSaida(idia, 6) - matrizSaida(idia, 7),
								modelo(0));
			matrizSaida(idia, 4) = inicializacao(icenario, 6) + matrizSaida(idia, 7) - matrizSaida(idia, 12);
			Rsup_tmp = ((inicializacao(icenario, 4) + precipitacao(icenario, idia) - matrizSaida(idia, 5) - matrizSaida(idia, 6) - matrizSaida(idia, 7)) - modelo(0));
			matrizSaida(idia, 2) = inicializacao(icenario, 5) + matrizSaida(idia, 5) - matrizSaida(idia, 8) - matrizSaida(idia, 9) -
								   matrizSaida(idia, 11) + std::max(0.0, Rsup_tmp);
			matrizSaida(idia, 3) = std::max(inicializacao(icenario, 3) + matrizSaida(idia, 8) - matrizSaida(idia, 10) - Emarg(icenario, idia), 0.0);
			matrizSaida(idia, 14) = matrizSaida(idia, 9) * area / 86.4;
			matrizSaida(idia, 15) = matrizSaida(idia, 10) * area / 86.4;
			matrizSaida(idia, 16) = matrizSaida(idia, 11) * area / 86.4;
			matrizSaida(idia, 17) = matrizSaida(idia, 12) * area / 86.4;
			matrizSaida(idia, 0) = matrizSaida(idia, 14) + matrizSaida(idia, 15) + matrizSaida(idia, 16) + matrizSaida(idia, 17);
			
			inicializacao(icenario, 4) = matrizSaida(idia, 1);
			inicializacao(icenario, 5) = matrizSaida(idia, 2);
			inicializacao(icenario, 6) = matrizSaida(idia, 4);
			inicializacao(icenario, 3) = matrizSaida(idia, 3);
		}
		listaSaida(icenario) = matrizSaida;
	}
	return listaSaida;
}
