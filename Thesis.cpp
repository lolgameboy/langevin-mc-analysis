#include <iostream>
#include <cstdlib>
#include <limits>
#include <vector>

#include "fmt/core.h"
#include "matplot/matplot.h"

#include "objects.h"
#include "utils.h"
#include "algorithms.h"
#include "experiments.h"
#include "Thesis.h"


using std::vector;
using std::pow;
using std::sqrt;
using std::cout;

const double REF_VALUE_1_1_SEED_0 = 0.5469523775068327;
const double REF_VALUE_1_10_SEED_0 = 0.5469623057974313;

vector<double> REF_IMAGE_4_10_SEED_0 = { 3.7015173673878336e-05, 0.000765472204080547, 0.012926500628362604, 0.049659201960078295, 0.21900159120798326, 0.32269015719264615, 0.11828525535777415, 0.012715394872559143, 0.0008362857837636259, 3.3935767680542344e-05 };



int main()
{
	// Perlin funct object
	PerlinFunct funct;
	funct.frequency = 4;
	funct.octaves = 10;
	const int bins = 10;
	const double stepsize = 0.01;

	StdNormalSampler stdSampler;

	srand(1);

	funct.frequency = 1;
	funct.octaves = 30;
	plotConvergence(funct, 0.01, 100, 100000);
	
	
	return 0;


	// Old experiments
	// 
	// 

	// Get ref value
	if (false) {
		std::cout << "N = " << 10000000 << ":\n\tMonte Carlo:\n\t";
		vector<double> pixels = naiveMonteCarlo(funct, 10000000, bins);
		for (double j = 0; j < bins; j++) {
			fmt::print(", {}", pixels[j]);
		}
		return 0;
	}

	for (double i = 0; i < 10; i++) {
		std::cout << i / 10 << ": " << funct.value(i / 10) << " ";
	}
	std::cout << "\n";
	for (double i = 0; i < 10; i++) {
		std::cout << i / 10 << ": " << funct.gradientlog(i / 10) << " ";
	}

	// Experiment 1
	if (false) {
		for (double x = 10; x < 20; x++)
		{
			int N = pow(2, x);
			std::cout << "N = " << N << ":\n\tMonte Carlo:\n\t";
			vector<double> pixels = naiveMonteCarlo(funct, N, bins);
			for (double j = 0; j < bins; j++) {
				//cout << "| " << pixels[j];
				//fmt::print("| {}", pixels[j]);
				//cout << "| " << abs(pixels[j] - prevPixels[j]);
				// Normalize pixel values
				pixels[j] /= N;
			}
			std::cout << "\nRMSE: " << rmse(pixels, REF_IMAGE_4_10_SEED_0, N);
			std::cout << "\n\tLangevin Monte Carlo:\n\t";
			pixels = langevinMonteCarlo(funct, N, stepsize, stdSampler, bins);
			for (double j = 0; j < bins; j++) {
				// Normalize pixel values
				pixels[j] /= N;
			}
			for (double j = 0; j < bins; j++) {
				//cout << "| " << pixels[j];
				//fmt::print("| {}", pixels[j]);
			}
			std::cout << "\nRMSE: " << rmse(pixels, REF_IMAGE_4_10_SEED_0, N);
			std::cout << "\n\Multi Level Langevin Monte Carlo:\n\t";
			pixels = multiLevelLangevinMonteCarlo(funct, N / 10, N, 0.02, stdSampler, bins);
			// No normalization needed, already present in above method.
			for (double j = 0; j < bins; j++) {
				//cout << "| " << pixels[j];
				//fmt::print("| {}", pixels[j]);
				//cout << "| " << abs(pixels[j] - prevPixels[j]);
			}
			std::cout << "\nRMSE: " << rmse(pixels, REF_IMAGE_4_10_SEED_0, N);
			fmt::print("\n");
		}
	}
	// Experiment 2
	if (false) {
		int Z = 10;
		vector<double> langevinErrors(10);
		vector<double> multiLevelErrors(10);
		for (int z = 0; z < Z; z++)
		{
			int N = 5000;
			double relError = 0;

			vector<double> pixels = langevinMonteCarlo(funct, N, stepsize, stdSampler, bins);
			for (double j = 0; j < bins; j++) {
				relError = abs((pixels[j] - REF_IMAGE_4_10_SEED_0[j]) / REF_IMAGE_4_10_SEED_0[j]);
			}
			langevinErrors[z] = relError;

			pixels = twoLevelLangevinMonteCarlo(funct, N, stepsize, stdSampler, bins);
			for (double j = 0; j < bins; j++) {
				relError = abs((pixels[j] - REF_IMAGE_4_10_SEED_0[j]) / REF_IMAGE_4_10_SEED_0[j]);
			}
			multiLevelErrors[z] = relError;
		}

		std::cout << "\n\tLangevin Monte Carlo:\n\t";
		for (double j = 0; j < Z; j++) {
			cout << " " << langevinErrors[j];
		}
		std::cout << "\n\Two Level Langevin Monte Carlo:\n\t";
		for (double j = 0; j < Z; j++) {
			cout << " " << multiLevelErrors[j];
		}
	}

	// Experiment 3
	if (false) {
		double N = 10000;
		double relError = 0;
		std::cout << "Multi Level Langevin Monte Carlo:\nBase:\n";

		vector<double> pixels = langevinMonteCarloBiased(funct, N, 0.02, stdSampler, bins);
		for (double j = 0; j < bins; j++) {
			// Normalize pixel values
			pixels[j] /= N;
			cout << "| " << pixels[j];
		}
		for (int i = 0; i < 50; i++) {
			fmt::print("\nLevel {} correction:\n", i);
			vector<double> correction = langevinMonteCarloRefinement(funct, N * pow(2, i), 0.02 * pow(2.0, -i), stdSampler, bins);
			for (double j = 0; j < bins; j++) {
				// Normalize correction values
				correction[j] /= N * pow(2, i);
				cout << "| " << correction[j];
				pixels[j] += correction[j];
			}

			fmt::print("\n\nFinal pixel values: ");
			for (double j = 0; j < bins; j++) {
				cout << "| " << pixels[j];
			}
		}

		fmt::print("\n\nFinal pixel values: ");
		for (double j = 0; j < bins; j++) {
			cout << "| " << pixels[j];
		}

		//std::cout << "\nRelative Error: " << relError;
		fmt::print("\n");
		std::cout << '\n';
	}
	
	// Experiment 4
	if (false) {
		double N = 10000;
		double relError = 0;
		std::cout << "\n\Two Level Langevin Monte Carlo:\n\t";
		vector<double> pixels = twoLevelLangevinMonteCarlo(funct, N / 3, stepsize, stdSampler, bins);
		for (double j = 0; j < bins; j++) {
			// Normalize pixel values
			pixels[j] /= N / 3;
		}
		for (double j = 0; j < bins; j++) {
			cout << "| " << pixels[j];
			//fmt::print("| {}", pixels[j]);
			//cout << "| " << abs(pixels[j] - prevPixels[j]);
			relError += abs((pixels[j] - REF_IMAGE_4_10_SEED_0[j]) / REF_IMAGE_4_10_SEED_0[j]);
		}
		std::cout << "\nRelative Error: " << relError;
	}
}