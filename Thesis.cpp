#include <iostream>
#include <cstdlib>
#include <cmath>
#include <limits>
#include <random>
#include "fmt/core.h"
#include "Thesis.h"
#include "PerlinNoise.hpp"

using namespace std;



const double REF_VALUE_1_1_SEED_0 = 0.5469523775068327;
const double REF_VALUE_1_10_SEED_0 = 0.5469623057974313;

const vector<double> REF_IMAGE_4_10_SEED_0 = { 3.7015173673878336e-05, 0.000765472204080547, 0.012926500628362604, 0.049659201960078295, 0.21900159120798326, 0.32269015719264615, 0.11828525535777415, 0.012715394872559143, 0.0008362857837636259, 3.3935767680542344e-05 };
class PerlinFunct;
class StdNormalSampler;
double naiveMonteCarloAverage(PerlinFunct& funct, double N);
vector<double> multiLevelLangevinMonteCarlo(PerlinFunct& funct, double base_N, double max_N, double stepsize, StdNormalSampler& stdSampler, int bins);
vector<double> twoLevelLangevinMonteCarlo(PerlinFunct& funct, double base_N, double stepsize, StdNormalSampler& stdSampler, int bins);
vector<double> naiveMonteCarlo(PerlinFunct& funct, double N, int bins);
vector<double> langevinMonteCarlo(PerlinFunct& funct, double N, double stepsize, StdNormalSampler& stdSampler, int bins);
vector<double> langevinMonteCarloBiased(PerlinFunct& funct, double N, double stepsize, StdNormalSampler& stdSampler, int bins);
vector<double> langevinMonteCarloRefinement(PerlinFunct& funct, double N, double stepsize, StdNormalSampler& stdSampler, int bins);
double randNum();
double normal_pdf(double x, double mean, double stddev);



class PerlinFunct {
public:
	siv::PerlinNoise perlin{ 0 };
	double frequency = 1;
	int octaves = 1;

	double value(double x) {
		// What to do with domain? Truncation makes gradient 0. 
		// Current fix: just let chain explore outside domain but do not count it to buckets.
		//if (x < 0 || x > 1) {
		//	return 0;
		//}
		double val = perlin.octave1D_01((x * frequency), octaves);
		val *= 2 * normal_pdf(x, 0.5, 0.1);
		return val;
	}

	double gradientlog(double x) {
		double grad = (log(value(x + 0.0000000001))
			- log(value(x - 0.0000000001))) / 0.0000000002;
		// Truncate gradients
		//if (grad > 10) {
		//	grad = 10;
		//}
		//else if (grad < -10) {
		//	grad = -10;
		//}

		return grad;
	}
};

class StdNormalSampler {
public:
	// Standard normal number generator
	std::random_device rd{};
	std::mt19937 gen{ rd() };
	std::normal_distribution<double> d{ 0, 1 };

	double sample() {
		return d(gen);
	}
};

int main()
{
	// Perlin funct object
	PerlinFunct funct;
	funct.frequency = 4;
	funct.octaves = 10;
	const int bins = 10;
	const double stepsize = 0.01;

	StdNormalSampler stdSampler;

	srand(0);

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
		for (double x = 11; x < 18; x++)
		{
			int N = pow(2, x);
			std::cout << "N = " << N << ":\n\tMonte Carlo:\n\t";
			vector<double> pixels = naiveMonteCarlo(funct, N, bins);
			double relError = 0;
			for (double j = 0; j < bins; j++) {
				cout << "| " << pixels[j];
				//fmt::print("| {}", pixels[j]);
				//cout << "| " << abs(pixels[j] - prevPixels[j]);
				relError += abs((pixels[j] - REF_IMAGE_4_10_SEED_0[j]) / REF_IMAGE_4_10_SEED_0[j]);
			}
			std::cout << "\nRelative Error: " << relError;
			std::cout << "\n\tLangevin Monte Carlo:\n\t";
			relError = 0;
			pixels = langevinMonteCarlo(funct, N, stepsize, stdSampler, bins);
			for (double j = 0; j < bins; j++) {
				// Normalize pixel values
				pixels[j] /= N;
			}
			for (double j = 0; j < bins; j++) {
				cout << "| " << pixels[j];
				//fmt::print("| {}", pixels[j]);
				relError += abs((pixels[j] - REF_IMAGE_4_10_SEED_0[j]) / REF_IMAGE_4_10_SEED_0[j]);
			}
			std::cout << "\nRelative Error: " << relError;
			relError = 0;
			std::cout << "\n\Two Level Langevin Monte Carlo:\n\t";
			pixels = twoLevelLangevinMonteCarlo(funct, N / 3, stepsize, stdSampler, bins);
			// No normalization needed, already present in above method.
			for (double j = 0; j < bins; j++) {
				cout << "| " << pixels[j];
				//fmt::print("| {}", pixels[j]);
				//cout << "| " << abs(pixels[j] - prevPixels[j]);
				relError += abs((pixels[j] - REF_IMAGE_4_10_SEED_0[j]) / REF_IMAGE_4_10_SEED_0[j]);
			}
			std::cout << "\nRelative Error: " << relError;
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
	if (true) {
		for (double x = 5; x < 6; x++) {
			double N = pow(10, x);
			double relError = 0;
			std::cout << "Multi Level Langevin Monte Carlo:\nBase:\n";

			vector<double> pixels = langevinMonteCarlo(funct, N, 0.01, stdSampler, bins);
			for (double j = 0; j < bins; j++) {
				// Normalize pixel values
				pixels[j] /= N;
				cout << "| " << pixels[j];
			}
			for (int i = 0; i < 50; i++) {
				fmt::print("\nLevel {} correction:\n", i);
				vector<double> correction = langevinMonteCarloRefinement(funct, N, 0.01 * pow(2.0, -i), stdSampler, bins);
				for (double j = 0; j < bins; j++) {
					// Normalize correction values
					correction[j] /= N;
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
		}
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

double normal_pdf(double x, double mean = 0.0, double stddev = 1.0) {
	const double pi = std::acos(-1);
	double exponent = -0.5 * std::pow((x - mean) / stddev, 2);
	double coefficient = 1.0 / (stddev * std::sqrt(2.0 * pi));
	return coefficient * std::exp(exponent);
}

vector<double> twoLevelLangevinMonteCarlo(PerlinFunct& funct, double base_N, double base_stepsize, StdNormalSampler& stdSampler, int bins) {
	vector<double> pixels = langevinMonteCarloBiased(funct, base_N, base_stepsize, stdSampler, bins);

	for (double j = 0; j < bins; j++) {
		// Normalize pixel values
		pixels[j] /= base_N;
	}

	vector<double> correction = langevinMonteCarloRefinement(funct, base_N, base_stepsize * 0.5, stdSampler, bins);
	for (double j = 0; j < bins; j++) {
		pixels[j] += correction[j] / base_N;
	}

	return pixels;
}

vector<double> multiLevelLangevinMonteCarlo(PerlinFunct& funct, double base_N, double max_N, double base_stepsize, StdNormalSampler& stdSampler, int bins) {
	vector<double> pixels = langevinMonteCarloBiased(funct, base_N, base_stepsize, stdSampler, bins);

	for (double j = 0; j < bins; j++) {
		// Normalize pixel values
		pixels[j] /= base_N;
	}

	for (int i = 0; i * base_N < max_N; i++) {
		vector<double> correction = langevinMonteCarloRefinement(funct, i * base_N, base_stepsize * pow(2, -i), stdSampler, bins);
		for (double j = 0; j < bins; j++) {
			pixels[j] += correction[j] / (i * base_N);
		}
	}

	return pixels;
}

vector<double> langevinMonteCarlo(PerlinFunct &funct, double N, double stepsize, StdNormalSampler &stdSampler, int bins) {
	// initial random sample
	double u = 0.5;

	vector<double> pixels(bins);

	for (double i = 0; i < N; i++) {
		// Sample brownian term
		double w = stdSampler.sample();
		double u_next = u + 0.5 * stepsize * funct.gradientlog(u) + sqrt(stepsize) * w;

		// Calculate transition probabilities (TODO can be simplified)
		double T_u_un = normal_pdf(u_next, u + 0.5 * stepsize * funct.gradientlog(u), sqrt(stepsize));
		double T_un_u = normal_pdf(u, u_next + 0.5 * stepsize * funct.gradientlog(u_next), sqrt(stepsize));

		// Acceptance ratio
		double a = min(1.0, (funct.value(u_next) / funct.value(u)) * (T_un_u / T_u_un));

		if (randNum() < a) {
			u = u_next;
		}

		int bin = floor(u * bins);
		if (bin == bins) {
			bin = bins - 1;
		}
		if (bin < 0) bin = 0;
		pixels[bin] += 1;
	}
	return pixels;
}

vector<double> langevinMonteCarloBiased(PerlinFunct& funct, double N, double stepsize, StdNormalSampler& stdSampler, int bins) {
	// Initial random sample
	double u = 0.5;

	vector<double> pixels(bins);

	for (double i = 0; i < N; i++) {
		// Sample brownian term
		double w = stdSampler.sample();
		double u_next = u + 0.5 * stepsize * funct.gradientlog(u) + sqrt(stepsize) * w;

		u = u_next;

		int bin = floor(u * bins);
		// Only count if valid bin
		if (bin >= bins || bin < 0) {

		}
		else {
			pixels[bin] += 1;
		}
	}
	return pixels;
}

vector<double> langevinMonteCarloRefinement(PerlinFunct& funct, double N, double stepsize, StdNormalSampler& stdSampler, int bins) {
	// initial random sample
	double u_C = 0.5;
	double u_F = u_C;

	vector<double> pixels(bins);

	for (double i = 0; i < N; i++) {
		// Sample brownian term
		double w1 = stdSampler.sample();
		double w2 = stdSampler.sample();
		double u_next_C = u_C + 0.5 * stepsize * funct.gradientlog(u_C) + sqrt(stepsize) * (w1 + w2);
		double u_next_F1 = u_F + 0.5 * stepsize / 2 * funct.gradientlog(u_F) + sqrt(stepsize/2) * w1;
		double u_next_F2 = u_next_F1 + 0.5 * stepsize / 2 * funct.gradientlog(u_next_F1) + sqrt(stepsize / 2) * w2;
		
		// Calculate transition probabilities (TODO can be simplified)
		double T_u_un = normal_pdf(u_next_C, u_C + 0.5 * stepsize * funct.gradientlog(u_C), sqrt(stepsize));
		double T_un_u = normal_pdf(u_C, u_next_C + 0.5 * stepsize * funct.gradientlog(u_next_C), sqrt(stepsize));

		// Acceptance ratio
		double a = min(1.0, (funct.value(u_next_C) / funct.value(u_C)) * (T_un_u / T_u_un));

		double T_u_un_F1 = normal_pdf(u_next_F1, u_F + 0.5 * stepsize * funct.gradientlog(u_F), sqrt(stepsize));
		double T_un_u_F1 = normal_pdf(u_F, u_next_F1 + 0.5 * stepsize * funct.gradientlog(u_next_F1), sqrt(stepsize));

		// Acceptance ratio
		double a_F1 = min(1.0, (funct.value(u_next_F1) / funct.value(u_F)) * (T_un_u_F1 / T_u_un_F1));

		double T_u_un_F2 = normal_pdf(u_next_F2, u_next_F1 + 0.5 * stepsize * funct.gradientlog(u_next_F1), sqrt(stepsize));
		double T_un_u_F2 = normal_pdf(u_next_F1, u_next_F2 + 0.5 * stepsize * funct.gradientlog(u_next_F2), sqrt(stepsize));

		// Acceptance ratio
		double a_F2 = min(1.0, (funct.value(u_next_F2) / funct.value(u_next_F1)) * (T_un_u_F2 / T_u_un_F2));

		if (randNum() < a) {
			u_C = u_next_C;
		}

		if (randNum() < a_F1) {
			if (randNum() < a_F2) {
				u_F = u_next_F2;
			}
			u_F = u_next_F1;
		}


		int binC = floor(u_next_C * bins);
		int binF = floor(u_next_F2 * bins);

		if (binC != binF) {
			if (!(binC < 0 || binC >= bins)) {
				pixels[binC] -= 1;
			}
			if (!(binF < 0 || binF >= bins)) {
				pixels[binF] += 1;
			}
		}

		u_F = u_C;

		//double interp = u_C + ((u_next_C - u_C) / 2);

		//binC = floor(interp * bins);
		//if (binC >= bins) {
		//	binC = bins - 1;
		//}
		//if (binC < 0) binC = 0;

		//binF = floor(u_next_F1 * bins);
		//if (binF >= bins) {
		//	binF = bins - 1;
		//}
		//if (binF < 0) binF = 0;

		//if (binC != binF) {
		//	pixels[binC] -= L;
		//	pixels[binF] += L;
		//}
	}
	return pixels;
}

vector<double> naiveMonteCarlo(PerlinFunct &funct, double N, int bins) {
	vector<double> pixels(bins);

	for (long long x = 0; x < N; x++) {
		double u = randNum();
		int bin = floor(u * bins);
		if (bin == bins) {
			bin = bins - 1;
		}
		pixels[bin] += funct.value(u) / N;
	}

	return pixels;
}

double naiveMonteCarloAverage(PerlinFunct& funct, double N) {
	double expectedValue = 0;
	for (long long x = 0; x < N; x++) {
		expectedValue += funct.value(randNum()) / N;
	}
	return expectedValue;
}

double randNum() {
	float r = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
	return static_cast <double> (r);
}