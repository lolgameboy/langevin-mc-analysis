#include <cmath>
#include <vector>
#include <iostream>

#include "fmt/core.h"

#include "algorithms.h"
#include "objects.h"
#include "utils.h"

using std::vector;
using std::pow;
using std::sqrt;
using std::min;
using std::cout;



vector<double> naiveMonteCarlo(PerlinFunct& funct, double N, int bins) {
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

vector<double> langevinMonteCarlo(PerlinFunct& funct, double N, double stepsize, StdNormalSampler& stdSampler, int bins) {
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
		pixels[bin] += 1 / N;
	}
	return pixels;
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
		vector<double> correction = langevinMonteCarloRefinement(funct, base_N, base_stepsize * pow(2, -i), stdSampler, bins);
		for (double j = 0; j < bins; j++) {
			pixels[j] += correction[j] / base_N;
		}
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
		double u_next_C = u_C + 0.5 * stepsize * funct.gradientlog(u_C) + sqrt(stepsize / 2) * (w1 + w2);
		double u_next_F1 = u_F + 0.5 * stepsize / 2 * funct.gradientlog(u_F) + sqrt(stepsize / 2) * w1;
		double u_next_F2 = u_next_F1 + 0.5 * stepsize / 2 * funct.gradientlog(u_next_F1) + sqrt(stepsize / 2) * w2;

		//// Calculate transition probabilities (TODO can be simplified)
		//double T_u_un = normal_pdf(u_next_C, u_C + 0.5 * stepsize * funct.gradientlog(u_C), sqrt(stepsize));
		//double T_un_u = normal_pdf(u_C, u_next_C + 0.5 * stepsize * funct.gradientlog(u_next_C), sqrt(stepsize));

		//// Acceptance ratio
		//double a = min(1.0, (funct.value(u_next_C) / funct.value(u_C)) * (T_un_u / T_u_un));

		//double T_u_un_F1 = normal_pdf(u_next_F1, u_F + 0.5 * stepsize * funct.gradientlog(u_F), sqrt(stepsize));
		//double T_un_u_F1 = normal_pdf(u_F, u_next_F1 + 0.5 * stepsize * funct.gradientlog(u_next_F1), sqrt(stepsize));

		//// Acceptance ratio
		//double a_F1 = min(1.0, (funct.value(u_next_F1) / funct.value(u_F)) * (T_un_u_F1 / T_u_un_F1));

		//double T_u_un_F2 = normal_pdf(u_next_F2, u_next_F1 + 0.5 * stepsize * funct.gradientlog(u_next_F1), sqrt(stepsize));
		//double T_un_u_F2 = normal_pdf(u_next_F1, u_next_F2 + 0.5 * stepsize * funct.gradientlog(u_next_F2), sqrt(stepsize));

		//// Acceptance ratio
		//double a_F2 = min(1.0, (funct.value(u_next_F2) / funct.value(u_next_F1)) * (T_un_u_F2 / T_u_un_F2));

		//if (randNum() < a) {
		//	u_C = u_next_C;
		//}

		//if (randNum() < a_F1) {
		//	if (randNum() < a_F2) {
		//		u_F = u_next_F2;
		//	}
		//	u_F = u_next_F1;
		//}

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

		u_C = u_next_C;
		u_F = u_next_F2;

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

vector<double> naiveMultiLevelMonteCarlo(PerlinFunct& funct, double base_N, double max_N, double factor, int bins) {
	vector<double> pixels(bins);
	vector<double> correction(bins);


	double total_N = 0;
	for (int i = 0; total_N < max_N; i++) {
		double N = base_N * pow(factor, i);

		for (long long x = 0; x < N; x++) {
			double u = randNum();
			int bin = floor(u * bins);
			if (bin == bins) {
				bin = bins - 1;
			}
			correction[bin] += (funct.value(u) - pixels[bin]) / N;
		}

		for (double j = 0; j < bins; j++) {
			cout << "| " << correction[j];
		}
		cout << "\n\n";

		for (int j = 0; j < pixels.size(); j++) {
			pixels[j] = (total_N / max_N) * pixels[j] + (1 - (total_N / max_N)) * correction[j];
		}

		total_N += N;
	}
	return pixels;
}