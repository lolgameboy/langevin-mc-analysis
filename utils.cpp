#include <random>
#include <vector>
#include <cmath>
#include <iostream>

#include "objects.h"
#include "utils.h"

using std::vector;

double randNum() {
	float r = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
	return static_cast <double> (r);
}

double normal_pdf(double x, double mean = 0.0, double stddev = 1.0) {
	const double pi = std::acos(-1);
	double exponent = -0.5 * std::pow((x - mean) / stddev, 2);
	double coefficient = 1.0 / (stddev * std::sqrt(2.0 * pi));
	return coefficient * std::exp(exponent);
}

double rmse(vector<double>& result, vector<double>& truth) {
	double error = 0;
	double K = result.size();
	for (int i = 0; i < K; i++) {
		error += std::pow(result[i] - truth[i], 2);
	}
	error /= K;
	error = sqrt(error);
	return error;
}

double naiveMonteCarloAverage(Funct& funct, double N) {
	double expectedValue = 0;
	for (long long x = 0; x < N; x++) {
		expectedValue += funct.value(randNum()) / N;
	}
	return expectedValue;
}

void printVector(vector<double> v) {
	for (int i = 0; i < v.size(); i++) {
		std::cout << v[i] << " |";
	}
}