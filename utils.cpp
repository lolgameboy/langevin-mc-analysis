#include <random>
#include <vector>
#include <cmath>
#include <iostream>

#include "objects.h"
#include "utils.h"
#include <Eigen/Dense>

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

double normal_pdf_Nd(VectorXd x, VectorXd mean, MatrixXd sigma) {
	//std::cout << "x: " << x << "\n";
	//std::cout << "mean: " << mean << "\n";
	//std::cout << "sigma: " << sigma << "\n";
	//std::cout << "result: " << pow(2 * EIGEN_PI, -sigma.diagonalSize() / 2) * pow(sigma.determinant(), -0.5) * exp(-0.5 * (x - mean).transpose() * sigma.inverse() * (x - mean)) << "\n";
	return pow(2 * EIGEN_PI, -sigma.diagonalSize() / 2) * pow(sigma.determinant(), -0.5) * exp(-0.5 * (x - mean).transpose() * sigma.inverse() * (x - mean));
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

double naiveMonteCarloAverageNd(NdFunct& funct, double N) {
	double expectedValue = 0;
	for (long long x = 0; x < N; x++) {
		expectedValue += funct.value(randVector(funct.dim)) / N;
	}
	return expectedValue;
}

void printVector(vector<double> v) {
	for (int i = 0; i < v.size(); i++) {
		std::cout << v[i] << " |";
	}
}

void printVectorAsMatrix(vector<double> v) {
	int side = sqrt(v.size());
	for (double i = 0; i < side; i++) {
		std::cout << "|";
		for (double j = 0; j < side; j++) {
			std::cout << v[i*side + j] << "|";
		}
		std::cout << "\n";
	}
}

VectorXd randVector(int dim) {
	VectorXd v = VectorXd(dim);
	for (int i = 0; i < dim; i++) {
		v(i) = randNum();
	}
	return v;
}

bool outsideDomain(VectorXd x) {
	for (int i = 0; i < x.size(); i++) {
		if (x(i) < 0 || x(i) > 1) {
			return true;
		}
	}
	return false;
}
