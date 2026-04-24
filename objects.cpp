#include <cmath>
#include <random>
#include <vector>
#include <iostream>

#include <Eigen/Dense>

#include "PerlinNoise.hpp"

#include "objects.h"
#include "utils.h"

using Eigen::VectorXd;

using std::vector;
using std::max;
using std::min;

double PerlinFunct::value(double x) {
	// What to do with domain? Truncation makes gradient 0. 
	// Current fix: just let chain explore outside domain but do not count it to buckets.
	// Update: kinda tricky. When using metropolis, it should just be 0 because then it will never be accepted. 
	// Although that could make it rather inefficient if the gradient points towards the edge. 
	// Maybe boundary conditions help?
	// It was disabled, but I enabled it again for experiments with regular langevin.
	if (x < 0 || x > 1) {
		return 0;
	}
	double val = perlin.octave1D_01((x * frequency), octaves);
	//val *= normal_pdf(x, 0.5, 0.1);
	val *= normal_pdf(x, 0.5, 0.2);
	return val;
}

// Todo domein functie en gradients ivm stepfunct eens herbekijken
double Funct::gradientlog(double x) {
	// Domain is only defined between 0 and 1
	// Todo I think this is not ideal, because what happens when x = 1.1 for example? Then second log is still outside domain
	double grad = (log(value(min(1.0, x + 0.0000000001)))
		- log(value(max(0.0 , x - 0.0000000001)))) / 0.0000000002;
	// Truncate gradients
	if (grad > 100) {
		grad = 100;
	}
	else if (grad < -100) {
		grad = -100;
	}

	return grad;
}

void PerlinFunct::setSeed(unsigned int s) {
	perlin.reseed(s);
}

StepFunct::StepFunct(int s, vector<double> vs) {
	steps = s;
	values = vs;
}

// TODO StepFunct gaat problemen krijgen met het logaritme
// Oplossing: stepwise linear function (of vormen de pieken daarvan ook een probleem omdat die niet C1 continu zijn?
double StepFunct::value(double x) {
	int bin = floor(x * steps);
	if (bin == steps) {
		bin = steps - 1;
	}
	if (bin < 0 || bin > steps) return 0;
	return values[bin];
}

DiffFunct::DiffFunct(Funct* f1, Funct* f2) {
	funct1 = f1;
	funct2 = f2;
}

double DiffFunct::value(double x) {
	return funct1->value(x) - funct2->value(x);
}

double StdNormalSampler::sample() {
	return d(gen);
}

VectorXd StdNormalSampler::sampleVector(int dim) {
	VectorXd v = VectorXd(dim);
	for (int i = 0; i < dim; i++) {
		v(i) = sample();
	}
	return v;
}

GaussianNd::GaussianNd(int dim, MatrixXd cov_matrix, bool limitDom) : NdFunct(dim, limitDom) {
	mean = VectorXd(dim);
	covar = cov_matrix;
	covar_inv = MatrixXd(dim, dim);

	for (int i = 0; i < dim; i++) {
		mean(i) = 0.5;
	}

	covar_inv = covar.inverse();
}

double GaussianNd::value(VectorXd x) {
	// Set function to 0 outside domain
	if (limitDomain) {
		for (int i = 0; i < x.size(); i++) {
			if (x(i) < 0 || x(i) > 1) {
				return 0;
			}
		}
	}
	
	// PDF of multivariate gaussian (without constant factor) from wikipedia
	return exp(-0.5 * (x - mean).transpose() * covar_inv * (x - mean));
}

VectorXd GaussianNd::gradientlog(VectorXd x) {
	return -(covar_inv * x);
}

ComplexNd::ComplexNd(int dim, int lam, bool limdom) : NdFunct(dim, limdom) {
	lambda = lam;
}

double ComplexNd::value(VectorXd x) {
	// Set function to 0 outside domain
	if (limitDomain) {
		for (int i = 0; i < x.size(); i++) {
			if (x(i) < 0 || x(i) > 1) {
				return 0;
			}
		}
	}

	double sum = 0;
	for (int i = 0; i < dim; i++) {
		sum += (pow(x(i), 2) / 2) + lambda * cos(x(i));
	}
	
	return sum;
}

VectorXd NdFunct::gradientlog(VectorXd x) {
	VectorXd grad(x.size());

	// Careful! Log of 0 is undefined. This can cause a problem outside the domain where everything is 0. 
	// Currently solved in langevin by instantly rejecting proposal outside domain.
	for (int i = 0; i < x.size(); i++) {
		x(i) = x(i) + 0.0000000001;
		double val1 = log(value(x));
		x(i) = x(i) - 0.0000000002;
		double val2 = log(value(x));
		grad(i) = (val1 - val2) / 0.0000000002;
		x(i) = x(i) + 0.0000000001;

		// Truncate gradients
		//if (grad > 100) {
		//	grad = 100;
		//}
		//else if (grad < -100) {
		//	grad = -100;
		//}
	}

	return grad;
}


