#pragma once

#include <vector>
#include <random>
#include <Eigen/Dense>

#include "PerlinNoise.hpp"

using Eigen::MatrixXd;
using Eigen::VectorXd;

class Funct {
public:
	virtual double value(double x) = 0;
	double gradientlog(double x);
};

class PerlinFunct : public Funct {
public:
	siv::PerlinNoise perlin{ 0 };
	double frequency = 1;
	int octaves = 1;

	double value(double x);
	void setSeed(unsigned int s);
};

// TODO StepFunct gaat problemen krijgen met het logaritme
class StepFunct : public Funct {
public:
	int steps = 1;
	std::vector<double> values;

	StepFunct(int steps, std::vector<double> values);
	double value(double x);
};

class DiffFunct : public Funct {
private:
	Funct* funct1;
	Funct* funct2;

public:
	DiffFunct(Funct* funct1, Funct* funct2);

	double value(double x);
};

class StdNormalSampler {
public:
	std::random_device rd{};
	std::mt19937 gen{ rd() };
	std::normal_distribution<double> d{ 0, 1 };

	double sample();
	VectorXd sampleVector(int dim);
};

class NdFunct {
public:
	int dim;

	bool limitDomain;

	NdFunct(int d, bool limdom) {
		dim = d;
		limitDomain = limdom;
	};

	virtual double value(VectorXd x) = 0;
	virtual VectorXd gradientlog(VectorXd x);
};

class GaussianNd : public NdFunct {
public:
	VectorXd mean;
	MatrixXd covar;
	MatrixXd covar_inv;

	GaussianNd(int d, MatrixXd cov_matrix, bool limitDom);

	double value(VectorXd x);
	virtual VectorXd gradientlog(VectorXd x);
};

class ComplexNd : public NdFunct {
public:
	double lambda;

	ComplexNd(int dim, int lam, bool limitDom);

	double value(VectorXd x);
};
