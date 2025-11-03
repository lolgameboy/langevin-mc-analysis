#pragma once

#include <vector>
#include <random>

#include "PerlinNoise.hpp"

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
};