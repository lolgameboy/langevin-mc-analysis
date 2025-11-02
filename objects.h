#pragma once

#include <random>

#include "PerlinNoise.hpp"

class PerlinFunct {
public:
	siv::PerlinNoise perlin{ 0 };
	double frequency = 1;
	int octaves = 1;

	double value(double x);
	double gradientlog(double x);
};

class StdNormalSampler {
public:
	std::random_device rd{};
	std::mt19937 gen{ rd() };
	std::normal_distribution<double> d{ 0, 1 };

	double sample();
};