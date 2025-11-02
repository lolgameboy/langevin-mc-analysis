#include <cmath>
#include <random>

#include "PerlinNoise.hpp"

#include "objects.h"
#include "utils.h"

double PerlinFunct::value(double x) {
	// What to do with domain? Truncation makes gradient 0. 
	// Current fix: just let chain explore outside domain but do not count it to buckets.
	//if (x < 0 || x > 1) {
	//	return 0;
	//}
	double val = perlin.octave1D_01((x * frequency), octaves);
	val *= 2 * normal_pdf(x, 0.5, 0.1);
	return val;
}

double PerlinFunct::gradientlog(double x) {
	double grad = (log(value(x + 0.0000000001))
		- log(value(x - 0.0000000001))) / 0.0000000002;
	// Truncate gradients
	//if (grad > 100) {
	//	grad = 100;
	//}
	//else if (grad < -100) {
	//	grad = -100;
	//}

	return grad;
}

double StdNormalSampler::sample() {
	return d(gen);
}