#include <cmath>
#include <random>
#include <vector>

#include "PerlinNoise.hpp"

#include "objects.h"
#include "utils.h"

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
	val *= 2 * normal_pdf(x, 0.5, 0.1);
	return val;
}

// Todo domein functie en gradients ivm stepfunct eens herbekijken
double Funct::gradientlog(double x) {
	// Domain is only defined between 0 and 1
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