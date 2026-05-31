#include "objects.h"

void plotConvergenceNd(
	int dim,
	double cov_factor,
	double stepsize, 
	int pixelSideSize, 
	double referenceN, 
	double beta, 
	double delta, 
	double diminAdapt1, 
	double diminAdapt2, 
	bool plotPrecond = false, 
	bool skipLowNCounts = false,
	double anisotropy_factor = 1,
	bool plotRefLine = true
);

void plotExperimentCorrelation(int dim, double cov_factor, double stepsize, int pixelSideSize, double referenceN, bool skipLowNCounts = true);