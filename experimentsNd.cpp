#include <vector>
#include <cmath>

#include "matplot/matplot.h"

#include <iostream>
#include "objects.h"
#include "utils.h"
#include "algorithmsNd.h"
#include "experimentsNd.h"

using std::vector;
using std::pow;
using std::sqrt;
using std::cout;

void plotConvergenceNd(int dim, double cov_factor, double stepsize, int pixelSideSize, double referenceN, double beta, double delta, double diminAdapt1, double diminAdapt2, bool plotPrecond, bool skipLowNCounts, double anisotropy_factor, bool plotRefLine) {
	// Create gaussian
	MatrixXd covar(dim, dim);
	for (int i = 0; i < dim; i++) {
		for (int j = 0; j < dim; j++) {
			covar(i, j) = pow(cov_factor, abs(i - j));
			// TODO? Not a parameter of experiment yet, but this is 'gaussian_width'. Essentially makes the gaussian wider or slimmer
			// Still remains a valid distribution though!
			// 0.1 is more difficult than 1 for naive MC because gaussian is big in center of screen and small at edge. 
			// A wider gaussian is more uniform over the entire domain.
			covar(i, j) *= 0.1;
			// Anisotropy
			if (i == 0) covar(i, j) *= anisotropy_factor;
		}
	}
	GaussianNd funct(dim, covar, true);
	
	vector<double> referencePixels = naiveMonteCarloNd(funct, referenceN, pixelSideSize);
	
	// CAREFUL! Normally, you multiply LMC by L to get actual function
	// Now, I turn it around and scale naive MC to integrate to 1. This makes it easier to compare between functions.
	// (i.e. it fixes the shifting issue from the logbook week 8 experiments)
	double L = naiveMonteCarloAverageNd(funct, referenceN);

	// SEE EXPLANATION ABOVE, REFERENCE PIXELS ALSO HAVE TO SCALE IN THIS APPROACH
	for (int j = 0; j < pixelSideSize * pixelSideSize; j++) {
		referencePixels[j] = referencePixels[j] / L;
	}

	int count = 0;
	while (pow(2, count) < referenceN / 10) {
		count++;
	}

	vector<double> Ns(count);
	for (double i = 0; i < count; i++)
	{
		Ns[i] = pow(2, i);
	}

	StdNormalSampler stdSampler;
	vector<double> errorsLMC(Ns.size());
	vector<double> errorsLMCPrecond(Ns.size());
	for (double i = 0; i < Ns.size(); i++)
	{
		double N = Ns[i];

		vector<double> pixels = langevinMonteCarloNd(funct, N, stepsize, stdSampler, pixelSideSize, beta, delta, false, false, diminAdapt1, diminAdapt2);
		vector<double> pixelsPrecond = langevinMonteCarloNd(funct, N, stepsize, stdSampler, pixelSideSize, beta, delta, true, false, diminAdapt1, diminAdapt2);
		cout << "\n " << L << "\n";
		// Langevin gives *proportional* distribution. Scale to get actual function.
		// UPDATE: Temporarily not anymore, see explanation above
		printVectorAsMatrix(pixels);
		// SEE EXPLANATION ABOVE
		//for (int j = 0; j < pixelSideSize * pixelSideSize; j++) {
		//	pixels[j] = pixels[j] * L;
		//	pixelsPrecond[j] = pixelsPrecond[j] * L;
		//}
		cout << "\n";

		errorsLMC[i] = rmse(pixels, referencePixels);
		errorsLMCPrecond[i] = rmse(pixelsPrecond, referencePixels);
	}

	vector<double> errorsMC(Ns.size());
	for (double i = 0; i < Ns.size(); i++)
	{
		double N = Ns[i];

		vector<double> pixels = naiveMonteCarloNd(funct, N, pixelSideSize);
		// SEE EXPLANATION ABOVE
		for (int j = 0; j < pixelSideSize * pixelSideSize; j++) {
			pixels[j] = pixels[j] / L;
		}
		errorsMC[i] = rmse(pixels, referencePixels);
	}

	vector<double> sqrtN(Ns.size());
	for (double i = 0; i < Ns.size(); i++) {
		sqrtN[i] = 1 / sqrt(Ns[i]);
	}

	matplot::figure()->size(1600, 1200);
	matplot::loglog(Ns, errorsMC)->display_name("Monte Carlo");
	matplot::hold("on");
	matplot::loglog(Ns, errorsLMC)->display_name("Langevin Monte Carlo");
	if (plotPrecond)
		matplot::loglog(Ns, errorsLMCPrecond)->display_name("Preconditioned Langevin Monte Carlo");
	//matplot::loglog(Ns, errorsMLMC)->display_name("Naive Multi Level Monte Carlo");
	//matplot::loglog(Ns, errorsHybrid)->display_name("Two Level Hybrid Monte Carlo");
	if (plotRefLine) matplot::loglog(Ns, sqrtN)->display_name("1/sqrt(N) reference");
	if (skipLowNCounts)
		matplot::xlim({ pow(pixelSideSize, 2) * 10, Ns[count - 1] });
	else 
		matplot::xlim({ 1, Ns[count - 1] });
	matplot::xlabel("N Samples in markov chain");
	matplot::ylabel("RMSE");
	std::string name = std::format("Convergence {}D-LMC, rho={}, stepsize={}, pixels={}x{}, beta={}", funct.dim, cov_factor, stepsize, pixelSideSize, pixelSideSize, beta);
	if (anisotropy_factor != 1) name += std::format(", gamma = {}", anisotropy_factor);
	if (diminAdapt1 != 0) name += std::format(", adaptationCoeff={}", diminAdapt1);
	matplot::title(name);
	matplot::legend();
	// Saving does not play nice. I save manually for now.
	//matplot::save("../../../figures/" + name + ".png");
	matplot::show();
}



void plotExperimentCorrelation(int dim, double cov_factor, double stepsize, int pixelSideSize, double referenceN, bool skipLowNCounts) {
	// Create gaussian
	MatrixXd covar(dim, dim);
	for (int i = 0; i < dim; i++) {
		for (int j = 0; j < dim; j++) {
			covar(i, j) = pow(0, abs(i - j));
			covar(i, j) *= 0.1;
		}
	}
	GaussianNd funct(dim, covar, true);

	vector<double> referencePixels = naiveMonteCarloNd(funct, referenceN, pixelSideSize);

	// CAREFUL! Normally, you multiply LMC by L to get actual function
	// Now, I turn it around and scale naive MC to integrate to 1. This makes it easier to compare between functions.
	// (i.e. it fixes the shifting issue from the logbook week 8 experiments)
	double L = naiveMonteCarloAverageNd(funct, referenceN);

	// SEE EXPLANATION ABOVE, REFERENCE PIXELS ALSO HAVE TO SCALE IN THIS APPROACH
	for (int j = 0; j < pixelSideSize * pixelSideSize; j++) {
		referencePixels[j] = referencePixels[j] / L;
	}

	int count = 0;
	while (pow(2, count) < referenceN / 10) {
		count++;
	}

	vector<double> Ns(count);
	for (double i = 0; i < count; i++)
	{
		Ns[i] = pow(2, i);
	}

	StdNormalSampler stdSampler;
	vector<double> errorsLMC(Ns.size());
	for (double i = 0; i < Ns.size(); i++)
	{
		double N = Ns[i];

		// dummy values
		float beta = 0;
		float delta = 0;
		float diminAdapt1 = 0;
		float diminAdapt2 = 0;

		vector<double> pixels = langevinMonteCarloNd(funct, N, stepsize, stdSampler, pixelSideSize, beta, delta, false, false, diminAdapt1, diminAdapt2);

		errorsLMC[i] = rmse(pixels, referencePixels);
	}

	vector<double> errorsMC(Ns.size());
	for (double i = 0; i < Ns.size(); i++)
	{
		double N = Ns[i];

		vector<double> pixels = naiveMonteCarloNd(funct, N, pixelSideSize);
		// SEE EXPLANATION ABOVE
		for (int j = 0; j < pixelSideSize * pixelSideSize; j++) {
			pixels[j] = pixels[j] / L;
		}
		errorsMC[i] = rmse(pixels, referencePixels);
	}

	vector<double> sqrtN(Ns.size());
	for (double i = 0; i < Ns.size(); i++) {
		sqrtN[i] = 1 / sqrt(Ns[i]);
	}


	matplot::figure()->size(1600, 1200);
	matplot::loglog(Ns, errorsMC)->display_name("Monte Carlo (rho=0)");
	matplot::hold("on");
	matplot::loglog(Ns, errorsLMC)->display_name("Langevin Monte Carlo (rho=0)");
	matplot::hold("on");
	matplot::loglog(Ns, sqrtN)->display_name("1/sqrt(N) reference");
	matplot::hold("on");


	MatrixXd covarCorr(dim, dim);
	// Gaussian with correlation
	for (int i = 0; i < dim; i++) {
		for (int j = 0; j < dim; j++) {
			covarCorr(i, j) = pow(cov_factor, abs(i - j));
			covarCorr(i, j) *= 0.1;
		}
	}

	GaussianNd functCorr(dim, covarCorr, true);

	referencePixels = naiveMonteCarloNd(functCorr, referenceN, pixelSideSize);
	L = naiveMonteCarloAverageNd(functCorr, referenceN);
	for (int j = 0; j < pixelSideSize * pixelSideSize; j++) {
		referencePixels[j] = referencePixels[j] / L;
	}

	for (double i = 0; i < Ns.size(); i++)
	{
		double N = Ns[i];

		// dummy values
		float beta = 0;
		float delta = 0;
		float diminAdapt1 = 0;
		float diminAdapt2 = 0;

		vector<double> pixels = langevinMonteCarloNd(functCorr, N, stepsize, stdSampler, pixelSideSize, beta, delta, false, false, diminAdapt1, diminAdapt2);

		errorsLMC[i] = rmse(pixels, referencePixels);
	}

	for (double i = 0; i < Ns.size(); i++)
	{
		double N = Ns[i];

		vector<double> pixels = naiveMonteCarloNd(functCorr, N, pixelSideSize);
		// SEE EXPLANATION ABOVE
		for (int j = 0; j < pixelSideSize * pixelSideSize; j++) {
			pixels[j] = pixels[j] / L;
		}
		errorsMC[i] = rmse(pixels, referencePixels);
	}

	matplot::loglog(Ns, errorsMC)->display_name(std::format("Monte Carlo (rho={})", cov_factor));
	matplot::hold("on");
	matplot::loglog(Ns, errorsLMC)->display_name(std::format("Langevin Monte Carlo (rho={})", cov_factor));
	matplot::hold("on");
	if (skipLowNCounts)
		matplot::xlim({ pow(pixelSideSize, 2) * 10, Ns[count - 1] });
	else
		matplot::xlim({ 1, Ns[count - 1] });
	matplot::xlabel("N Samples in markov chain");
	matplot::ylabel("RMSE");
	std::string name = std::format("Convergence {}D-LMC, stepsize={}, pixels={}x{}", functCorr.dim, stepsize, pixelSideSize, pixelSideSize);
	matplot::title(name);
	matplot::legend();
	matplot::show();
}