#include <vector>
#include <cmath>

#include "fmt/core.h"
#include "matplot/matplot.h"

#include "objects.h"
#include "utils.h"
#include "algorithms.h"
#include "experiments.h"

using std::vector;
using std::pow;
using std::sqrt;
using std::cout;

void plotConvergence(PerlinFunct& funct, double stepsize, int bins, double referenceN) {
	vector<double> referencePixels = naiveMonteCarlo(funct, referenceN, bins);
	double L = naiveMonteCarloAverage(funct, referenceN);

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

		vector<double> pixels = langevinMonteCarlo(funct, N, stepsize, stdSampler, bins);
		cout << "\n " << L << "\n";
		// Langevin gives *proportional* distribution. Scale to get actual function.
		for (int j = 0; j < bins; j++) {
			pixels[j] = pixels[j] * L;
			cout << pixels[j] << " | ";
		}
		cout << "\n";

		errorsLMC[i] = rmse(pixels, referencePixels);
	}

	vector<double> errorsMC(Ns.size());
	for (double i = 0; i < Ns.size(); i++)
	{
		double N = Ns[i];

		vector<double> pixels = naiveMonteCarlo(funct, N, bins);

		errorsMC[i] = rmse(pixels, referencePixels);
	}

	vector<double> errorsMLMC(Ns.size());
	for (double i = 0; i < Ns.size(); i++)
	{
		double N = Ns[i];

		//vector<double> pixels = naiveMultiLevelMonteCarlo(funct, N / 1.8, N, 0.5, bins, false);

		//errorsMLMC[i] = rmse(pixels, referencePixels);
	}

	vector<double> errorsHybrid(Ns.size());
	for (double i = 0; i < Ns.size(); i++)
	{
		double N = Ns[i];

		//vector<double> pixels = twoLevelHybrid(funct, N, bins, stepsize, stdSampler);

		// errorsHybrid[i] = rmse(pixels, referencePixels);
	}

	vector<double> sqrtN(Ns.size());
	for (double i = 0; i < Ns.size(); i++) {
		sqrtN[i] = 1/sqrt(Ns[i]);
	}

	matplot::figure()->size(1600, 1200);
	matplot::loglog(Ns, errorsMC)->display_name("Monte Carlo");
	matplot::hold("on");
	matplot::loglog(Ns, errorsLMC)->display_name("Langevin Monte Carlo");
	//matplot::loglog(Ns, errorsMLMC)->display_name("Naive Multi Level Monte Carlo");
	//matplot::loglog(Ns, errorsHybrid)->display_name("Two Level Hybrid Monte Carlo");
	matplot::loglog(Ns, sqrtN)->display_name("1/sqrt(N) reference");
	matplot::legend();
	matplot::show();
}

void plotMonteCarlo(Funct &funct, double N, int bins) {
	vector<double> x(bins);

	for (double i = 0; i < bins; i++) {
		x[i] = i / bins;
	}

	vector<double> pixels = naiveMonteCarlo(funct, N, bins);
	printVector(pixels);
	std::cout << "\n";

	matplot::figure_handle fig = matplot::figure();
	fig->size(1600, 1200);
	plotFunct(funct, (double)bins * 10, fig);

	matplot::plot(x, pixels);
	matplot::legend();
	matplot::show();
}

void plotLangevin(Funct& funct, double N, int bins, double stepsize) {
	StdNormalSampler stdSampler;
	vector<double> pixels = langevinMonteCarlo(funct, N, 0.01, stdSampler, bins);

	printVector(pixels);
	std::cout << "\n";

	// Calculate L to correctly assess accuracy.
	double L = naiveMonteCarloAverage(funct, N);
	vector<double> x(bins);
	vector<double> scaled_pixels(bins);

	for (double i = 0; i < bins; i++) {
		x[i] = i / bins;
		scaled_pixels[i] = pixels[i] * L;
	}

	// Log the integral
	double integral = 0;
	for (double i = 0; i < bins; i++) {
		// Should intergrate to 1
		integral += pixels[i] / bins;
		std::cout << integral;
		std::cout << "\n";
	}

	matplot::figure_handle fig = matplot::figure();
	fig->size(1600, 1200);
	plotFunct(funct, (double)bins * 10, fig);
	matplot::stairs(x, pixels);
	matplot::stairs(x, scaled_pixels);
	matplot::legend();
	matplot::show();
}

void plotFunct(Funct& funct, int pointCount, matplot::figure_handle fig) {
	vector<double> func_y(pointCount);
	vector<double> func_x(pointCount);
	
	for (double i = 0; i < func_x.size(); i++) {
		func_x[i] = i / func_x.size();
		func_y[i] = funct.value(func_x[i]);
	}
	
	matplot::plot(func_x, func_y)->display_name("Target function");
	matplot::hold("on");
}