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

void plotConvergence(PerlinFunct& funct, double stepsize, int bins, double referenceN) {
	vector<double> referencePixels = naiveMonteCarlo(funct, referenceN, bins);

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

		errorsLMC[i] = rmse(pixels, referencePixels, N);
	}

	vector<double> errorsMC(Ns.size());
	for (double i = 0; i < Ns.size(); i++)
	{
		double N = Ns[i];

		vector<double> pixels = naiveMonteCarlo(funct, N, bins);

		errorsMC[i] = rmse(pixels, referencePixels, N);
	}

	vector<double> errorsMLMC(Ns.size());
	for (double i = 0; i < Ns.size(); i++)
	{
		double N = Ns[i];

		vector<double> pixels = naiveMultiLevelMonteCarlo(funct, N / 1.8, N, 0.5, bins, false);

		errorsMLMC[i] = rmse(pixels, referencePixels, N);
	}

	vector<double> errorsHybrid(Ns.size());
	for (double i = 0; i < Ns.size(); i++)
	{
		double N = Ns[i];

		vector<double> pixels = twoLevelHybrid(funct, N, bins, stepsize, stdSampler);

		errorsHybrid[i] = rmse(pixels, referencePixels, N);
	}

	vector<double> sqrtN(Ns.size());
	for (double i = 0; i < Ns.size(); i++) {
		sqrtN[i] = 1/sqrt(Ns[i]);
	}

	vector<double> inverseN(Ns.size());
	for (double i = 0; i < Ns.size(); i++) {
		inverseN[i] = 1 / Ns[i];
	}

	matplot::figure()->size(1600, 1200);
	matplot::loglog(Ns, errorsMC)->display_name("Monte Carlo");
	matplot::hold("on");
	matplot::loglog(Ns, errorsLMC)->display_name("Langevin Monte Carlo");
	matplot::loglog(Ns, errorsMLMC)->display_name("Naive Multi Level Monte Carlo");
	matplot::loglog(Ns, errorsHybrid)->display_name("Two Level Hybrid Monte Carlo");
	matplot::loglog(Ns, sqrtN)->display_name("1/sqrt(N) reference");
	matplot::loglog(Ns, inverseN)->display_name("1/N reference");
	matplot::legend();
	matplot::show();
}