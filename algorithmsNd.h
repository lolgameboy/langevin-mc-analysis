#include <vector>
#include "objects.h"

using std::vector;

vector<double> naiveMonteCarloNd(NdFunct& funct, double N, double pixelSideSize);
vector<double> langevinMonteCarloNd(NdFunct& funct, double N, double stepsize, StdNormalSampler& stdSampler, double pixelSideSize, double beta, double delta, bool precond, bool log=false, double diminAdapt1=1, double diminAdapt2=1);