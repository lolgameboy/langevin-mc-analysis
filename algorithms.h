#pragma once

#include <vector>

#include "objects.h"

using std::vector;

vector<double> naiveMonteCarlo(Funct& funct, double N, int bins);
vector<double> langevinMonteCarlo(Funct& funct, double N, double stepsize, StdNormalSampler& stdSampler, int bins);
vector<double> multiLevelLangevinMonteCarlo(PerlinFunct& funct, double base_N, double max_N, double stepsize, StdNormalSampler& stdSampler, int bins);
vector<double> twoLevelLangevinMonteCarlo(PerlinFunct& funct, double base_N, double stepsize, StdNormalSampler& stdSampler, int bins);
vector<double> langevinMonteCarloBiased(PerlinFunct& funct, double N, double stepsize, StdNormalSampler& stdSampler, int bins);
vector<double> langevinMonteCarloRefinement(PerlinFunct& funct, double N, double stepsize, StdNormalSampler& stdSampler, int bins);
vector<double> naiveMultiLevelMonteCarlo(PerlinFunct& funct, double base_N, double max_N, double factor, int bins, bool log = false);
vector<double> twoLevelHybrid(Funct& funct, double N, double bins, double stepsize, StdNormalSampler& stdSampler);