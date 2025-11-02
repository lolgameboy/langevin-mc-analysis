#pragma once

#include <vector>

#include "objects.h"

using std::vector;

vector<double> naiveMonteCarlo(PerlinFunct& funct, double N, int bins);
vector<double> langevinMonteCarlo(PerlinFunct& funct, double N, double stepsize, StdNormalSampler& stdSampler, int bins);
vector<double> multiLevelLangevinMonteCarlo(PerlinFunct& funct, double base_N, double max_N, double stepsize, StdNormalSampler& stdSampler, int bins);
vector<double> twoLevelLangevinMonteCarlo(PerlinFunct& funct, double base_N, double stepsize, StdNormalSampler& stdSampler, int bins);
vector<double> langevinMonteCarloBiased(PerlinFunct& funct, double N, double stepsize, StdNormalSampler& stdSampler, int bins);
vector<double> langevinMonteCarloRefinement(PerlinFunct& funct, double N, double stepsize, StdNormalSampler& stdSampler, int bins);
vector<double> naiveMultiLevelMonteCarlo(PerlinFunct& funct, double base_N, double max_N, double factor, int bins);