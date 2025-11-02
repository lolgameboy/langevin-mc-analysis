#pragma once

#include <vector>

#include "objects.h"

using std::vector;

double rmse(vector<double>& result, vector<double>& truth, double N);
double randNum();
double normal_pdf(double x, double mean, double stddev);
double naiveMonteCarloAverage(PerlinFunct& funct, double N);