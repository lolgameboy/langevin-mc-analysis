#pragma once

#include <vector>

#include "objects.h"

using std::vector;

double rmse(vector<double>& result, vector<double>& truth);
double randNum();
double normal_pdf(double x, double mean, double stddev);
double naiveMonteCarloAverage(Funct& funct, double N);
void printVector(vector<double> v);