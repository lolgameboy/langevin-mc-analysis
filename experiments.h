#pragma once

#include "objects.h"

void plotConvergence(PerlinFunct& funct, double stepsize, int bins, double referenceN = 10000000);
void plotMonteCarlo(Funct& funct, double N, int bins);
void plotLangevin(Funct& funct, double N, int bins, double stepsize);
void plotFunct(Funct& funct, int pointCount, matplot::figure_handle fig);