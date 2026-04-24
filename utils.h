#pragma once

#include <vector>
#include <Eigen/Dense>

#include "objects.h";

double rmse(std::vector<double>& result, std::vector<double>& truth);
double randNum();
double normal_pdf(double x, double mean, double stddev);
double normal_pdf_Nd(VectorXd x, VectorXd mean, MatrixXd sigma);
double naiveMonteCarloAverage(Funct& funct, double N);
double naiveMonteCarloAverageNd(NdFunct& funct, double N);
void printVector(std::vector<double> v);
void printVectorAsMatrix(std::vector<double> v);
VectorXd randVector(int dim);
bool outsideDomain(VectorXd x);