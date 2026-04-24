#include <cmath>
#include <vector>
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Cholesky>

#include <unsupported/Eigen/MatrixFunctions>

#include "algorithmsNd.h"
#include "objects.h"
#include "utils.h"

using std::vector;

vector<double> naiveMonteCarloNd(NdFunct& funct, double N, double pixelSideSize) {
	vector<double> pixels(pixelSideSize * pixelSideSize);
	vector<double> sampleCounts(pixelSideSize*pixelSideSize);

	for (long long x = 0; x < N; x++) {
		VectorXd u = randVector(funct.dim);
		int xbin = floor(u(0) * pixelSideSize);
		int ybin = floor(u(1) * pixelSideSize);
		if (xbin == pixelSideSize) { xbin = pixelSideSize - 1; }
		if (ybin == pixelSideSize) { ybin = pixelSideSize - 1; }

		int bin = ybin * pixelSideSize + xbin;
		pixels[bin] += funct.value(u);
		sampleCounts[bin] += 1;
	}

	for (int j = 0; j < pixelSideSize * pixelSideSize; j++) {
		if (sampleCounts[j] > 0) pixels[j] /= sampleCounts[j];
	}

	return pixels;
}


vector<double> langevinMonteCarloNd(NdFunct& funct, double N, double stepsize, StdNormalSampler& stdSampler, double pixelSideSize, double beta, double delta, bool precond, bool log, double diminAdapt1, double diminAdapt2) {
	// initial random uniform sample
	VectorXd u = randVector(funct.dim);

	vector<double> pixels(pixelSideSize*pixelSideSize);

	MatrixXd G = MatrixXd(funct.dim, funct.dim).setZero();

	for (double i = 0; i < N; i++) {
		// Calculate preconditioning matrix
		MatrixXd diag = (funct.gradientlog(u) * funct.gradientlog(u).transpose()).diagonal().asDiagonal().toDenseMatrix();
		G = beta * G + (1 - beta) * diag;								
		// OMG, bij i=0 geeft dit probleem
		double adaptCoeff = (1 / pow(i + 1, diminAdapt1));
		MatrixXd M = ((1 - adaptCoeff) * MatrixXd::Identity(funct.dim, funct.dim) + adaptCoeff * G.llt().matrixL().toDenseMatrix()).inverse();

		if (!precond) {
			// Disable preconditioning
			M = MatrixXd::Identity(funct.dim, funct.dim);
		}
		// Sample brownian term
		VectorXd w = stdSampler.sampleVector(funct.dim);
		VectorXd u_next = u + 0.5 * stepsize * M * funct.gradientlog(u) + sqrt(stepsize) * M.llt().matrixL().toDenseMatrix() * w;
		
		//if u_next is outside domain, u stays the same.
		if (!outsideDomain(u_next)) {
			if (log) std::cout << "Next sample: " << u_next << "\n";

			// Calculate transition probabilities (TODO can be simplified)
			double T_u_un = normal_pdf_Nd(u_next, u + 0.5 * stepsize * M * funct.gradientlog(u), stepsize * M);
			if (log) std::cout << "T_u_un: " << T_u_un << "\n";
			double T_un_u = normal_pdf_Nd(u, u_next + 0.5 * stepsize * M * funct.gradientlog(u_next), stepsize * M);
			if (log) std::cout << "T_un_u: " << T_un_u << "\n";

			// Acceptance ratio
			double a = std::min(1.0, (funct.value(u_next) / funct.value(u)) * (T_un_u / T_u_un));
			if (log) std::cout << "Acceptance chance: " << a << "\n";

			if (randNum() < a) {
				u = u_next;
			}
		}
		int xbin = floor(u(0) * pixelSideSize);
		int ybin = floor(u(1) * pixelSideSize);
		if (xbin == pixelSideSize) { xbin = pixelSideSize - 1; }
		if (ybin == pixelSideSize) { ybin = pixelSideSize - 1; }

		int bin = ybin * pixelSideSize + xbin;

		pixels[bin] += 1 / N * pixelSideSize * pixelSideSize; // Important realization: * bins is nessecary to ensure result integrates to 1. 
																// Think of uniform dist from 0 to 1 with 10 bins as example.
		if(log) std::cout << "END OF ITER\n\n";
	}
	return pixels;
}