	#pragma once

#include "ObjectiveFunction.h"
#include "GradientDescentFunctionMinimizer.h"

/**
	use Newton's method to optimize a function. p will store the final value that minimizes the function, and its initial value
	is used to start the optimization method.

	Task: find p that minimize f(p). This means that df/dp(p) = 0.
	df/dp(p+dp) ~ df/dp(p) + d/dp(df/dp) * dp = 0 ==> -df/dp(p) = d/dp(df/dp) * dp
	Iterating the above, will hopefully get p that minimizes f.
*/
class NewtonFunctionMinimizer : public GradientDescentFunctionMinimizer {
public:
	NewtonFunctionMinimizer(int maxIterations = 100, double solveResidual = 0.0001, int maxLineSearchIterations = 15)
		: GradientDescentFunctionMinimizer(maxIterations, solveResidual, maxLineSearchIterations) {	}

	virtual ~NewtonFunctionMinimizer() {}

protected:
	// The search direction is given by -Hinv * g
	virtual void computeSearchDirection(ObjectiveFunction *function, const VectorXd &x, VectorXd& dx) {

		// Ex 1.3
		Eigen::SparseMatrix<double> A(x.size(), x.size());
		A.setZero();
		function->addHessianEntriesTo(hessianEntries, x);

		A.setFromTriplets(hessianEntries.begin(), hessianEntries.end());

		function->addGradientTo(dx, x);
		dx = -1 * dx;

		if (H.rows() == 0)
			H = A;
		else 
			H += A;

		Eigen::SimplicialLDLT< Eigen::SparseMatrix<double>> solver;
		solver.compute(H);
		VectorXd p = solver.solve(dx);

		dx = p;

	}



public:
	SparseMatrixd H;
	std::vector<Tripletd> hessianEntries;
};
