
#pragma once

#include "ObjectiveFunction.h"

class RosenbrockFunction : public ObjectiveFunction {
public:

	RosenbrockFunction() {
		a = 1; b = 100;
	}

	virtual double computeValue(const VectorXd& x) {

		// Ex 1.1
		// return f(x)
		double X = x[0], Y = x[1];
		return pow(a - X, 2) + b * pow(Y - pow(X, 2), 2);
	}

	virtual void addGradientTo(VectorXd& grad, const VectorXd& x) {

		// Ex 1.1
		// write df/dx in `grad`
		double X = x[0];
		double Y = x[1];
		grad[0] = -4 * b*X*(-X * X + Y) - 2 * (a - X);
		grad[1] = 2 * b*(Y - X * X);

	}

	virtual void addHessianEntriesTo(std::vector<Tripletd>& hessianEntries, const VectorXd& x) {

		// Ex 1.2
		// write d^2f/dx^2 in `hessianEntries`
		double X = x[0];
		double Y = x[1];
		double dxx, dxy, dyx, dyy;
		dxx = 2 + 4 * b*(-Y + 3 * X*X); // index 0,0
		dxy = -4 * b*X; // index 0,1
		dyx = -4 * b*X; // index 1,0
		dyy = 2 * b; // index 1,1

		hessianEntries.push_back(Tripletd(0, 0, dxx));
		hessianEntries.push_back(Tripletd(0, 1, dxy));
		hessianEntries.push_back(Tripletd(1, 0, dyx));
		hessianEntries.push_back(Tripletd(1, 1, dyy));
	}

	double a, b;
};
