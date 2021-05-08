#pragma once

#include "ObjectiveFunction.h"

#include <string>
#include <iostream>
#include <cmath>
#include <cfloat>

class GradientDescentFunctionMinimizer{
public:
	GradientDescentFunctionMinimizer(int maxIterations=100, double solveResidual=1e-5, int maxLineSearchIterations=15)
		: maxIterations(maxIterations), solveResidual(solveResidual), maxLineSearchIterations(maxLineSearchIterations){
	}

	virtual ~GradientDescentFunctionMinimizer(){}

	int getLastIterations() { return lastIterations; }

	virtual bool minimize(ObjectiveFunction *function, VectorXd &x){

		//number of parameters...
		int N = (int) x.size();
		resize(xi, N);
		resize(dx, N);
		resize(gradient, N);

		xi = x;

		bool optimizationConverged = false;
		int i=0;
		for(; i < maxIterations; i++) {
			computeSearchDirection(function, xi, dx);

			if (dx.norm() < solveResidual){
				optimizationConverged = true;
				break;
			}

			doLineSearch(function, dx, xi);
		}

		lastIterations = i;

		//p now holds the parameter values at the start of the iteration...
		x = xi;

		//and done!
		return optimizationConverged;
	}

protected:
	// Since the gradient of a function gives the direction of steepest descent, all one needs to do is go in that direction...

	virtual void computeSearchDirection(ObjectiveFunction *function, const VectorXd &x, VectorXd& dx) {

		// Ex. 1.1
		function->addGradientTo(dx, x);
		dx = -dx; // Search direction is the opposite of the gradient value

	}

	virtual void doLineSearch(ObjectiveFunction *function, const VectorXd& dx, VectorXd& xi)
	{

		// Ex. 1.1
		double step_size = 1;
		double scaling_factor = 0.5;
		int i = 0;
		VectorXd upd_x;

		upd_x = xi + step_size * dx;
		while ((function->computeValue(upd_x) >= function->computeValue(xi)) && i < maxLineSearchIterations ) {
			step_size = scaling_factor * step_size;
			upd_x = xi + step_size * dx;
			i++;
		}

		xi = upd_x;

	//	std::cout << "f(xi): " << function->computeValue(xi) << "( " << xi.row(0) << ", " << xi.row(1) << " )" << std::endl;
	//	std::cout << "dx: " << dx.norm() << std::endl;
	}	

protected:
	double solveResidual = 1e-5;
	int maxIterations = 100;
	int maxLineSearchIterations = 15;

	VectorXd xi, dx, gradient;

	// some stats about the last time `minimize` was called
	int lastIterations = -1;
};
