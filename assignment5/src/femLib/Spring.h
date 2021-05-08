#pragma once

#include "Element.h"

/**
			This class implements the interface for an elementary energy unit. As a function of deformed, undeformed,
			and other parameters, such as boundary conditions, each class that extends this one will define a potential energy.
			The deformed energy depends on a number of nodes.
*/
class Spring : public Element {

public:
	Spring(const std::array<int, 2> &nodeIndices, const VectorXd &X)
		: nodeIndices(nodeIndices) {
	}
	virtual ~Spring() {}

	// Returns the number of nodes this unit depends on
	virtual int getNumNodes() const {
		return 2;
	}
	// Returns the global index of node `i`
	virtual int getNodeIndex(int i) const {
		return nodeIndices[i];
	}

	// Returns the element's mass
	virtual double getMass() const {
		return 0;
	}

	// Returns the energy value given deformed `x` and undeformed `X` state
	virtual double getEnergy(const VectorXd& x, const VectorXd& X) {

		// Ex 1.2
		// Task: Given `x` and `X`, return the spring energy.

		// Some notes:
		// `x` and `X` contain the current and rest positions of all
		// nodes. You can extract the position of e.g. node 0 like this:
		// Vector2d x1 = getVertex(0, x);
		// or to get the rest position of node 0:
		// Vector X1 = getVertex(0, X);
		// The spring stiffness is stored in the variable `k`.

		double energy ,L,l;

		Vector2d xi = getNodePos(0, x);
		Vector2d xj = getNodePos(1, x);
		Vector2d Xi = getNodePos(0, X);
		Vector2d Xj = getNodePos(1, X);

		L = (Xi - Xj).norm();
		l = (xi - xj).norm();
		energy = 0.5*k*(pow((l - L), 2));
	//	energy = 0.5*k*pow((l / L) - 1, 2)*L;
		return energy;
	}

	// Adds the gradient to `grad` given deformed `x` and undeformed `X` state
	virtual void addEnergyGradientTo(const VectorXd& x, const VectorXd& X, VectorXd& grad) {

		// Ex 1.2
		// Task: Given `x` and `X`, add the gradient of the spring energy to `grad`.

		// Again, you can extract the position of e.g. node 0 like this:
		// Vector2d x1 = getVertex(0, x);
		// and the spring stiffness is stored in `k`.

		// Remember that `grad` is a vector of size 2*N, where N is the total
		// number of nodes in the system. Make sure you are writing to the
		// correct location in `grad`. To get the global index of node 0 of
		// this spring, use this function:
		// int globalIndex0 = getNodeIndex(0);
		// or for node 1
		// int globalIndex1 = getNodeIndex(1);
		Vector2d x0 = getNodePos(0, x);
		Vector2d x1 = getNodePos(1, x);
		Vector2d X0 = getNodePos(0, X);
		Vector2d X1 = getNodePos(1, X);

		double l = (x0 - x1).norm();
		double L = (X0 - X1).norm();

		int globalIndex0 = getNodeIndex(0);
		int globalIndex1 = getNodeIndex(1);

		Vector2d LL; LL[0] = LL[1] = L;
		Vector2d d = (k / l) * ((x0 - x1) - LL);
		grad[2 * globalIndex0] = d[0];
		grad[2 * globalIndex0 + 1] = d[1];
		d = (k / l) * ((x1 - x0) - LL);
		grad[2 * globalIndex1] = d[0];
		grad[2 * globalIndex1 + 1] = d[1];


	}

	// Adds the hessian entries to `hesEntries` given deformed `x` and undeformed `X` state
	virtual void addEnergyHessianTo(const VectorXd& x, const VectorXd& X, std::vector<Tripletd>& hesEntries) {

		// Ex 1.4
		// Task: Given `x` and `X`, add the hessian of the spring energy to `hesEntries`.
		Vector2d x_0 = getNodePos(0, x);
		Vector2d x_1 = getNodePos(1, x);
		Vector2d X0 = getNodePos(0, X);
		Vector2d X1 = getNodePos(1, X);
		double L = (X0 - X1).norm();
		double x1 = x_0[0];
		double x2 = x_1[0];
		double y1 = x_0[1];
		double y2 = x_1[1];
		int globIndex1 = getNodeIndex(0);
		int globIndex2 = getNodeIndex(1);
		std::vector<double> hessian;

		double l = (x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2);

		// x1 subs:
		// x1
		double val = (k*(pow((x1 - x2), 2) + (x1 - x2)*(L - x1 + x2) + pow(y1 - y2, 2)) / ((l)*(sqrt(l))));
		hesEntries.push_back(Tripletd(2 * globIndex1, 2 * globIndex1, val));
		// y1
		val = (-k * ((y2 - y1)*(L - x1 + x2)) / ((l)*(sqrt(l))));
		hesEntries.push_back(Tripletd(2 * globIndex1, 2 * globIndex1 + 1, val));
		// x2
		val = (-k * (pow((x1 - x2), 2) + (x1 - x2)*(L - x1 + x2) + pow(y1 - y2, 2)) / ((l)*(sqrt(l))));
		hesEntries.push_back(Tripletd(2 * globIndex1, 2 * globIndex2, val));
		// y2
		val = (-k * ((y1 - y2)*(L - x1 + x2)) / ((l)*(sqrt(l))));
		hesEntries.push_back(Tripletd(2 * globIndex1, 2 * globIndex2 + 1, val));

		// y1 subs:
		//x1
		val = (-k * ((x2 - x1)*(L - y1 + y2)) / ((l)*(sqrt(l))));
		hesEntries.push_back(Tripletd(2 * globIndex1 + 1, 2 * globIndex1, val));
		//y1
		val = (k*(pow((x1 - x2), 2) - (y2 - y1)*(L - y1 + y2) + pow(-y1 + y2, 2)) / ((l)*(sqrt(l))));
		hesEntries.push_back(Tripletd(2 * globIndex1 + 1, 2 * globIndex1 + 1, val));
		//x2
		val = (-k * ((x1 - x2)*(L - y1 + y2)) / ((l)*(sqrt(l))));
		hesEntries.push_back(Tripletd(2 * globIndex1 + 1, 2 * globIndex2, val));
		//y2
		val = (-k * (pow((x1 - x2), 2) + (y1 - y2)*(L - y1 + y2) + pow(y1 - y2, 2)) / ((l)*(sqrt(l))));
		hesEntries.push_back(Tripletd(2 * globIndex1 + 1, 2 * globIndex2 + 1, val));

		// x2 subs
		//x1
		val = (-k * (pow((x2 - x1), 2) + (x2 - x1)*(L - x2 + x1) + pow(y1 - y2, 2)) / ((l)*(sqrt(l))));
		hesEntries.push_back(Tripletd(2 * globIndex2, 2 * globIndex1, val));
		//y1
		val = (-k * ((y2 - y1)*(L + x1 - x2)) / ((l)*(sqrt(l))));
		hesEntries.push_back(Tripletd(2 * globIndex2, 2 * globIndex1 + 1, val));
		//x2
		val = (k * (pow((x1 - x2), 2) - (x1 - x2)*(L + x1 - x2) + pow(y1 - y2, 2)) / ((l)*(sqrt(l))));
		hesEntries.push_back(Tripletd(2 * globIndex2, 2 * globIndex2, val));
		//y2
		val = (-k * ((y1 - y2)*(L + x1 - x2)) / ((l)*(sqrt(l))));
		hesEntries.push_back(Tripletd(2 * globIndex2, 2 * globIndex2 + 1, val));

		//y2 subs
		//x1
		val = (-k * ((x2 - x1)*(L + y1 - y2)) / ((l)*(sqrt(l))));
		hesEntries.push_back(Tripletd(2 * globIndex2 + 1, 2 * globIndex1, val));
		//y1
		val = (-k * (pow((x1 - x2), 2) + (y2 - y1)*(L - y2 + y1) + pow(-y1 + y2, 2)) / ((l)*(sqrt(l))));
		hesEntries.push_back(Tripletd(2 * globIndex2 + 1, 2 * globIndex1 + 1, val));
		//x2
		val = (-k * ((x1 - x2)*(L + y1 - y2)) / ((l)*(sqrt(l))));
		hesEntries.push_back(Tripletd(2 * globIndex2 + 1, 2 * globIndex2, val));
		//y2
		val = (k * (pow((x1 - x2), 2) - (y1 - y2)*(L + y1 - y2) + pow(y1 - y2, 2)) / ((l)*(sqrt(l))));
		hesEntries.push_back(Tripletd(2 * globIndex2 + 1, 2 * globIndex2 + 1, val));


	}

protected:
	// the collection of nodes that define the triangle element
	std::array<int, 2> nodeIndices;
	// spring stiffness
	double k = 20.0;
};

