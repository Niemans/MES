#include "Node.h"

Node::Node() {
	T = nullptr;
	x = 0;
	y = 0;
	edge = 0;
}

Node::Node(double t0, double tend, double dt, double T0)
{
	int nr = (tend - t0) / dt + 1;
	T = new double[nr];
	T[0] = T0;
	for (int i = 1; i < nr; i++) {
		T[i] = 0;
	}
	x = 0;
	y = 0;
	edge = 0;
}


Node::~Node()
{
	if (T != nullptr) {
		delete[] T;
	}
}
