#include "gauss.h"

struct Node
{
    double x, y;
    bool edge;
    double* T;

    Node();
    Node(double t0, double tend, double dt, double T0);

    ~Node();
};
