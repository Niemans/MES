#include "Element.h"

struct Grid
{
    double H, B;
    long int nH, nB, nN, nE;
    Node** node_tab;
    Element** elem_tab;
    double k_wsp, alpha, ro, c;
    double** H_global;
    double* P_global;
    double** C_global;
    double T0, Tot;

    int nR;
    double  t_start, t_end, dt;

    Grid(double, double, int, int, double, double, double, double, double, double, double, double, double, int, int);
    Grid(string, int, int, double);

    ~Grid();

private: void getParameters(string pathname, int gnr, int dnr);
public:
    void globalHmatrix();
    void globalPmatrix();
    void globalCmatrix();

    void results();

    int* printElementID(int nr); //wyœwietla konkretny element z jego punktami
    void printNodesPosition(); //wyœwietla konkretne wspó³rzêdne wszystkich punktów
    void printNodesID(); //wyœwietla ca³oœæ
    void printNodeEdging();
    void printHglobal();
    void printPglobal();
    void printCglobal();
    void printResults();
    void printMinMax();
};
