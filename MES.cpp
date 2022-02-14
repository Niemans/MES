#include "grid.h"

const int gaussNumber = 2; //2,3
const int dimentions = 2; //dzia³a na 100% tylko 2D

double  k = 25, a = 300, c = 700, ro = 7800; //wsp. k (przewodzenie ciep³a), alfa (konwekcyjnej wymiany ciep³a), c (pojemnoœæ cieplna), ro (gêstoœæ)
double tau = 1, tend = 100; //tau = dt (czas), tend czas koñcowy
double T0 = 100, Tot = 1200; //temp startowa, Temp otoczenia

//const double b = 0.025, h = 0.025;
//const int nB = 2, nH = 2;
//const double b = 0.1, h = 0.2;
//const int nB = 4, nH = 5;
const double b = 0.1, h = 0.1;
const int nB = 31, nH = 31;

double** Gauss::wN = nullptr;
double** Gauss::xN = nullptr;

int main()
{
    cout << fixed << setprecision(3);
    Gauss::data();
    Grid grid(b, h, nB, nH, k, a, c, ro, 0, tend, tau, T0, Tot, gaussNumber, dimentions);
    
    grid.printNodesPosition();
    grid.printNodesID();
    grid.printNodeEdging();  
    //grid.elem_tab[0]->N.printdNdE();
    //grid.elem_tab[0]->N.printdNdn();
    //grid.elem_tab[0]->N.printEdgeN();
    //grid.elem_tab[0]->N.printAreaN();

    Node* nodetab[4]; // 4 bo liczba node'ów w lelemencie to 4

    for (int i = 0; i < grid.nE; i++) //i to numer elementu
    {
        for (int j = 0; j < dimentions * dimentions; j++) //punkty potrzebne do stworzenia jakobianów dla danego elementu
        {
            nodetab[j] = grid.node_tab[grid.elem_tab[i][0].ID_tab[j]];
        }
        //cout << "------------------------------------------------------\n";
        //cout << "\t\t\t\tElement nr " << i << endl << endl;
        for (int j = 0; j < gaussNumber * gaussNumber; j++) //j to punkty po gaussie
        {
            grid.elem_tab[i]->Jakobian2D(j, nodetab);
            //cout << "Punkt lokalny nr " << j << endl;
            //grid.elem_tab[i]->printJakobian2x2(j,0);
            //grid.elem_tab[i]->printJakobian2x2(j,1);
            //cout << "det i 1/det:\n";
            //cout << setprecision(8) << grid.elem_tab[i]->det[j] << "   " << 1 / grid.elem_tab[i]->det[j] << endl << endl;
        }
        grid.elem_tab[i]->Hmatrix(k);
        //grid.elem_tab[i]->printH();
        grid.elem_tab[i]->Cmatrix(c,ro);
        //grid.elem_tab[i]->printC();
        
        double length[4]; //4 œciany ma ka¿dy z elementów w 2D
        //cout << "dlugosci odcinka:" << endl;
        for (int j = 0; j < 4; j++)
        {
            if (grid.node_tab[grid.elem_tab[i]->ID_tab[(j + 3) % 4]]->edge == 0 ||
                grid.node_tab[grid.elem_tab[i]->ID_tab[ j      % 4]]->edge == 0)
            {
                length[j] = 0;
            }
            else
            {  //sqrt((x1-x2)^2 + (y1-y2)^2)
                length[j] = sqrt
                (   pow(grid.node_tab[grid.elem_tab[i]->ID_tab[(j + 3) % 4]]->x -
                        grid.node_tab[grid.elem_tab[i]->ID_tab[j % 4]]->x, 2) +
                    pow(grid.node_tab[grid.elem_tab[i]->ID_tab[(j + 3) % 4]]->y -
                        grid.node_tab[grid.elem_tab[i]->ID_tab[j % 4]]->y, 2));
            }
           // cout << length[j] << endl;
        }

        grid.elem_tab[i]->edgeHmatrix(length,grid.alpha);
        //grid.elem_tab[i]->printH_edge();
        //grid.elem_tab[i]->printHbc();
        grid.elem_tab[i]->edgePmatrix(length, grid.alpha, Tot);
        //grid.elem_tab[i]->printPbc();
    }
    cout << "------------------------------------------------------\n";
    cout << "\t\t\t\tGlobalne " << endl << endl;
    grid.globalHmatrix();
    //grid.printHglobal();
    grid.globalPmatrix();
    //grid.printPglobal();
    grid.globalCmatrix();
    //grid.printCglobal();

    grid.results();
    //grid.printResults();
    grid.printMinMax();

    Gauss::deldata();
}

