#include "Node.h"

struct Ntables //jest napisany dla 4 punkt�w, nie dla 9!!!
{
    int Gnr; //liczba punkt�w gaussowskich
    int Dnr; //liczba wymiar�w
    int edgenr;
    int pointnr;

    double** dNdE;
    double** dNdn;
    double** edgeN;
    double** areaN;

    Ntables(int, int);
    ~Ntables();

    /* dla gaussnr r�wnego 0
                 n
                 |
                (+1)
                 |
            p4   |  p3
    --(-1)--------------(+1)----E
            p1   |  p2
                 |
                (-1)
    */

    void numbers(); //oblicza dN/du z N=(1+-E)(1+-n)/4 dla ka�dego z punkt�w, gdzie u to E lub n
    void edgePoints(); //oblicza N na �ciankach
    void areaPoints(); //obliczba N w pubktach na elemencie

    void printdNdE(); //print
    void printdNdn(); //print
    void printEdgeN(); //print
    void printAreaN(); //print

};