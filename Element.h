#include "Ntables.h"

struct Element //�cianki
{
    int number;
    int Gnr;
    int Dnr;
    Ntables N;

    int* ID_tab;        //zazwyczaj 4, bo 4 wierzcho�ki
    double* det;        //det jakobian�w
    double** J;         //tablica Jakobian�w, gdzie ka�dy Jakobian jest macierz�
    double** J_inv;     //tablica Jakobian�w^(-1), gdzie ka�dy Jakobian^(-1) jest macierz�
    double** H;         // macierz H elementu
    double*** H_edge;   //macierze H dla �cian zewn�trznych
    double** Hbc;       //suma Hbc
    double* Pbc;        //macierz 4x1 -> wektor P elementu
    double** C;         //macierz c elementu

    Element(int,int, int);
    ~Element();

    void det2D(int jnr); //jnr - numer jakobianu 
    void Jakobian2D(int jnr, Node** node); //jnr ^, Ntables zbi�r N/dE i N/dn, tablice x i y punkt�w
    void Hmatrix(double k_wsp); //wsp�czynnik k(T)
    void edgeHmatrix(double edgeLength[], double alpha); //d�ugo�� kraw�dzi 1 elementu, wsp. alfa //oblicza H_edge i Hbc
    void edgePmatrix(double edgeLength[], double alpha, double tot); //d�ugo�ci kraw�dzi 1 elementu, wsp. alfa i temp. Tot;
    void Cmatrix(double c_wsp, double ro); //sta�e pobierane na pocz�tku

    void printJakobian2x2(int jnr, bool inv); //print konkretnego jakobiana
    void printH(); //print macierz H danej elementu
    void printH_edge(); //print macierz H_edge danego elementu
    void printHbc(); //print Hbc
    void printPbc(); //print Pbc
    void printC(); //print macierzy C
};