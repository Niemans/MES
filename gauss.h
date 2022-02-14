#include <iostream>
#include <math.h>
#include <iomanip>
#include <fstream>
#include <string>

using namespace std;

struct Gauss
{
    static double** xN; //miejsca x w tablece Gaussa
    static double** wN; //d�ugo�ci odcinka w

    static void data(); // przypisanie danych do zaalokowanych tablic
    static double* systemOfEquationsGauss(double** matrix, double* vector, int size); //rozwi�zywanie uk�ad�w r�wna� za pomoc� eliminacji gaussa
    static void deldata();  //zwolnienie pami�ci
};