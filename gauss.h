#include <iostream>
#include <math.h>
#include <iomanip>
#include <fstream>
#include <string>

using namespace std;

struct Gauss
{
    static double** xN; //miejsca x w tablece Gaussa
    static double** wN; //d³ugoœci odcinka w

    static void data(); // przypisanie danych do zaalokowanych tablic
    static double* systemOfEquationsGauss(double** matrix, double* vector, int size); //rozwi¹zywanie uk³adów równañ za pomoc¹ eliminacji gaussa
    static void deldata();  //zwolnienie pamiêci
};