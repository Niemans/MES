#include "Ntables.h"



//mes 3  tableka N/dE i N/dn

Ntables::Ntables(int gaussnr, int dimensionnr)
{
    Gnr = gaussnr;
    Dnr = dimensionnr;
    pointnr = (int)(pow(2, dimensionnr));
    edgenr = dimensionnr * (int)pow(2,dimensionnr-1);


    dNdE = new double* [Gnr * Gnr];
    dNdn = new double* [Gnr * Gnr];
    areaN = new double* [Gnr * Gnr];
    for (int i = 0; i < Gnr * Gnr; i++)
    {
        dNdE[i] = new double[pointnr];
        dNdn[i] = new double[pointnr];
        areaN[i] = new double [pointnr];

    }

    edgeN = new double* [Gnr * edgenr];
    for (int i = 0; i < Gnr * edgenr; i++)
    {
        edgeN[i] = new double[pointnr];
    }


    numbers();
    edgePoints();
    areaPoints();
}

Ntables::~Ntables()
{
    for (int i = 0; i < Gnr * Gnr; i++)
    {
        delete[] dNdE[i];
        delete[] dNdn[i];
        delete[] areaN[i];
    }
    for (int i = 0; i < Gnr * 4; i++)
    {
        delete[] edgeN[i];
    }
    delete[] dNdE;
    delete[] dNdn;
    delete[] edgeN;
    delete[] areaN;
}

//N=(1+-E)(1+-n)/4
void Ntables::numbers()
{
    for (int i = 0; i < Gnr * Gnr; i++) 
    {
        if (Gnr == 2)
        {
            //N/dE = (+-1)*(1+-n)/4
            dNdE[i][0] = (-1.) * (1. - Gauss::xN[0][i / 2]) / 4.;       //--
            dNdE[i][1] =         (1. - Gauss::xN[0][i / 2]) / 4.;       //+-
            dNdE[i][2] =         (1. + Gauss::xN[0][i / 2]) / 4.;       //++
            dNdE[i][3] = (-1.) * (1. + Gauss::xN[0][i / 2]) / 4.;       //-+

            //N/dn = (1+-E)*(+-1)/4
            dNdn[i][0] = (1. - Gauss::xN[0][i % 3 > 0]) * (-1.) / 4.;   //--
            dNdn[i][1] = (1. + Gauss::xN[0][i % 3 > 0]) * (-1.) / 4.;   //+-
            dNdn[i][2] = (1. + Gauss::xN[0][i % 3 > 0]) / 4.;           //++
            dNdn[i][3] = (1. - Gauss::xN[0][i % 3 > 0]) / 4.;           //-+
        }
        else if (Gnr >= 3)
        {
            //N/dE = (+-1)*(1+-n)/4
            dNdE[i][0] = (-1.) * (1. - Gauss::xN[Gnr - 2][i % Gnr]) / 4.;       //--
            dNdE[i][1] =         (1. - Gauss::xN[Gnr - 2][i % Gnr]) / 4.;       //+-
            dNdE[i][2] =         (1. + Gauss::xN[Gnr - 2][i % Gnr]) / 4.;       //++
            dNdE[i][3] = (-1.) * (1. + Gauss::xN[Gnr - 2][i % Gnr]) / 4.;       //-+

            dNdn[i][0] = (1. - Gauss::xN[Gnr - 2][i / Gnr]) * (-1.) / 4.;       //--
            dNdn[i][1] = (1. + Gauss::xN[Gnr - 2][i / Gnr]) * (-1.) / 4.;       //+-
            dNdn[i][2] = (1. + Gauss::xN[Gnr - 2][i / Gnr])         / 4.;       //++
            dNdn[i][3] = (1. - Gauss::xN[Gnr - 2][i / Gnr])         / 4.;       //-+
        }
        else
        {
            cout << "Aktualnie sa tylko przewidziane 2 i 3 punktowe metody calkowania Gaussa";
            exit(-1);
        }
    }
}

void Ntables::edgePoints()
{
    int help = 0;
    for (int i = 0; i < Gnr * edgenr; i++)
    {
        for (int j = 0; j < pointnr; j++)
        {
            edgeN[i][j] = 0;
        }
    }

    for (int i = 0; i < Gnr * edgenr; i++)
    {
        if (i != 0 && i % Gnr == 0)
        {
            help++;
        }

        for (int j = 0; j < pointnr; j++)
        {
            edgeN[i][(j + help) % pointnr] = (1. - Gauss::xN[Gnr - 2][i % Gnr]) / 2.;
            edgeN[i][(j + 1 + help) % pointnr] = 1 - edgeN[i][(j + help) % pointnr];

            for (int k = 2; k < pointnr; k++)
            {
                edgeN[i][(j + k + help) % (pointnr)] = 0;
            }
        }
    }  
}

void Ntables::areaPoints()
{
    //N=(1+-E)(1+-n)/4
    for (int i = 0; i < Gnr * Gnr; i++)
    {
        
        if (Gnr == 2)
        {
            areaN[i][0] = (1. - Gauss::xN[0][i % 3 > 0]) * (1. - Gauss::xN[0][i / 2]) / 4.; //--
            areaN[i][1] = (1. + Gauss::xN[0][i % 3 > 0]) * (1. - Gauss::xN[0][i / 2]) / 4.; //+-
            areaN[i][2] = (1. + Gauss::xN[0][i % 3 > 0]) * (1. + Gauss::xN[0][i / 2]) / 4.; //++
            areaN[i][3] = (1. - Gauss::xN[0][i % 3 > 0]) * (1. + Gauss::xN[0][i / 2]) / 4.; //+-

        }
        else if (Gnr >= 3)
        {
            areaN[i][0] = (1. - Gauss::xN[Gnr - 2][i / Gnr]) * (1. - Gauss::xN[Gnr - 2][i % Gnr]) / 4.; //--
            areaN[i][1] = (1. + Gauss::xN[Gnr - 2][i / Gnr]) * (1. - Gauss::xN[Gnr - 2][i % Gnr]) / 4.; //+-
            areaN[i][2] = (1. + Gauss::xN[Gnr - 2][i / Gnr]) * (1. + Gauss::xN[Gnr - 2][i % Gnr]) / 4.; //++
            areaN[i][3] = (1. - Gauss::xN[Gnr - 2][i / Gnr]) * (1. + Gauss::xN[Gnr - 2][i % Gnr]) / 4.; //-+
        }
        else
        {
            cout << "Aktualnie sa tylko przewidziane 2 i 3 punktowe metody calkowania Gaussa";
            exit(-1);
        }
    }
}



void Ntables::printdNdE()
{
   cout << "dN/dE:" << endl;
   for (int i = 0; i < Gnr * Gnr; i++)
   {
       cout << dNdE[i][0] << "   " << dNdE[i][1] << "   " << dNdE[i][2] << "   " << dNdE[i][3] << endl;
   }
   cout << endl;
}

void Ntables::printdNdn()
{
    cout << "dN/dn:" << endl;
    for (int i = 0; i < Gnr * Gnr; i++)
    {
        cout << dNdn[i][0] << "   " << dNdn[i][1] << "   " << dNdn[i][2] << "   " << dNdn[i][3] << endl;
    }
    cout << endl; 
}

void Ntables::printEdgeN()
{
    cout << "EdgeNs:" << endl;
    for (int i = 0; i < Gnr * edgenr; i++)
    {
        cout << edgeN[i][0] << " " << edgeN[i][1] << " " << edgeN[i][2] << " " << edgeN[i][3] << " ";
        cout << endl;
    }
    cout << endl;
}

void Ntables::printAreaN()
{
    cout << "AreaNs:" << endl;
    for (int i = 0; i < Gnr * Gnr; i++)
    {
        cout << areaN[i][0] << " " << areaN[i][1] << " " << areaN[i][2] << " " << areaN[i][3] << " ";
        cout << endl;
    }
}



