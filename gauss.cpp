#include "gauss.h"

void Gauss::data() // przypisanie danych do zaalokowanych tablic
{
    xN = new double* [2];
    xN[0] = new double[2];
    xN[1] = new double[3];

    wN = new double* [2];
    wN[0] = new double[2];
    wN[1] = new double[3];

    xN[0][0] = -1. / sqrt(3.); xN[0][1] = 1. / sqrt(3.);
    xN[1][0] = -sqrt(3. / 5.); xN[1][1] = 0; xN[1][2] = sqrt(3. / 5.);

    wN[0][0] = 1.; wN[0][1] = 1.;
    wN[1][0] = 5. / 9.; wN[1][1] = 8. / 9.; wN[1][2] = 5. / 9.;
}

double* Gauss::systemOfEquationsGauss(double** matrix, double* vector, int size) //HxT + P = 0
{
	if (size < 1)
	{
		cout << "matrixes does not exist, its size is less then 1" << '\a' << endl;
		return nullptr;
	}

	for (int i = 0; i < size; i++)
	{
		if (matrix[i][i] == 0)
		{
			cout << "values in diagonal are 0" << '\a' << endl;
			return nullptr;
		}
	}

	double* temperature_tab = new double[size];
	double** help_tab = new double* [size];
	for (int i = 0; i < size; i++)
		help_tab[i] = new double[size + 1];

	for (int i = 0; i < size; i++)
		for (int j = 0; j < size + 1; j++)
		{
			help_tab[i][j] = matrix[i][j];
			if (j == size)
				help_tab[i][j] = vector[i];
		}

	for (int i = 0; i < size - 1; i++)
		for (int j = i + 1; j < size; j++)
		{
			double help = help_tab[j][i];
			for (int k = i; k < size + 1; k++)
				help_tab[j][k] = help_tab[j][k] - ((help / help_tab[i][i]) * help_tab[i][k]);
		}



	for (int i = size - 1; i >= 0; i--)
	{
		temperature_tab[i] = help_tab[i][size];
		for (int j = size - 1; j >= i; j--)
		{
			if (j != i) temperature_tab[i] -= help_tab[i][j] * temperature_tab[j];
			else temperature_tab[i] = temperature_tab[i] / help_tab[i][j];
		}
	}

	for (int i = 0; i < size; i++)
	{
		delete[] help_tab[i];
	}
	delete[] help_tab;

	return temperature_tab;
}

void Gauss::deldata()  //zwolnienie pamiêci
{
    delete[] xN[1];
    delete[] xN[0];
    delete[] xN;

    delete[] wN[1];
    delete[] wN[0];
    delete[] wN;
}

