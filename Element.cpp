#include "Element.h"

Element::Element(int enr, int gnr, int dnr): N(gnr,dnr)
{
	number = enr;
	Gnr = gnr;
	Dnr = dnr;

	ID_tab = new int[N.pointnr]; //element ma 4 wierzcho³ki
	det = new double[Gnr * Gnr];

	J = new double* [Gnr * Gnr];
	J_inv = new double* [Gnr * Gnr];
	for (int i = 0; i < Gnr * Gnr; i++)
	{
		J[i] = new double[Dnr * Dnr];
		J_inv[i] = new double[Dnr * Dnr];
	}

	H = new double* [N.pointnr];
	C = new double* [N.pointnr];
	for (int i = 0; i < N.pointnr; i++)
	{
		H[i] = new double[N.pointnr];
		C[i] = new double[N.pointnr];
	}
	
	H_edge = new double** [N.edgenr];
	for (int i = 0; i < N.edgenr; i++)
	{
		H_edge[i] = new double* [N.pointnr];
		for (int j = 0; j < N.pointnr; j++)
		{
			H_edge[i][j] = new double [N.pointnr];
		}
	}

	Hbc = new double* [N.pointnr];
	for (int i = 0; i < N.pointnr; i++)
	{
		Hbc[i] = new double[N.pointnr];
	}
	Pbc = new double[N.pointnr];
}

Element::~Element()
{
	for (int i = 0; i < Gnr * Gnr; i++)
	{
		delete[] J[i];
		delete[] J_inv[i];
	}
	for (int i = 0; i < Dnr * Dnr; i++)
	{
		delete[] H[i];
		delete[] C[i];
	}
	for (int i = 0; i < N.edgenr; i++)
	{
		for (int j = 0; j < N.pointnr; j++)
		{
			delete[] H_edge[i][j];
		}
		delete[] H_edge[i];
	}
	for (int i = 0; i < N.pointnr; i++)
	{
		delete[] Hbc[i];
	}
	delete[] H_edge;
	delete[] Hbc;
	delete[] J;
	delete[] J_inv;
	delete[] H;
	delete[] ID_tab;
	delete[] det;
	delete[] Pbc;
	delete[] C;
}

void Element::det2D(int gnr) //gnr - numer jakobianu 
{
	det[gnr] = J[gnr][0] * J[gnr][3] - J[gnr][1] * J[gnr][2];
}

void Element::Jakobian2D(int jnr, Node** node)
{

	for (int i = 0; i < Dnr * Dnr; i++)
	{
		J[jnr][i] = 0;
		J_inv[jnr][i] = 0;
	}


	for (int i = 0; i < Dnr * Dnr; i++)
	{
		J[jnr][0] += N.dNdE[jnr][i] * node[i]->x;
		J[jnr][1] += N.dNdE[jnr][i] * node[i]->y;

		J[jnr][2] += N.dNdn[jnr][i] * node[i]->x;
		J[jnr][3] += N.dNdn[jnr][i] * node[i]->y;
	}

	det2D(jnr);

	//jakobian odwrotny
	J_inv[jnr][0] = J[jnr][3] * 1 / det[jnr];
	J_inv[jnr][1] = -J[jnr][1] * 1 / det[jnr];
	J_inv[jnr][2] = -J[jnr][2] * 1 / det[jnr];
	J_inv[jnr][3] = J[jnr][0] * 1 / det[jnr];
}

void Element::Hmatrix(double k_wsp)
{
	double** dNdx = new double* [Gnr * Gnr];
	double** dNdy = new double* [Gnr * Gnr];
	double*** sideH = new double** [Gnr * Gnr];
	for (int i = 0; i < Gnr * Gnr; i++)
	{
		dNdx[i] = new double[Dnr * Dnr];
		dNdy[i] = new double[Dnr * Dnr];
		sideH[i] = new double*[Dnr * Dnr];
		for (int j = 0; j < Dnr * Dnr; j++)
		{
			sideH[i][j] = new double[Dnr * Dnr];
		}
	}

	//macierze N/dx i N/dy
	for (int i = 0; i < Gnr * Gnr; i++)
	{
		for (int j = 0; j < Dnr * Dnr; j++)
		{
			dNdx[i][j] =  N.dNdE[i][j] * J_inv[i][0] + N.dNdn[i][j] * J_inv[i][1];
			dNdy[i][j] =  N.dNdE[i][j] * J_inv[i][2] + N.dNdn[i][j] * J_inv[i][3];
		}

	}

	
	/*cout << "N/dx:" << endl;
	for (int i = 0; i < Gnr * Gnr; i++)
	{
		for (int j = 0; j < Dnr * Dnr; j++)
		{
			cout << dNdx[i][j] << " ";
		}
		cout << endl;
	}
	cout << endl;
	cout << "N/dy:" << endl;
	for (int j = 0; j < Gnr * Gnr; j++)
	{
		for (int i = 0; i < Dnr * Dnr; i++)
		{
			cout << dNdy[j][i] << " ";
		}
		cout << endl;
	}
	cout << endl;*/

	//H poboczne
	for (int k = 0; k < Gnr * Gnr; k++)
	{
		for (int i = 0; i < Dnr * Dnr; i++)
		{
			for (int j = 0; j < Dnr * Dnr; j++)
			{
				sideH[k][i][j] = k_wsp * (dNdx[k][i] * dNdx[k][j] + dNdy[k][i] * dNdy[k][j]) * det[j];
			}
		}
	}
	
	
	//for (int i = 0; i < Gnr * Gnr; i++)
	//{
	//	cout << "H nr: " << i << endl;
	//	for (int j = 0; j < Dnr * Dnr; j++)
	//	{
	//		for (int z = 0; z < Dnr * Dnr; z++)
	//		{
	//			cout << sideH[i][j][z] << " ";
	//		}
	//		cout << endl;
	//	}
	//	cout << endl;
	//}
	//cout << endl;

	//czyszczenie pamiêci
	for (int i = 0; i < Dnr * Dnr; i++)
	{
		for (int j = 0; j < Dnr * Dnr; j++)
		{
			H[i][j] = 0;
		}
	}

	//H g³ówne
	for (int i = 0; i < Dnr * Dnr; i++)
	{
		for (int j = 0; j < Dnr * Dnr; j++)
		{
			for (int k = 0; k < Gnr * Gnr; k++)
			{
				H[i][j] += sideH[k][i][j] * Gauss::wN[Gnr - 2][k / Gnr] * Gauss::wN[Gnr - 2][k % Gnr];
			}
		}
	}


	//usuwanie pamiêci
	for (int i = 0; i < Gnr * Gnr; i++)
	{
		for (int j = 0; j < Dnr * Dnr; j++)
		{
			delete[] sideH[i][j];
		}

		delete[] dNdx[i];
		delete[] dNdy[i];
		delete[] sideH[i];
	}
	delete[] dNdx;
	delete[] dNdy;
	delete[] sideH;

}

void Element::edgeHmatrix(double edgeLength[], double alpha)
{
	//new
	double* edge_length = new double[N.edgenr];
	for (int i = 0; i < N.edgenr; i++)
	{
		edge_length[i] = edgeLength[i] / 2.;
	}

	double*** sideH = new double**[Gnr * N.edgenr];
	for (int i = 0; i < Gnr * N.edgenr; i++)
	{
		sideH[i] = new double*[N.pointnr];
		for (int j = 0; j < N.pointnr; j++)
		{
			sideH[i][j] = new double[N.pointnr];
		}
	}

	//side H i zerowanie H_edge
	for (int i = 0; i < Gnr * N.edgenr; i++)
	{
		for (int j = 0; j < N.pointnr; j++)
		{
			for (int k = 0; k < N.pointnr; k++)
			{
				sideH[i][j][k] = alpha * N.edgeN[i][j] * N.edgeN[i][k] * Gauss::wN[Gnr - 2][i % Gnr];
				if (i < N.edgenr) H_edge[i][j][k] = 0;
			}
		}
	}

	//for (int i = 0; i < Gnr * N.edgenr; i++)
	//{
	//	cout << "sideH nr: " << i << endl;
	//	for (int j = 0; j < N.pointnr; j++)
	//	{
	//		for (int k = 0; k < N.pointnr; k++)
	//		{
	//			cout << sideH[i][j][k] << " ";
	//		}
	//		cout << endl;
	//	}
	//	cout << endl;
	//}
	//cout << endl;

	//H_edge = wszystkie side H jego œciany (jest ich Gnr) * edge_length
	for (int i = 0; i < N.edgenr * Gnr; i++)
	{
		for (int j = 0; j < N.pointnr; j++)
		{
			for (int k = 0; k < N.pointnr; k++)
			{
				H_edge[i/Gnr][j][k] += (sideH[i][j][k] * edge_length[i/Gnr]);
			}
		}
	}

	//Hbc = suma H_edge
	for (int i = 0; i < N.pointnr; i++)
	{
		for (int j = 0; j < N.pointnr; j++)
		{
			Hbc[i][j] = 0;
			for (int k = 0; k < N.edgenr; k++)
			{
				Hbc[i][j] += H_edge[k][i][j];
			}
		}
	}

	//delete
	for (int i = 0; i < Gnr * N.edgenr; i++)
	{
		for (int j = 0; j < N.pointnr; j++)
		{
			delete[] sideH[i][j];
		}
		delete sideH[i];
	}
	delete[] sideH;
	delete[] edge_length;
}

void Element::edgePmatrix(double edgeLength[], double alpha, double tot)
{
	double* edge_length = new double[N.edgenr];
	for (int i = 0; i < N.edgenr; i++)
	{
		edge_length[i] = edgeLength[i] / 2.;
	}

	double** sideP = new double* [Gnr * N.edgenr];
	for (int i = 0; i < Gnr * N.edgenr; i++)
	{
		sideP[i] = new double[N.pointnr];
	}


	//cout << "Side P:" << endl;
	for (int i = 0; i < Gnr * N.edgenr; i++)
	{
		for (int j = 0; j < N.pointnr; j++)
		{
			sideP[i][j] = alpha * N.edgeN[i][j] * tot * Gauss::wN[Gnr - 2][i % Gnr] * edge_length[i/Gnr];
			//cout << sideP[i][j] << " ";
		}
		//cout << endl;
	}
	//cout << endl;

	for (int i = 0; i < N.pointnr; i++)
	{
		Pbc[i] = 0;
	}

	for (int i = 0; i < N.edgenr * Gnr; i++)
	{
		for (int k = 0; k < N.pointnr; k++)
		{
			Pbc[k] += sideP[i][k];
		}
	}

	for (int i = 0; i < Gnr * N.edgenr; i++)
		delete[] sideP[i];
	delete[] sideP;
	delete[] edge_length;
}

void Element::Cmatrix(double c_wsp, double ro)
{
	double*** helpC = new double** [Gnr * Gnr];
	for (int i = 0; i < Gnr * Gnr; i++)
	{
		helpC[i] = new double* [N.pointnr];
		for (int j = 0; j < N.pointnr; j++)
		{
			helpC[i][j] = new double[N.pointnr];
		}
	}

	//cout << "HelpC:" << endl;
	for (int k = 0; k < Gnr * Gnr; k++)
	{
		for (int i = 0; i < N.pointnr; i++)
		{
			for (int j = 0; j < N.pointnr; j++)
			{
				helpC[k][i][j] = c_wsp * ro * (N.areaN[k][i] * N.areaN[k][j]) * det[j];
				//cout << setw(8) << helpC[k][i][j] << " ";
			}
			//cout  << endl;
		}
		//cout << endl;
	}
	//cout << endl;

	for (int i = 0; i < N.pointnr; i++)
	{
		for (int j = 0; j < N.pointnr; j++)
		{
			C[i][j] = 0;
		}
	}
		
	//C g³ówne
	for (int i = 0; i < N.pointnr; i++)
	{
		for (int j = 0; j < N.pointnr; j++)
		{
			for (int k = 0; k < Gnr * Gnr; k++)
			{
				C[i][j] += helpC[k][i][j] * Gauss::wN[Gnr - 2][k / Gnr] * Gauss::wN[Gnr - 2][k % Gnr];
			}
		}
	}

	for (int i = 0; i < Gnr * Gnr; i++)
	{
		for (int j = 0; j < N.pointnr; j++)
		{
			delete[] helpC[i][j];
		}

		delete[] helpC[i];
	}
	delete[] helpC;

}


void Element::printJakobian2x2(int jnr, bool inv)
{
	int pom = 0;
	int old_pom = 0;

	if (inv == 0)
	{
		cout << "J:\n";
		cout << J[jnr][0] << "  " << J[jnr][1] << endl;
		cout << J[jnr][2] << "  " << J[jnr][3] << endl;
	}
	else
	{
		cout << "J inv:\n";
		cout << J[jnr][0] << "  " << J[jnr][1] << endl;
		cout << J[jnr][2] << "  " << J[jnr][3] << endl;
	}
	cout << endl;
}

void Element::printH()
{
	cout << "H: " << endl;

	for (int i = 0; i < N.pointnr; i++)
	{
		for (int j = 0; j < N.pointnr; j++)
		{
			cout << H[i][j] << " ";
		}
		cout << endl;
	}
	cout << endl;
}

void Element::printH_edge()
{
	for (int i = 0; i < N.edgenr; i++)
	{
		cout << "H_edge nr: " << i << endl;
		for (int j = 0; j < N.pointnr; j++)
		{
			for (int k = 0; k < N.pointnr; k++)
			{
				cout << H_edge[i][j][k] << " ";
			}
			cout << endl;
		}
		cout << endl;
	}
	cout << endl;
}

void Element::printHbc()
{
	cout << "Hbc: " << endl;
	for (int i = 0; i < N.pointnr; i++)
	{
		for (int j = 0; j < N.pointnr; j++)
		{
			cout << Hbc[i][j] << " ";
		}
		cout << endl;
	}
	cout << endl;
}

void Element::printPbc()
{
	cout << "Pbc: " << endl;
	for (int i = 0; i < N.pointnr; i++)
	{
		cout << Pbc[i] << " ";
	}
	cout << endl;
}

void Element::printC()
{
	cout << "C: " << endl;

	for (int i = 0; i < N.pointnr; i++)
	{
		for (int j = 0; j < N.pointnr; j++)
		{
			cout << C[i][j] << " ";
		}
		cout << endl;
	}
	cout << endl;
}
