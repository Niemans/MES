 #include "grid.h"

Grid::Grid(double b = 0.1, double h = 0.2, int nb = 4, int nh = 5, double kwsp = 30., double alfa = 25., double specHeat =700., double dens = 7800., double tstart = 0, double tend = 500, double tau = 50, double T0 = 100, double ToT = 1200., int gnr = 2, int dnr = 2)
{
    //wpisanie wartoœci
    B = b;
    H = h;
    nB = nb;
    nH = nh;
    nN = nh * nb;
    nE = (nh - 1) * (nb - 1);
    k_wsp = kwsp;
    alpha = alfa;
    c = specHeat;
    ro = dens;
    this->T0 = T0;
    this->Tot = ToT;

    t_start = tstart;
    t_end = tend;
    dt = tau;
    nR = (t_end - t_start) / dt + 1;

    node_tab = new Node*[nN];
    for(int i = 0; i < nN; i++) {
        node_tab[i] = new Node(tstart, tend, tau, T0);
    }
    elem_tab = new Element*[nE];
    for (int i = 0; i < nE; i++)
    {
        elem_tab[i] = new Element(i, gnr, dnr);
    }

    H_global = new double* [nN];
    C_global = new double* [nN];
    for (int i = 0; i < nN; i++)
    {
        H_global[i] = new double[nN];
        C_global[i] = new double[nN];
    }
    P_global = new double[nN];

    //wpisywanie wartoœci dla danych punktów
    int new_line = 0;
    for (int i = 0; i < nN; i++)
    {
        node_tab[i]->edge = (i < nH || i > nN - nH || i % nH == 0 || i % nH == nH - 1);

        new_line += 1 * (i % nH == 0 && i != 0);
        node_tab[i]->x = new_line * B / (nB - 1.);
        node_tab[i]->y = (i % nH) * H / (nH - 1.);
    }

    //wpisywanie ID dla danych elementów
    new_line = 0;
    for (int i = 0; i < nE; i++)
    {
        new_line += 1 * (i % (nH - 1) == 0 && i != 0);
        elem_tab[i]->ID_tab[0] = i + new_line;
        elem_tab[i]->ID_tab[1] = i + nH + new_line;
        elem_tab[i]->ID_tab[2] = i + nH + new_line + 1;
        elem_tab[i]->ID_tab[3] = i + new_line + 1;
    }
}

Grid::Grid(string pathname, int gnr = 2, int dnr = 2, double tstart = 0) 
{
    getParameters("ads", gnr, dnr);
    t_start = tstart;
    nR = (t_end - t_start) / dt + 1;

    H_global = new double* [nN];
    C_global = new double* [nN];
    for (int i = 0; i < nN; i++)
    {
        H_global[i] = new double[nN];
        C_global[i] = new double[nN];
    }
    P_global = new double[nN];

}

Grid::~Grid()
{
    for (int i = 0; i < nE; i++)
    {
        delete elem_tab[i];
    }
    delete[] elem_tab;

    for (int i = 0; i < nN; i++)
    {
        delete[] H_global[i];
        delete[] C_global[i];
        delete node_tab[i];
    }
    delete[] H_global;
    delete[] P_global;
    delete[] C_global;
    delete[] node_tab;
}

void Grid::getParameters(string name, int gnr, int dnr) {
    //pocz¹tkowe zmienne
    string myText;
    int bc_nr = 0;
    double info[10] = { 0 };

    //pierwsze otwarcie pliku - sprawdzanie liczby danego fragmentu
    ifstream MyReadFile(name);
    while (getline(MyReadFile, myText) && myText[0] >= 'A' && myText[0] <= 'z') {}
    while (getline(MyReadFile, myText) && myText[1] != 'B') {}
    while (getline(MyReadFile, myText))
        for (int i = 0; i <= myText.size(); i++)
            if (myText[i] == ',' || myText[i] == '\0')
                bc_nr++;
    MyReadFile.close();


    //drugi odczyt pliku - zapis do alokowanych tablic
    ifstream MyReadFile2(name);

    int nr = 0;
    while (getline(MyReadFile2, myText) && myText[0] >= 'A' && myText[0] <= 'z') {
        for (unsigned int i = 0; i < myText.size(); i++)
            if (myText[i] >= '0' && myText[i] <= '9')
                info[nr] = info[nr] * 10 + myText[i] - 48;
        nr++;
    }

    nN = info[8];
    nE = info[9];
    double** node = new double* [nN];
    for (int i = 0; i < nN; i++) {
        node[i] = new double[2]{ 0,0 };
    }
    int** elem = new int* [nE];
    for (int i = 0; i < nE; i++) {
        elem[i] = new int[4]{ 0,0,0,0 };
    }
    int* bc = new int[bc_nr];
    for (int i = 0; i < bc_nr; i++) {
        bc[i] = 0;
    }

    nr = 0;
    string number = "";
    int coma = 0;

    while (getline(MyReadFile2, myText) && myText[1] != 'E') {
        for (int i = 0; i < myText.size(); i++) {
            if ((coma == 1 || coma == 2) && (myText[i] >= '0' && myText[i] <= '9' || myText[i] <= '.')) {
                number += myText[i];
            }

            if (myText[i] == ',') {
                if (coma == 1) {
                    node[nr][0] = stod(number);
                    number = "";
                }
                coma++;
            }
            if (i == myText.size() - 1) {
                node[nr][1] = stod(number);
                number = "";
            }
        }
        nr++;
        coma = 0;
    }

    nr = 0;
    number = "";
    coma = 0;

    while (getline(MyReadFile2, myText) && myText[1] != 'B') {
        for (int i = 0; i < myText.size(); i++) {
            if ((coma != 0) && (myText[i] >= '0' && myText[i] <= '9')) {
                number += myText[i];
            }

            if (myText[i] == ',') {
                if (coma != 0) {
                    elem[nr][coma - 1] = stoi(number);
                    number = "";
                }
                coma++;
            }
            if (i == myText.size() - 1) {
                elem[nr][3] = stoi(number);
                number = "";
            }
        }
        nr++;
        coma = 0;
    }

    number = "";
    coma = 0;

    while (getline(MyReadFile, myText)) {
        for (int i = 0; i < myText.size(); i++) {
            if (myText[i] >= '0' && myText[i] <= '9') {
                number += myText[i];
            }

            if (myText[i] == ',') {
                bc[coma] = stoi(number) - 1;
                number = "";
                coma++;
            }
            if (i == myText.size() - 1) {
                bc[coma] = stoi(number) - 1;
                number = "";
            }
        }
        coma = 0;
    }

    MyReadFile.close();


    //zapis do tablicy Node'ów i elementów 
    this->t_end = info[0];
    this->dt = info[1];
    this->k_wsp = info[2];
    this->alpha = info[3];
    this->Tot = info[4];
    this->T0 = info[5];
    this->ro = info[6];
    this->c = info[7];
    this->nN = info[8];
    this->nE = info[9];

    node_tab = new Node * [nN];
    for (int i = 0; i < nN; i++) {
        node_tab[i] = new Node(this->t_start, this->t_end, this->dt, T0);
        node_tab[i]->x = node[i][0];
        node_tab[i]->y = node[i][1];
    }
    for (int i = 0; i < bc_nr; i++) {
        node_tab[bc[i]]->edge = 1;
    }

    elem_tab = new Element * [nE];
    for (int i = 0; i < nE; i++)
    {
        elem_tab[i] = new Element(i, gnr, dnr);
        for (int j = 0; j < 4; j++) { //4 wierzcho³ki
            elem_tab[i]->ID_tab[j] = elem[i][j];
        }
    }
}


void Grid::globalHmatrix()
{
    for (int i = 0; i < nN; i++)
    {
        for (int j = 0; j < nN; j++)
        {
            H_global[i][j] = 0;
        }
    }

    for (int i = 0; i < nE; i++)
    {
        for (int j = 0; j < 4; j++) //4 to N.pointnr, bo liczba ID 
        {
            for (int k = 0; k < 4; k++) //4 N.pointnr,  bo liczba ID 
            {
                H_global[elem_tab[i]->ID_tab[j]][elem_tab[i]->ID_tab[k]] += elem_tab[i]->H[j][k] + elem_tab[i]->Hbc[j][k];
            }
        }
    }
}

void Grid::globalPmatrix()
{
    for (int j = 0; j < nN; j++)
    {
        P_global[j] = 0;
    }

    for (int i = 0; i < nE; i++)
    {
        for (int j = 0; j < 4; j++) //4 to N.pointnr, bo liczba ID
        {
            P_global[elem_tab[i]->ID_tab[j]] += elem_tab[i]->Pbc[j];
        }
    }
}

void Grid::globalCmatrix()
{
    for (int i = 0; i < nN; i++)
    {
        for (int j = 0; j < nN; j++)
        {
            C_global[i][j] = 0;
        }
    }

    for (int i = 0; i < nE; i++)
    {
        for (int j = 0; j < 4; j++) //4 to N.pointnr, bo liczba ID 
        {
            for (int k = 0; k < 4; k++) //4 N.pointnr,  bo liczba ID 
            {
                C_global[elem_tab[i]->ID_tab[j]][elem_tab[i]->ID_tab[k]] += elem_tab[i]->C[j][k];
            }
        }
    }
}

void Grid::results()
{
    for (int i = 0; i < nN; i++) //H = H + C/dt
    {
        for (int j = 0; j < nN; j++)
        {
            H_global[i][j] += (C_global[i][j] / dt);
        }
    }

    double* newP = new double[nN];
    for (int k = 1; k < nR; k++) {
        //cout << "pom nr:" << k << endl;
        for (int i = 0; i < nN; i++) //P = P + (C/dt) * T0
        {
            newP[i] = P_global[i];
            for (int j = 0; j < nN; j++) {
                newP[i] += (C_global[i][j] / dt) * node_tab[j]->T[k - 1];
            }
            //cout << newP[i] << endl;
        }
       // cout << endl;

        double* T1 = Gauss::systemOfEquationsGauss(H_global, newP, nN);
        for (int i = 0; i < nN; i++) {
            node_tab[i]->T[k] = T1[i];
        }
        delete[] T1;
    }

    delete[] newP;
}


int* Grid::printElementID(int nr) //wyœwietla konkretny element z jego punktami
{
    if (nr >= nE || nr < 0)
    {
        cout << "No element has this number";
        return nullptr;
    }
    else
    {
        cout << "nodes' IDs around element no. " << nr << endl;
        cout << setw(static_cast<streamsize>(trunc(log10(nN))) + 1);
        cout << elem_tab[nr]->ID_tab[3] << " " << elem_tab[nr]->ID_tab[2];
        cout << endl << setw(static_cast<streamsize>(trunc(log10(nN))) + 1);
        cout << elem_tab[nr]->ID_tab[0] << " " << elem_tab[nr]->ID_tab[1];
        cout << endl;
    }
    return nullptr;
}

void Grid::printNodesPosition() //wyœwietla konkretne wspó³rzêdne wszystkich punktów
{
    cout << "Position of nodes (x,y)" << endl;
    for (int i = nH - 1; i >= 0; i--)
    {
        for (int j = 0; j < nB; j++)
        {
            cout << "(" << node_tab[j * nH + i]->x << "," << node_tab[j * nH + i]->y << ") ";
        }
        cout << endl << endl;
    }
}

void Grid::printNodesID() //wyœwietla ca³oœæ
{
    cout << "nodes' IDs" << endl;
    for (int i = nH - 1; i >= 0; i--)
    {
        for (int j = 0; j < nB; j++)
        {
            cout << setw(static_cast<streamsize>(trunc(log10(nN))) + 1);
            cout << j * nH + i << " ";
        }
        cout << endl << endl;
    }
}

void Grid::printNodeEdging()
{
    cout << "1 is edge, 0 are insides:" << endl;
    for (int i = nH - 1; i >= 0; i--)
    {
        for (int j = 0; j < nB; j++)
        {
            cout <<  node_tab[j * nH + i]->edge << "  ";
        }
        cout << endl << endl;
    }
}

void Grid::printHglobal()
{
    cout << "H globalne: " << endl;
    for (int i = 0; i < nN; i++)
    {
        for (int j = 0; j < nN; j++)
        {
            cout << setw(7) << H_global[i][j];
        }
        cout << endl;
    }
    cout << endl;
}

void Grid::printPglobal()
{
    cout << "P globalne: " << endl;

    for (int i = 0; i < nN; i++)
    {
        cout << P_global[i] << endl;
    }
    cout << endl;
}

void Grid::printCglobal()
{
    cout << "C globalne: " << endl;
    for (int i = 0; i < nN; i++)
    {
        for (int j = 0; j < nN; j++)
        {
            cout << setw(10) << C_global[i][j];
        }
        cout << endl;
    }
    cout << endl;
}

void Grid::printResults()
{
    for (int i = 0; i < nR*dt; i++) {
        cout << "time: " << i << endl;
        for (int j = 0; j < nN; j++) {
            cout << setw(3) << j << ". " << setw(8) << node_tab[j]->T[i] << "   ";
            if (j%(15) == 0) {
                cout << endl;
            }

        }
        cout << endl;
    }
    cout << endl;
}

void Grid::printMinMax()
{
    double min;
    double max;
    cout << "Step  &  Minimum & Maximum: "<< endl;
    for (int i = 0; i < nR; i++) {
        min = node_tab[0]->T[i];
        max = node_tab[0]->T[i];
        for (int j = 0; j < nN; j++) {
            if (min > node_tab[j]->T[i]) min = node_tab[j]->T[i];
            if (max < node_tab[j]->T[i]) max = node_tab[j]->T[i];
        }
        cout << setw(5) << setprecision(1) << i*dt << "  " << setprecision(15)  << min << "  " << max << endl;
    }
    cout << endl;
}
