# include <cstdlib>
# include <iostream>
# include <fstream>
# include <iomanip>
# include <cmath>
# include <ctime>
#include <string>
#include<vector>
#include <ctime>

using namespace std;

void rk4vec(double* xOut, double* yOut, double* TraceJ, double t0, double x0, double y0, double dt, double T,
    double A, double W, int n_system, double B_system, double a_system);

void calcLyapunov(double* L1out, double* L2out, double x0, double y0, double t0,
    double dt, double dT, double Trun, double tTransient, double A, double W,
    int n_system, double B_system, double a_system);

void writeFile(string fileName, vector <double> AVec, vector <double> WVec,
    vector <double> L1Vec, vector <double> L2Vec);

//double dx_dt(double t, double x, double y);
//double dy_dt(double t, double x, double y);

int main()
{

    cout << "Enter the set number:";

    int iSet;
    cin >> iSet;

    cout << "\n\n--------------START RUNNING------------------------\n";

    string readFileName = "inputFile_Set" + to_string(iSet) + ".txt";
    string wtireFileName = "outputFile_Set" + to_string(iSet) + ".txt";



    vector <string> vecInput;

    int iLine = 0;
    std::string delimiter = "=";


    std::ifstream filein(readFileName);
    for (std::string line; std::getline(filein, line); )
    {
        vecInput.push_back(line.substr(1 + line.find(delimiter)));
    }


    // system:
    int n_system = stoi(vecInput[0]);
    double B_system = stod(vecInput[1]);
    double a_system = stod(vecInput[2]);

    double x0 = stod(vecInput[3]);
    double y0 = stod(vecInput[4]);

    double t0 = stod(vecInput[5]);
    double dt = stod(vecInput[6]);
    double dT = stod(vecInput[7]);
    double tTransient = stod(vecInput[8]);
    double Trun = stod(vecInput[9]);

    double Amin = stod(vecInput[10]);
    double Amax = stod(vecInput[11]);
    double Ares = stod(vecInput[12]);

    double Wmin = stod(vecInput[13]);
    double Wmax = stod(vecInput[14]);
    double Wres = stod(vecInput[15]);

    int nWriteInterval = stoi(vecInput[16]);

    time_t tstart, tend;
    tstart = time(0);

    int nA = int((Amax - Amin) / Ares) + 1;
    int nW = int((Wmax - Wmin) / Wres) + 1;

    double W = 1.09; double A = 1.254; // Chaotic


    //writing info of simulation in the result file:

    ofstream fileOut;
    fileOut.open(wtireFileName);
    fileOut << "Lyapunov Exponents of the forced system with Bessel-like damping.\n";
    fileOut << " System info: (n, B, a) = (" << n_system << ", " << B_system << ", " << a_system << ")\n";
    fileOut << "(x0, y0, t0)= (" << x0 << ", " << y0 << ", " << t0 << ")\n";
    fileOut << "Runtime for calculations=" << Trun << ";  Transinet time=" << tTransient << ";  dT=" << dT << ";    dt=" << dt << endl;
    fileOut << "(Amin, Ares, Amax)= (" << Amin << ", " << Ares << ", " << Amax << "); number of cases:" << nA << "\n";
    fileOut << "(Wmin, Wres, Wmax)= (" << Wmin << ", " << Wres << ", " << Wmax << "); number of cases:" << nW << "\n";
    fileOut << "*******************************************************************************\n";
    fileOut << "A,      W,      L1,     L2\n";
    fileOut << "----------------------------------------\n";
    fileOut.close();

    double L1, L2;

    vector <double> AVec;
    vector <double> WVec;
    vector <double> L1Vec;
    vector <double> L2Vec;


    for (int iA = 0; iA < nA; iA++)
    {
        for (int iW = 0; iW < nW; iW++)
        {
            A = Amin + iA * Ares;
            W = Wmin + iW * Wres;

            calcLyapunov(&L1, &L2, x0, y0, t0,
                dt, dT, Trun, tTransient, A, W, n_system, B_system, a_system);

            AVec.push_back(A);
            WVec.push_back(W);
            L1Vec.push_back(L1);
            L2Vec.push_back(L2);


            if (AVec.size() == nWriteInterval)
            {

                cout << "A=" << A << ", W=" << W << ";   ";
                cout << "Writing " << nWriteInterval << " rows of results ...\n";

                writeFile(wtireFileName, AVec, WVec,
                    L1Vec, L2Vec);

                AVec.clear();
                WVec.clear();
                L1Vec.clear();
                L2Vec.clear();

            }
        }
    }

    cout << "A=" << A << ", W=" << W << ";   ";
    cout << "Writing " << AVec.size() << " rows of results ...\n\n";

    writeFile(wtireFileName, AVec, WVec,
        L1Vec, L2Vec);

    tend = time(0);

    cout << "---------------------------------------------------\n";
    cout << "Elapsed time: " << difftime(tend, tstart) << " second(s)." << endl;


    cout << "Finish!\n";
    return 0;

}


void calcLyapunov(double* L1out, double* L2out, double x0, double y0, double t0,
    double dt, double dT, double Trun, double tTransient, double A, double W,
    int n_system, double B_system, double a_system)

{
    double d0 = 1.e-8;
    double dx, dy, d1;
    double x0Dev, y0Dev;

    double xOut, yOut;
    double xOutDev, yOutDev;
    double trJmean, trJTemp;

    double L1, L2, L1mean, L2mean;

    L1mean = 0;
    L2mean = 0;
    trJmean = 0;

    rk4vec(&xOut, &yOut, &trJTemp, t0, x0, y0, dT, tTransient,
        A, W, n_system, B_system, a_system);     //Transient Response

    t0 = t0 + tTransient;
    x0 = xOut; y0 = yOut;

    x0Dev = x0;
    y0Dev = y0 - d0;

    int nSteps = int(Trun / dT);

    for (int i = 1; i <= nSteps; i++)
    {
        rk4vec(&xOutDev, &yOutDev, &trJTemp, t0, x0Dev, y0Dev, dt, dT,
            A, W, n_system, B_system, a_system);
        rk4vec(&xOut, &yOut, &trJTemp, t0, x0, y0, dt, dT,
            A, W, n_system, B_system, a_system);

        t0 = t0 + dT;
        x0 = xOut; y0 = yOut;

        dx = xOut - xOutDev;
        dy = yOut - yOutDev;
        d1 = sqrt(dx * dx + dy * dy);
        L1 = log(d1 / d0) / dT;
        L1mean += L1;

        x0Dev = x0 - (d0 / d1) * dx;
        y0Dev = y0 - (d0 / d1) * dy;

        trJmean += trJTemp;
    }

    *L1out = L1mean / nSteps;
    *L2out = (trJmean - L1mean) / nSteps;
}

void rk4vec(double* xOut, double* yOut, double* TraceJ, double t0, double x0, double y0, double dt, double T,
    double A, double W, int n_system, double B_system, double a_system)
{
    //****************************************************************************


    double x = x0;
    double y = y0;

    double ti = t0;

    double k1, k2, k3, k4;
    double L1, L2, L3, L4;

    double trJ = 0; // local variable: trJ = B*besselj(n,x+a);


    int n = int(T / dt);

    for (int i = 0; i < n; i++)
    {
        /*
        k1 = dx_dt(ti, x, y);
        L1 = dy_dt(ti, x, y);


        k2 = dx_dt(ti + dt / 2, x + dt * k1 / 2, y + dt * L1 / 2);
        L2 = dy_dt(ti + dt / 2, x + dt * k1 / 2, y + dt * L1 / 2);

        k3 = dx_dt(ti + dt / 2, x + dt * k2 / 2, y + dt * L2 / 2);
        L3 = dy_dt(ti + dt / 2, x + dt * k2 / 2, y + dt * L2 / 2);

        k4 = dx_dt(ti + dt, x + dt * k3, y + dt * L3);
        L4 = dy_dt(ti + dt, x + dt * k3, y + dt * L3);

        x = x + 1. / 6. * dt * (k1 + 2. * k2 + 2. * k3 + k4);
        y = y + 1. / 6. * dt * (L1 + 2. * L2 + 2. * L3 + L4);
        */


        k1 = y;
        L1 = -x + B_system * y * _jn(n_system, x + a_system) +
            A * sin(W * ti);

        k2 = y + dt * L1 / 2;
        L2 = -(x + dt * k1 / 2) + B_system * (y + dt * L1 / 2) * _jn(n_system, (x + dt * k1 / 2) + a_system) +
            A * sin(W * (ti + dt / 2));

        k3 = y + dt * L2 / 2;
        L3 = (-(x + dt * k2 / 2) + B_system * (y + dt * L2 / 2) * _jn(n_system, (x + dt * k2 / 2) + a_system)) +
            A * sin(W * (ti + dt / 2));

        k4 = y + dt * L3;
        L4 = (-(x + dt * k3) + B_system * (y + dt * L3) * _jn(n_system, (x + dt * k3) + a_system)) +
            A * sin(W * (ti + dt));

        x = x + 1. / 6. * dt * (k1 + 2. * k2 + 2. * k3 + k4);
        y = y + 1. / 6. * dt * (L1 + 2. * L2 + 2. * L3 + L4);

        ti = ti + dt;

        trJ += B_system * _jn(n_system, x + a_system);

    }

    *xOut = x;
    *yOut = y;
    *TraceJ = trJ / n;

}

/*
double dx_dt(double t, double x, double y)
{
    return y;
}

double dy_dt(double t, double x, double y)
{
    //return -x + muVandepol*y*(1-x*x) + AVandepol*sin(wVandepol*t);
    return -x + B_system* y * _jn(n_system, x + a_system)+AA*sin(ww*t);

}
*/


void writeFile(string fileName, vector <double> AVec, vector <double> WVec,
    vector <double> L1Vec, vector <double> L2Vec)
{
    int nRows = AVec.size();

    //double Lmax, Lmin;

    ofstream fileOut;
    fileOut.open(fileName, std::ios_base::app);

    for (int i = 0; i < nRows; i++)
    {

        //Lmax =max(L1Vec[i], L2Vec[i]);
        //Lmin = min(L1Vec[i], L2Vec[i]);

        fileOut << AVec[i] << "," << WVec[i] << ","
            << L1Vec[i] << "," << L2Vec[i] << endl;
    }

    fileOut.close();
}

