#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>
#include <armadillo>
#include <random>

using namespace std;
using namespace arma;
// inline function for periodic boundary conditions
inline int periodic(int i, int limit, int add){
    return (i+limit+add) % (limit);
}

// function to initialize energy, spin matrix and magnetization
void makeMatrix(imat &spins, int matrixlength, string spindirection,
                double &E, double &M){
    // setting up spin matrix and initial magnetization
    if (spindirection == "up"){
        for (int i=0; i<matrixlength; i++){
            for (int j=0; j<matrixlength; j++){
                spins(i,j) = 1;
            }
        }
    } else if (spindirection == "down"){
        for (int i=0; i<matrixlength; i++){
            for (int j=0; j<matrixlength; j++){
                spins(i,j) = -1;
            }
        }
    } else if (spindirection == "random"){
        for (int i=0; i<matrixlength; i++){
            for (int j=0; j<matrixlength; j++){
                double r = rand() / ((double) numeric_limits<int>::max());
                if (r <= 0.5){
                    spins(i,j) = 1;
                } else {
                    spins(i,j) = -1;
                }
            }
        }
    } else { cout << "Incorrect spindirection setting" << endl; }
    for (int i=0; i<matrixlength; i++){
        for (int j=0; j<matrixlength; j++){
            M += (double) spins(i,j);
            E -= (double) spins(i,j)*
                            (spins(i ,periodic(j,matrixlength,-1))
                           + spins(periodic(i,matrixlength,-1), j));
        }
    }
}


void ising_prob(string filename_data, string filename_E,
                imat spins, int length, double Energy, double Magnet,
                double temp, int MCcycles, int threshold){

    int accepted=0;
    double beta = 1./temp;
    std::random_device rd;
    std::mt19937_64 gen(rd());
    // Set up the uniform distribution for x \in [[0, 1]
    std::uniform_real_distribution<double> RandomNumberGenerator(0.0,1.0);

    vec diffE = zeros<mat>(17);
    for( int de =-8; de <= 8; de+=4) diffE(de+8)
              = exp(-de/temp);
    vec exval                   = zeros<vec>(4);
    vec exval_per_mc            = zeros<vec>(4);
    vec exval_per_mc_per_spin   = zeros<vec>(4);

    // start to datafile
    ofstream outfile;
    outfile.open(filename_data);
    outfile << "#msc E E^2 |M| M^2 accepted" << endl;

    ofstream E_file;
    E_file.open(filename_E);
    E_file << "#E" << endl;

    // here be dragons - montecarlo of metropolis algorithm
    for (int cycle=1; cycle < MCcycles; cycle++){
        //metropolis(matrix, beta, length, Energy, Magnet, accepted);

        for (int x = 0; x<length; x++){
            for (int y=0; y<length; y++){
                int i = (int) (RandomNumberGenerator(gen)*(double)length);
                int j = (int) (RandomNumberGenerator(gen)*(double)length);
                int dE = 2*spins(i,j)*(spins( periodic(i,length, 1), j)
                                      +spins( periodic(i,length,-1), j)
                                      +spins( i, periodic(j,length, 1))
                                      +spins( i, periodic(j,length,-1)) );
                if (RandomNumberGenerator(gen)<=diffE(dE+8)){
                    spins(i,j) *= -1.0;
                    Energy += (double) dE;
                    Magnet += 2*spins(i, j);
                }
            }
        }

        exval(0) += Energy;
        exval(1) += Energy*Energy;
        exval(2) += fabs(Magnet);
        exval(3) += Magnet*Magnet;

        if (((100000*cycle) % (MCcycles) == 0) || (cycle == MCcycles)){
            outfile << cycle;
            for (int i=0; i<4; i++){
                outfile << " " << exval(i)/cycle;
            }
//            cout << "Iterative stuff: " << endl;
//            cout << "<E>    = " << exval(0)/cycle << endl;
//            cout << "<E**2> = " << exval(1)/cycle << endl;
//            cout << "<M>    = " << exval(2)/cycle << endl;
//            cout << "<M**2> = " << exval(3)/cycle << endl;
            outfile << " " << accepted << endl;
        }
        if (cycle >= threshold){
            E_file << Energy << endl;
        }
    }
    // parameters averaged over montecarlocycles and spins
    exval_per_mc = exval / MCcycles;
    exval_per_mc_per_spin = exval_per_mc / (length*length);

    // closing
    outfile.close();
    E_file.close();

    // print mean parameters per particle
    cout << "Mean energy per MC cycle:         E     = " <<
            exval_per_mc(0) << endl;
    cout << "Mean energy per MC cycle:         E**2  = " <<
            exval_per_mc(1) << endl;
    cout << "Specific heat:                    Cv    = " <<
            (exval_per_mc(1) - exval_per_mc(0)*exval_per_mc(0))/
            (temp*temp) << endl;
    cout << "Mean magnetization per MC cycle: |M|    = " <<
            exval_per_mc(2) << endl;
    cout << "Mean magnetization per MC cycle: |M**2| = " <<
            exval_per_mc(3) << endl;
    cout << "Susceptibility:                   X     = " <<
            (exval_per_mc(3) - exval_per_mc(2)*exval_per_mc(2))/
            (temp) << endl;
}

int main(int argc, char *argv[]){
    //srand(time(NULL)); // seeding of random numbers

    // engage
    double temp=1.;
    int L=2, MCcycles=1e+5, MCthreshold=1e+5;
    string datafilename, matrixargument="up", energyfilename;

    if (argc == 5){
        L               = atoi(argv[1]);
        MCcycles        = atoi(argv[2]);
        temp            = atof(argv[3]);
        datafilename    = argv[4];
    } if (argc == 6){
        L               = atoi(argv[1]);
        MCcycles        = atoi(argv[2]);
        temp            = atof(argv[3]);
        datafilename    = argv[4];
        matrixargument  = argv[5];
    }
    // set up system
    imat spinsystem(L,L);
    double Energy = 0, Magnet = 0;
    makeMatrix(spinsystem, L, matrixargument, Energy, Magnet);

    cout << "Starting Ising model simulation." << endl;

    // metropolis
    // probability run
    datafilename = "..\\data\\Prob_L"+to_string(L)
                          +"_mc"+to_string(MCcycles)
                           +"_T"+to_string((int) (100*temp)) // scaled
                        +"_spin"+matrixargument + ".dat";
    energyfilename = "..\\data\\Energyprob_L"+to_string(L)
                              +"_mc"+to_string(MCcycles)
                               +"_T"+to_string((int) (100*temp))// scaled
                            +"_spin"+matrixargument + ".dat";
    cout << "saved data to file path" << endl
         << datafilename << endl
         << energyfilename << endl;
    ising_prob(datafilename, energyfilename, spinsystem, L, Energy, Magnet,
               temp, MCcycles, MCthreshold);
}
