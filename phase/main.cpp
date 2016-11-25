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
    cout << "Lattice orientation: " << spindirection << endl;
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
    // initializing randomness
    std::random_device rd;
    std::mt19937_64 gen(rd());
    // Set up the uniform distribution for x \in [[0, 1]
    std::uniform_real_distribution<double> RandomNumberGenerator(0.0,1.0);

    // difference in energies, made accessible
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
    cout << "Is file open?: " << outfile.is_open() << endl;
    cout << filename_data << endl;

    ofstream E_file;
    E_file.open(filename_E);
    E_file << "#E sigma_E" << endl;

    // here be dragons - montecarlo of metropolis algorithm
    for (int cycle=1; cycle <= MCcycles; cycle++){

        // flipping one spins, trying to get a lower state
        for (int x = 0; x<length; x++){
            for (int y=0; y<length; y++){
                // deciding on random position in the matrix
                int i = (int) (RandomNumberGenerator(gen)*(double)length);
                int j = (int) (RandomNumberGenerator(gen)*(double)length);
                // finding the energy differential of this spin's position
                int dE = 2*spins(i,j)*(spins( periodic(i,length, 1), j)
                                      +spins( periodic(i,length,-1), j)
                                      +spins( i, periodic(j,length, 1))
                                      +spins( i, periodic(j,length,-1)) );
                if (RandomNumberGenerator(gen)<=diffE(dE+8)){
                    spins(i,j) *= -1.0;
                    Energy += (double) dE;
                    Magnet += 2*spins(i, j);
                    accepted += 1;
                } // adds the new energy and magnetization in case this
            }     //has lower energy than the previous.
        }
        exval(0) += Energy;
        exval(1) += Energy*Energy;
        exval(2) += fabs(Magnet);
        exval(3) += Magnet*Magnet;

        if (((1000*cycle) % (MCcycles) == 0) || (cycle == MCcycles)){
            outfile << cycle;
            for (int i=0; i<4; i++){
                outfile << " " << exval(i)/cycle;
            }
            outfile << " " << accepted << endl;
        }
        if (cycle >= threshold){
            double norm = 1./((double) cycle);
            double en_it = exval(0);
            double en2_it = exval(1);

            en_it *= norm;
            en2_it *= norm;

            E_file << Energy << " " <<
                      (en2_it - ( en_it*en_it) )
                   << endl;
        }
    }
    // parameters averaged over montecarlocycles and spins
    exval_per_mc = exval / MCcycles;
    exval_per_mc_per_spin = exval_per_mc / (length*length);
    // print mean parameters per particle
    cout << "Mean energy per MC cycle:         E     = " <<
            exval_per_mc(0) << endl;
    cout << "Mean energy per MC cycle:         E**2  = " <<
            exval_per_mc(1) << endl;
    cout << "Mean magnetization per MC cycle: |M|    = " <<
            exval_per_mc(2) << endl;
    cout << "Mean magnetization per MC cycle: |M**2| = " <<
            exval_per_mc(3) << endl;
    cout << "Specific heat:                    Cv    = " <<
            (exval_per_mc(1) - exval_per_mc(0)*exval_per_mc(0))/
            (temp*temp)<< endl;
    cout << "Susceptibility:                   X     = " <<
            (exval_per_mc(3) - exval_per_mc(2)*exval_per_mc(2))/
            (temp) << endl;

    // closing
    outfile.close();
    E_file.close();

}

int main(int argc, char *argv[]){
    // engage MPI Magic
//    MPI_Init(NULL, NULL); //being conscious about environment is good
//    int numberOfProcesses;// get no. of processes
//    MPI_Comm_size(MPI_COMM_WORLD, &numberOfProcesses);
//    int processID;        // get rank of process
//    MPI_Comm_rank(MPI_COMM_WORLD, &processID);

    double temp=2., finaltemp=2.3; // as per exercise
    int L[5], MCcycles=1e+5, MCthreshold=1.3e+4; // yes, yes, i know.
    string datafilename, matrixargument="up", energyfilename;

    L[0] = 20; L[1] = 40; L[2] = 60; L[3] = 100; L[4] = 140;

    for (int i=0; i<=4; i++){
        // going over all the sizes of L

        for (double itemp=temp; itemp<=finaltemp; itemp+=0.0125){
            cout << endl;
            // should make around 8 examples

            imat spinsystem(L[i],L[i]);     // needs to reset the matrix
            double Energy = 0, Magnet = 0;  //and its values
            makeMatrix(spinsystem, L[i], matrixargument, Energy, Magnet);

            cout << "MCcycle tot: " << MCcycles << endl;
            cout << "Lattice len: " << L[i] << endl;
            cout << "Temperature: " << itemp << endl;
            cout << "Spin orient: " << matrixargument << endl;
            datafilename = "..\\data\\L"+to_string(L[i])
                                        +"\\Prob_L"+to_string(L[i])
                                  +"_mc"+to_string(MCcycles)
                                   +"_T"+to_string((int) (100*itemp))
                                                        // scaled
                                +"_spin"+matrixargument + ".dat";
            energyfilename = "..\\data\\L"+to_string(L[i])+
                                        "\\Energyprob_L"+to_string(L[i])
                                      +"_mc"+to_string(MCcycles)
                                       +"_T"+to_string((int) (100*itemp))
                                                        // scaled
                                    +"_spin"+matrixargument + ".dat";
            // probability run
            ising_prob(datafilename, energyfilename, spinsystem, L[i],
                       Energy, Magnet, itemp, MCcycles, MCthreshold);
        }
    }
//    MPI_Finalize();

}
