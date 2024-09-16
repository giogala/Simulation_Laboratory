/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#ifndef __System__
#define __System__

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <armadillo>
#include <stdlib.h> //exit
#include "particle.h"
#include "random.h"

using namespace std;
using namespace arma;

class System {

private:
    const int _ndim = 3;  // Dimensionality of the system
    bool _restart;        // Flag indicating if the simulation is restarted
    int _sim_type;        // Type of simulation (e.g., Lennard-Jones, Ising)
    int _npart;           // Number of particles
    int _nblocks;         // Number of blocks for block averaging
    int _nsteps;          // Number of simulation steps in each block
    int _nattempts;       // Number of attempted moves
    int _naccepted;       // Number of accepted moves
    double _temp, _beta;  // Temperature and inverse temperature
    double _rho, _volume; // Density and volume of the system
    double _r_cut;        // Cutoff radius for pair interactions
    double _delta;        // Displacement step for particle moves
    double _J, _H;        // Parameters for the Ising Hamiltonian
    double _t, _n;        // Parameters for equilibration: temperature of 1st sim and # of sims
    double *_ts;          // Pointer for equilibration: temperatures of sims
    int _k;               // Index of Equilibration sim
    string inp;           // directory with input files
    string out;           // directory with output files
    string inp_f;         // name of input file located in inp directory
    vec _side;            // Box dimensions
    vec _halfside;        // Half of box dimensions
    Random _rnd;          // Random number generator
    field <Particle> _particle; // Field of particle objects representing the system
    vec _fx, _fy, _fz;    // Forces on particles along x, y, and z directions
        // Properties
    int _nprop; // Number of properties being measured
    bool _measure_penergy, _measure_kenergy, _measure_tenergy;// Flags for measuring   ifferent energies
    bool _measure_temp, _measure_pressure, _measure_gofr;     // Flags for measuring   emperature, pressure, and radial dist. function
    bool _measure_magnet, _measure_cv, _measure_chi;          // Flags for measuring   agnetization, heat capacity, and susceptibility
    int _index_penergy, _index_kenergy, _index_tenergy;       // Indices for accessing     nergy-related properties in vec _measurement
    int _index_temp, _index_pressure, _index_gofr;            // Indices for accessing     emperature, pressure, and radial dist. function
    int* _index_t;
    int _index_magnet, _index_cv, _index_chi;                 // Indices for accessing     agnetization, heat capacity, and susceptibility
    int _n_bins;           // Number of bins for radial distribution function
    double _bin_size;      // Size of bins for radial distribution function
    double _vtail, _ptail; // Tail corrections for energy and pressure
    vec _block_av;         // Block averages of properties
    vec _global_av;        // Global averages of properties
    vec _global_av2;       // Squared global averages of properties
    vec _average;          // Average values of properties
    vec _measurement;      // Measured values of properties
public: // Function declarations
    int get_nbl();              // Get the number of blocks
    int get_nsteps();           // Get the number of steps in each block
    void SetInputf(string inpf);// Set the input file name
    void initialize(string inp_d,string out_d,string ran_d,string f_inp);          // Initialize system properties
    void initialize_properties();// Initialize properties for measurement
    void finalize();            // Finalize system and clean up
    void write_configuration(); // Write final system configuration to XYZ file
    void write_XYZ(int nconf);  // Write system configuration in XYZ format on the fly
    void write_velocities();    // Write final particle velocities to file
    void read_configuration();  // Read system configuration from file
    void initialize_velocities(int k);// Initialize particle velocities
    void step();                // Perform a simulation step
    void block_reset(int blk);  // Reset block averages
    void measure();             // Measure properties of the system
    void averages(int blk);     // Compute averages of properties
    double error(double acc, double acc2, int blk); // Compute error
    void move(int part);        // Move a particle
    bool metro(int part);       // Perform Metropolis acceptance-rejection step
    double pbc(double position, int i); // Apply periodic boundary conditions for coordinates
    int pbc(int i);             // Apply periodic boundary conditions for spins
    void Verlet();              // Perform Verlet integration step
    double Force(int i, int dim); // Calculate force on a particle along a dimension
    double Boltzmann(int i, bool xnew); // Calculate Boltzmann factor for Metropolis acceptance

};

#endif // __System__

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
