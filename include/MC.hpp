#include <iostream>
#include <armadillo>
#include <random>

class IsingModel
{
private:         

public:
    unsigned int seed = 593640;

    int L;                  
    int N;                  
    arma::Mat<int> state;
    double T;

    double e_beta;

    int total_energy;
    int total_magnetisation;

    std::map<int, double> pmap;

    IsingModel(int lattice_size, double temperature);
    
    void MCMC_cycle();

    double sample_energy();

    double sample_magnetisation();

    void set_first_sample(std::string setting);

    Mat<int> IsingModel::return_lattice();

};