#include "MC.hpp"

using namespace std;
using namespace arma;


//Constructor
IsingModel::IsingModel(int lattice_size, double temperature)
{
    L = lattice_size;
    N = L*L;
    state.set_size(L,L);

    T = temperature;
    
    e_beta = exp(-1./T);

    //filling in values for probability distribution ratio map
    pmap[-8] = exp(-8)*e_beta;
    pmap[-4] = exp(-4)*e_beta;
    pmap[0] = e_beta;
    pmap[4] = exp(4)*e_beta;
    pmap[8] = exp(8)*e_beta;


}

//calling this method runs one Monte Carlo cycle
void IsingModel::MCMC_cycle()
{
    // unsigned int seed = chrono::system_clock::now().time_since_epoch().count();


    for(int i=0; i<N; i++)
    {
        mt19937 generator;

        generator.seed(seed+i);

        //Choose a uniformly random spin in the state lattice
        uniform_int_distribution<int> uni_dist_int(0,N-1);
        int spin_i = uni_dist_int(generator);

    
        //Metropolis algo for accept/reject:
        //details on method to find the ratio of probabilities can be found in the rapport

        int neighbour_energy = state((spin_i+L)%N) + state((spin_i+N-L)%N) + state((spin_i+1)%L,spin_i/L) + state((spin_i+L-1)%L,spin_i/L);
        int energy_change = state(spin_i)*neighbour_energy*2; 
        double prob_ratio = pmap[energy_change];


        double A = min(1., prob_ratio);

        uniform_real_distribution<double> uni_dist_real(0,1);
        double r = uni_dist_real(generator);

        if(r <= A)
        {
            state(spin_i) *= -1;
            total_energy += energy_change;
            total_magnetisation += 2*state(spin_i);
        }
    }
}


double IsingModel::sample_energy()
{
    // int total_energy = 0;

    return static_cast<double>(total_energy)/N;

}

double IsingModel::sample_magnetisation()
{

    return abs(static_cast<double>(total_magnetisation))/N;
}

//Could be part of constructor
void IsingModel::set_first_sample(string setting)
{
    if(setting == "mix")
    {
        mt19937 generator;

        generator.seed(seed);

        uniform_int_distribution<int> uni_dist_01(0,1);

        for(int i=0; i<N; i++)
        {
            state(i) = uni_dist_01(generator)*2 - 1;
        }

    }
    
    if(setting == "uniform")
    {
        state.fill(1);
    }

    //find initial total energy
    total_energy = 0;
    for(int i=0; i<N; i++)
    {
        total_energy += state(i)*(state((i+L)%N)+state((i+1)%L, i/L));
    }

    //find initial total magnetisation
    total_magnetisation = accu(state);
}

Mat<int> IsingModel::return_lattice()
{
    return state;
}


int main()
{
    IsingModel magnet(50, 200);
    magnet.set_first_sample("mix");
    int N_samples = 500;

    vec Energy(N_samples);
    vec Magnetisation(N_samples);

    for(int i=0; i<N_samples; i++)
    {
        magnet.MCMC_cycle();
        Energy(i) = magnet.sample_energy();
        Magnetisation(i) = magnet.sample_magnetisation();
    }

    Energy.save("Energy_samples.txt", raw_ascii);
    Magnetisation.save("Mag_samples.txt", raw_ascii);

    double average_E = mean(Energy);
    double average_M = mean(Magnetisation);
    double variance_E = var(Energy);
    double variance_M = var(Magnetisation);

    return 0;

}
