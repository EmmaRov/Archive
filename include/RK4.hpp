
#include <iostream>
#include <armadillo>
#include <vector>

class Particle
{
public:
    Particle(int charge, int mass, arma::vec position, arma::vec velocity);

    int c;
    int m;
    arma::vec r;
    arma::vec v;
    
    arma::vec kv1;    arma::vec kr1;
    arma::vec kv2;    arma::vec kr2;
    arma::vec kv3;    arma::vec kr3;
    arma::vec kv4;    arma::vec kr4;
    
};

arma::vec forces(std::vector<Particle>& particles, int index, double t);

int RK4(std::vector<Particle>& particles, double t, double h);

int ForwardEuler(std::vector<Particle>& particles, double t, double h);


