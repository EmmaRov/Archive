
#include <iostream>
#include <armadillo>
#include <vector>

class Particle
{
public:
    Particle(int c, int m, arma::vec r, arma::vec v);

    int c;
    int m;
    vec r;
    vec v;
    
    vec kv1;    vec kr1;
    vec kv2;    vec kr2;
    vec kv3;    vec kr3;
    vec kv4;    vec kr4;
    
};

vec forces(vector<Particle>& particles, int index, double t);

int RK4(vector<Particle>& particles, double t, double h);

int ForwardEuler(vector<Particle>& particles, double t, double h);


