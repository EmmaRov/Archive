#include "RK4.hpp"

using namespace std;
using namespace arma;


//definition of constructor
Particle::Particle(int charge, int mass, vec position, vec velocity)
{
    c = charge;
    m = mass;
    r = position;
    v = velocity;
}


//test force that only accounts for interactve forces between particles
vec forces(vector<Particle>& particles, int index, double t)
{
    vec tot_force = {0., 0., 0.};
    //Particle-particle interaction
    for(int i = 0; i<particles.size(); i++)
    {   
        if (i != index)
        {
            vec R = particles.at(index).r - particles[i].r ;
            tot_force += particles[i].c * R/(pow(arma::norm(R), 3));
        }
    }    

    //Lorentz force (could add eletric field, but this is only a test)
    vec vel = particles.at(index).v;
    vec B_force = {vel.at(1), vel.at(0), 0.}; //result of the cross product if B_0 = 1
    tot_force += B_force;

    return tot_force;
}

//Calling function evolves the system 1 timestep
//t is what time, h is length of timestep
int RK4(vector<Particle>& particles, double t, double h)
{
    //copying vector of state before we begin time-evolution
    vector<Particle> safe_particles = particles;
    //declaring all the k-values
    

    //first step in RK4
    for(int i = 0; i<particles.size(); i++)
    {
        particles[i].kv1 = forces(particles, i, t)/particles[i].m;
        particles[i].kr1 = particles[i].v;

        particles[i].v = safe_particles[i].v + particles[i].kv1*h/2.;
        particles[i].r = safe_particles[i].r + particles[i].kr1*h/2.;
    }

    //second step of RK4
    for(int i = 0; i<particles.size(); i++)
    {
        particles[i].kv2 = forces(particles, i, t+(h/2.))/particles[i].m;
        particles[i].kr2 = particles[i].v;

        particles[i].v = safe_particles[i].v + particles[i].kv2*h/2.;
        particles[i].r = safe_particles[i].r + particles[i].kr2*h/2.;
    }

    //third step of RK4
    for(int i = 0; i<particles.size(); i++)
    {
        particles[i].kv3 = forces(particles, i, t+(h/2.))/particles[i].m;
        particles[i].kr3 = particles[i].v;
        
        particles[i].v = safe_particles[i].v + particles[i].kv3*h;
        particles[i].r = safe_particles[i].r + particles[i].kr3*h;
    }

    //fourth step of RK4
    for(int i = 0; i<particles.size(); i++)
    {
        particles[i].kv4 = forces(particles, i, t+h)/particles[i].m;
        particles[i].kr4 = particles[i].v;
    
        
        //Final update, now we have evolved one timestep
        particles[i].v = safe_particles[i].v + 
        h*(particles[i].kv1 + 2.*particles[i].kv2 + 2.*particles[i].kv3 + particles[i].kv4)/6.;
        particles[i].r = safe_particles[i].r + 
        h*(particles[i].kr1 + 2.*particles[i].kr2 + 2.*particles[i].kr3 + particles[i].kr4)/6.;
    }

    return 0;
}

int ForwardEuler(vector<Particle>& particles, double t, double h)
{
    vector<Particle> safe_particles = particles;

    for(int i = 0; i<particles.size(); i++)
    {
        vec a = forces(particles, i, t)/safe_particles[i].m;
        particles[i].v = safe_particles[i].v + h*a;
        particles[i].r = safe_particles[i].r + h*safe_particles[i].v;
    }
    

    return 0;
}


int main()
{
    vector<Particle> particles;
    vec zero = {0.,0.,0.};
    vec one = {1.,0.,1.};
    vec two = {1.,-1.,0.};
    vec three = {0.,1.,0.};
    Particle p1(1, 1, one, two);
    Particle p2(1, 1, two, three);
    Particle p3(1, 1, three, one);
    particles.push_back(p1);
    particles.push_back(p2);
    particles.push_back(p3);

    
    int N = 1e5;
    double T = 1.;
    double h = T/N;
    vec t = linspace(0,T,N);

    mat log_1(N, 3, fill::zeros);
    mat log_2(N, 3, fill::zeros);
    mat log_3(N, 3, fill::zeros);
    
    for (double i= 0; i<N; i++)
    {
        log_1.row(i) = particles[0].r.t();
        log_2.row(i) = particles[1].r.t();
        log_3.row(i) = particles[2].r.t();
        RK4(particles, t(i), h);
    }
    log_1.save("part_1_data.txt", raw_ascii);
    log_2.save("part_2_data.txt", raw_ascii);
    log_3.save("part_3_data.txt", raw_ascii);

    return 0;
}

