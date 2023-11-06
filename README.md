# Some of my earlier projects
This repostitory is a collection of my previous projects that I am most happy about. Most of them are parts of  bigger projects that I have tackled in collaboration with others. For this reason they might not be very useful on their own, or general enough to be directly applied to other use. However they complete in their application of useful algorithms and data structure. 


## spectral_lines.py
This is the oldest program in this folder, and was written only about a year after I learned my first programming language; Python. The problem we had to solve, was to find spectral lines in a very big and noisy data-set. In short, I had to scan the 3-dimentional parameter space and use the chi-squared method to compare data with potential spectral line profiles. 
The problem was set up to scan the whole parameter space, but I wanted to experiment, so I decided to get creative, and invented my own algorithm (it probably already existed, but it was new to me). There are more details in the assosiated(paper name), but in short summary: Instead of scanning the whole parameter space, I make an initial guess of two parameters, say A and B and scan only the last parameter C, for the best fit. Then i fix C and B, and start scanning the A-space. Fix this and go to work on the B space. For every such scan i need only run over a 1D space, rather than 3D. This made for an impressive run-time, but at a loss of accuracy (but no loss of precision). Since we in this case could be rather confident in our initial guess, the algorithm held more or less the same accuracy as my peers who scanned the wole parameter space. 

As mentioned, there is a more detaild description of the algorithm and project()

## wavelet_analysis.py


## Jacobi.cpp


## RK4.cpp
RK4.cpp, with its assosiated header RK4.hpp is a part of a bigger project to which I contributed with a numerical integrator. In this program i apply the fourth order Runge-Kutta method of integration on particles flying around, affected by some force. In this spesific case, the force acting was dependent on the position of all other particles in the system. This made the implementation of RK4 a bit more tricky, as all of the particles had to be evolved in time in perfect tandem. The probelm was solved, and the program tested and proved to work perfectly when merged with the rest of the code structure (not included here, as most of it is not my code). 
RK4.cpp is tested in main() with three particles and a dimentionless force to show that it works. 

RK4_plotting.py is a short python script used for plotting the data it produces.

How to build and run:
$ g++ src/RK4.cpp -I include -o RK4.exe -larmadillo
$ ./RK4.exe
$ python3 RK4_plotting.py 

## MC.cpp
In this project we used Markov-chain Monte Carlo algorithm with the Metropolis algorithm to study the Ising-model. The program contains a IsingModel class, where instances of this class were run in parallel to sample millions of samples from the studied probability distributions. The class contains 4 relevant methods and a constructor. set_first_sample() lets you decide which state you start in, the following samples will be dependent on this choice, so a better start yields better results. sample_energy() and sample_magnetisation() returns the current state. And MC_cycle() runs a Monte Carlo cycle. In the end of the file is a main() function showcasing an example of how these methods can be used.

How to build and run:
$ g++ src/MC.cpp -I include -o MC.exe -larmadillo
$ ./MC.exe

There is also a python program that can be usde for some crude plotting (was initialy only used for testing) called MC_plotting.py  (remember to add that)