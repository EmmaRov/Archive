import numpy as np
from tqdm import trange
import matplotlib.pyplot as plt
import time
import tabulate


#We create .npy files instead of .txt because these are much faster to load.
#But we only have to do this once, so we comment it out. From now on we load the .npy files
"""
txtdata = np.loadtxt('spectrum_seed70_600nm_3000nm.txt')
np.save('spectrumbutbetter.npy', txtdata)
"""

data = np.transpose(np.load('spectrumbutbetter.npy'))
# print(np.shape(data))
#(2, 10000000)

"""
txtsigma = np.loadtxt('sigma_noise.txt')
np.save('noisebutbetter.npy', txtsigma)
"""

noise = np.transpose(np.load('noisebutbetter.npy'))
# print(np.shape(noise))
#(2, 10000000)


# plt.plot(data[0],data[1])
# plt.show()

k = const.k_B               #Boltzmann constant
u = const.m_p               #atomic mass unit
c = const.c                 #speed of light

#spectral lines of common gasses
O2a = 632   ;O2b = 690   ;O2c = 760   ;H2Oa = 720  ;H2Ob = 820  ;H2Oc = 940
CO2a = 1400 ;CO2b = 1600 ;CH4a = 1660 ;CH4b = 2200 ;COa = 2340   ;N2Oa = 2870

#masses of common gasses
O2mass = 32*u   ;H2Omass = 18*u ;CO2mass = 44*u ;CH4mass = 16*u ;COmass = 28*u  ;N2Omass = 44*u


lines = np.array([O2a , O2b, O2c, H2Oa, H2Ob, H2Oc, CO2a, CO2b, CH4a, CH4b, COa, N2Oa])
masses = np.array([O2mass, O2mass, O2mass, H2Omass, H2Omass, H2Omass, CO2mass, CO2mass, CH4mass, CH4mass, COmass, N2Omass])
titles = [r'$O_2$ 632nm spectre', r'$O_2$ 690nm spectre', r'$O_2$ 760nm spectre', r'$H_2O$ 720nm spectre', r'$H_2O$ 820nm spectre', r'$H_2O$ 940nm spectre',
r'$CO_2$ 1400nm spectre', r'$CO_2$ 1600nm spectre', r'$CH_4$ 1660nm spectre', r'$CH_4$ 2200nm spectre', r'CO 2340nm spectre', r'$N_2O$ 2870nm spectre']


#each instance of this class is going to be an absorption line we are looking for.
class findspecline:
    def __init__(self, lambda_0, mass, data):
        #the wavelength we are centered around
        self.lambda_0 = lambda_
        #mass of the molecule we are looking for
        self.mass = mass
        #maximum maximum redshift
        self.maxshift = (10e3*self.lambda_0)/c + self.lambda_0
        #maximum blueshift
        self.minshift = -(10e3*self.lambda_0)/c + self.lambda_0
        #first index of data that we look at
        self.minrange = np.argmin(abs(data[0]-self.minshift))
        #last index of data that we look at
        self.maxrange = np.argmin(abs(data[0]-self.maxshift))
        #the range of indexes that we look at
        self.range = data[0,self.minrange:self.maxrange]
        #all the data
        self.data = data


        #Shows only the segment around lambda_0 without the fitted line.
    def showrange(self):

        plt.plot(self.range,self.data[1,self.minrange:self.maxrange])
        plt.show()

    """
    the next three methods (update-*) all have the same structure.
    First we define sigma (note that this sigma is not related to the noise) which will be used in the model
    Next we define 2d-arrays for evaluating all the wavelengths in the segment while also varying the parameter.
    Then we produce our models (F-*) that we comapare to the data with chi^2-minimization (*_chi).
    Finally we find which model in F-* that fitted best and keep the index of the best value of that paremeter as best-*
    """

    #Varying lambda while holding the two other parameters constant
    def updatelambda(self):
        lsigma = self.lambda_0/c * np.sqrt((k*self.temp[self.besttemp])/self.mass)
        Rl, Lambdav = np.meshgrid(self.range, self.lambdav)
        Fl = 1 + (self.F_min[self.bestFmin]-1) * np.exp(- ((Rl - Lambdav)**2)/(2*lsigma**2))
        l_chi = np.sqrt(np.sum(((self.Data - Fl)/noise[1,self.minrange:self.maxrange])**2, axis=1))
        self.bestlambda = np.where(l_chi == np.min(l_chi))[0][0]

    #Varying the temp while holding the two other parameters constant
    def updatetemp(self):
        tsigma = self.lambda_0/c * np.sqrt((k*self.temp)/self.mass)
        Rt , S = np.meshgrid(self.range, tsigma)
        F_temp = 1 + (self.F_min[self.bestFmin]-1) * np.exp(- ((Rt-self.lambdav[self.bestlambda])**2)/(2*S**2))
        temp_chi = np.sqrt(np.sum(((self.Data - F_temp)/noise[1,self.minrange:self.maxrange])**2, axis=1))
        self.besttemp = np.where(temp_chi == np.min(temp_chi))[0][0]

    #Varying F_min while holding the two other parameters constant
    def updateFmin(self):
        fsigma = self.lambda_0/c * np.sqrt((k*self.temp[self.besttemp])/self.mass)
        Rf, F_MIN = np.meshgrid(self.range, self.F_min)
        F_depth = 1 + (F_MIN-1) * np.exp(-((Rf - self.lambdav[self.bestlambda])**2)/(2*fsigma**2))
        f_chi = np.sqrt(np.sum(((self.Data - F_depth)/noise[1,self.minrange:self.maxrange])**2, axis=1))
        self.bestFmin = np.where(f_chi == np.min(f_chi))[0][0]


    #Finds the best fitted model after n iterations of the Gibbs optimizer with initial guesses.
    def bestmodel(self, N, n, title, show=True):
        self.N = N     #grid resolution for parameters
        self.n = n     #number of iterations of optimizer
        self.title = title
        self.Data = np.tile(self.data[1,self.minrange:self.maxrange], ((self.N,1)))

        # Range of each parameter space
        self.temp = np.linspace(150,450,self.N)
        self.F_min = np.linspace(0.7,0.9,self.N)
        self.lambdav = np.linspace(self.minshift, self.maxshift, self.N)

        # Indexes of initial guess of parameters
        self.besttemp = int(self.N/2)
        self.bestFmin = 0
        self.bestlambda = np.argmin(abs(self.lambdav - self.lambda_0))

        #for each iteration of this loop we update the "best" indexes for each of the parameters.
        #After enough iterations, the three parameters converge, these are the values used in the final model.
        for i in range(self.n):
            self.updatelambda()
            self.updatetemp()
            self.updateFmin()

        #final model
        self.bestsigma = self.lambda_0/c * np.sqrt((k*self.temp[self.besttemp])/self.mass)
        self.F = 1 + (self.F_min[self.bestFmin] - 1) * np.exp(- ((self.range - self.lambdav[self.bestlambda])**2)/(2*self.bestsigma**2))


        print(f'When trying to find wavelengt {self.lambda_0}, the best model had the following parameters:')
        print(f'Wavelenght at minima: {self.lambdav[self.bestlambda]}\n Tempreature: {self.temp[self.besttemp]} \n Flux at minima: {self.F_min[self.bestFmin]}')
        print()

        if show:
            plt.plot(self.range,self.data[1,self.minrange:self.maxrange], label='Data')
            plt.plot(self.range,self.F, label='Fitted line')
            plt.xlabel('Wavelenght [nm]', fontsize=16)
            plt.ylabel('Flux [-]', fontsize=16)
            plt.title(self.title, fontsize=16)
            plt.legend(prop={'size':14})
            # plt.savefig(f'{self.lambda_0}_specter.png')
            plt.show()



#test function to prove concept and test accuracy
def test(lambda0, T, fmin, speed, mass, show=True): #all the ingredients in order to make a lineprofile
    dummylambda = lambda0 #the wavelength of our fabricated signal without shifts
    maxshiftdummy = (10e3*dummylambda)/c + dummylambda #max blueshift
    minshiftdummy = -(10e3*dummylambda)/c + dummylambda #max redshift

    minrange = np.argmin(abs(data[0]-minshiftdummy)) #first index of data to look at
    maxrange = np.argmin(abs(data[0]-maxshiftdummy)) #last index of data to look at
    drange = data[0,minrange:maxrange]               #range of data to look at

    dummyshift = (speed*dummylambda)/c + dummylambda  #shifting the fabricated dummy wavelength
    dummysigma = dummylambda/c * np.sqrt(k*T/(mass))  #Finding standard deviation for dummy lineprofile

    #make the lineprofile
    dummyF = 1 + (fmin-1) * np.exp(- ((drange - dummyshift)**2)/(2*dummysigma**2))
    disortedF = np.zeros(len(dummyF)) #empty array to fill with distorted values of 'dummyF'

    #we want the lineprofile to be scrambled in the same way as the real data is,
    #so that the chi-minimization in the class works with this dummy
    for i in range(len(disortedF)):
        disortedF[i] = np.random.normal(dummyF[i],noise[1,minrange+i])

    #shows the lineprofile we made and the data we will feed the main algorithm
    if show:
        plt.plot(drange, disortedF)
        plt.plot(drange, dummyF)
        # plt.savefig('testline.png')
        plt.show()

    return np.array([drange,disortedF])

#is this elegant and generalized? no. Do i have the patience to change this? also no.
def runtest(N,l):
    diff = np.zeros((N,3))
    dummyshift = (4e3*l)/c + l
    for i in range(N):
        dummydata = test(l, 355, 0.8, 4e3, 20*u, False)
        dummy = findspecline(l, 20*u, dummydata)
        dummy.bestmodel(1000, 10, 'Test spectrum', False)
        finaltemps = dummy.temp[dummy.besttemp]
        finalfmins = dummy.F_min[dummy.bestFmin]
        finallambdas = dummy.lambdav[dummy.bestlambda]
        diff[i] = np.array([abs(dummyshift-finallambdas)/dummyshift, abs(finaltemps-355)/355, abs(finalfmins-0.8)/0.8])
    print(f'mean relative error for lambda {np.mean(diff[:,0])}')
    print(f'mean relative error for temp {np.mean(diff[:,1])}')
    print(f'mean relative error for F_min {np.mean(diff[:,2])}')

# runtest(50,800)

t1 = time.time()
for i in range(len(lines)):
    O = findspecline(lines[i], masses[i], data)
    O.bestmodel(1000, 10, titles[i], False)
    # finaltemps[i] = O.temp[O.besttemp]
    # finalfmins[i] = O.F_min[O.bestFmin]
    # finallambdas[i] = O.lambdav[O.bestlambda]
t2 = time.time()

# print(f'The simulation ran in {t2-t1} seconds')

def latextable():
    finaltemps = np.zeros(len(lines))
    finallambdas = np.zeros(len(lines))
    finalfmins = np.zeros(len(lines))
    table = np.transpose(np.array([titles, finallambdas, finaltemps, finalfmins], dtype = 'object'))
    print(tabulate.tabulate(table, tablefmt = 'latex_raw'))

# latextable()
