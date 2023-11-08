import numpy as np
import matplotlib.pyplot as plt
import scipy.io.wavfile as wavfile
from time import time


# hvilken = input('Hva skal du kjøre?')
hvilken = 'anton'

#kutter signalet fra første punkt lyden starter
def cuttit(double_signal):
    signal = double_signal[:,0]
    start = np.argwhere(abs(signal)>1000)[0][0] #støy er under 350
    cut_signal = signal[start:]
    return cut_signal

#produserer .wav fil for å prudusere 1sek puls med bestemt frekvens
def lagelyd(Hz):
    dur = np.linspace(0, 1, fs_rate)
    data = np.sin(dur*Hz*2*np.pi)
    wavfile.write(f'{Hz}Hz_1s_impuls.wav', fs_rate, data.astype(np.float32))

#tar input-variabelen 'hvilken' og preparerer dataene
def velg(hvilken):
    if hvilken == 'anton':
        #Anton som sier hei
        fs_rate, doubble_signal = wavfile.read("proffessoranton.wav")
        signal = doubble_signal[90000:152000,0]
    if hvilken == 'test':
        #test
        t0 = np.linspace(0,10,100000)
        signal = np.sin(2*np.pi*200*t0)
        dt = t0[1] - t0[0]
        fs_rate = 1/dt
    if hvilken == '147':
        fs_rate, double_signal = wavfile.read('147Hz.wav')
        signal = cuttit(double_signal)
    if hvilken == '397':
        fs_rate, double_signal = wavfile.read('397Hz.wav')
        signal = cuttit(double_signal)
    if hvilken == '550':
        fs_rate, double_signal = wavfile.read('550Hz.wav')
        signal = cuttit(double_signal)
    if hvilken == '792':
        fs_rate, double_signal = wavfile.read('792Hz.wav')
        signal = cuttit(double_signal)
    if hvilken == '1612':
        fs_rate, double_signal = wavfile.read('1612Hz.wav')
        signal = cuttit(double_signal)

    return fs_rate, signal


#Gjør wavelet analyse med Morlet wavelet og plotter
#Gjør også noen tids-tester for å sammenlikne med andre metoder
def wavelet_analyse(K, S, t, omega):

    #Frekvensene vi er interresert i å analysere
    omega_a = np.linspace(100, 750, 100)*2*np.pi

    #meshgrid-aksene vi skal ha resultatet på
    aa, tt = np.meshgrid(omega_a, t, indexing = 'ij')

    t2 = time()

    #Fourier transformert Morlet wavelet
    Psi = 2*(np.exp(-(K*(omega-aa)/aa)**2) - np.exp(-K**2) * np.exp(-(K*omega/aa)**2))

    t3 = time()

    print(f'time for fouriertransform: {t1-t0} \n time for generating Psi matrix: {t3-t2}')

    t4 = time()
    #resultatmatrise som følge av konvolusjonsteoremet
    result = abs(np.fft.ifft(S*Psi))
    t5 = time()

    print(f'time to produce result matrix: {t5-t4}')


    plt.contourf(tt, aa/(2*np.pi), result, cmap = 'hot', levels = 25)
    plt.xlabel('time [s]')
    plt.ylabel('frequency [Hz]')
    plt.show()



#velger en fil og gjør dataene klare
fs_rate, signal = velg(hvilken)

dt = 1/fs_rate
N = len(signal)
T = dt*N
t = np.linspace(0,T,N)

#Parameter for bredden til morlet wavelet
K = 6

t0 = time()

S = np.fft.fft(signal)
omega = np.fft.fftfreq(N, dt)*2*np.pi

t1 = time()

wavelet_analyse(K, S, t, omega)
