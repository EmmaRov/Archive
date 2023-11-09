import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt

#download and prepare data
data = sio.loadmat('data.mat')  

#coordinates 
x = data.get('x')
y = data.get('y')

#velocity component values
u = data.get('u')
v = data.get('v')

#limit between gas and liquid (found beforehand)
xit = data.get('xit')
yit = data.get('yit')

X = np.transpose(x) 
Y = np.transpose(y)
U = np.transpose(u)
V = np.transpose(v)

#find 'infinitesimals'
dx = x[0,1] - x[0,0]
dy = y[1,0] - y[0,0]


#Test function 
def sizecheck():
    try: 
        assert x.shape == (201,194), 'x matrix wrong size'
        assert y.shape == (201,194), 'y matrix wrong size'
        assert u.shape == (201,194), 'u matrix wrong size'
        assert v.shape == (201,194), 'v matrix wrong size'
        assert xit.shape == (1,194), 'xit vector wrong size'
        assert yit.shape == (1,194), 'yit vector wrong size'
        assert np.all(np.diff(x) == dx), 'x is not linearly spaced'
        assert np.all(np.diff(y, axis=0) == dy), 'y not linearly spaced'
        # assert np.all(y[-1,:] - y[0,:] == 100), 'nei'
    except AssertionError as msg:
        print(msg)
    

def speed():
    H = np.sqrt(U**2+V**2) #scalar field for speed

    #Plotting for speeds from 600 mm/s
    fig, axs = plt.subplots(2)
    axs[0].set_title('speed from 600 mm/s')
    div_plot1 = axs[0].contourf(X, Y, H, levels=np.linspace(600, 4200, 40),cmap='jet') 
    plt.colorbar(div_plot1, ax=axs[0])
    axs[0].plot(xit[0], yit[0],'k')
    axs[0].set_xlabel('x-axis [mm]')
    axs[0].set_ylabel('y-axis [mm]')

    #Plotting for speeds under 1100 mm/s
    axs[1].set_title('speed under 1100 mm/s')
    div_plot2 = axs[1].contourf(X,Y,H, levels=np.linspace(0, 1100, 60),cmap='jet') 
    plt.colorbar(div_plot2, ax=axs[1])
    axs[1].plot(xit[0], yit[0], 'k')
    axs[1].set_xlabel('x-axis [mm]')
    axs[1].set_ylabel('y-axis [mm]')

    fig.tight_layout()
    plt.show()  



#Function for plotting a square at given x coordinate 
def box(y2,y1):   

    x34 = X[34][0]
    x69 = X[69][0]
    plt.plot([x34,x69],[Y[0][y2],Y[0][y2]],'b')
    plt.plot([x69,x69],[Y[0][y2],Y[0][y1]], 'g')
    plt.plot([x34,x69],[Y[0][y1],Y[0][y1]], 'r')
    plt.plot([x34,x34],[Y[0][y2],Y[0][y1]], 'k')

#Function for plotting velocity as a vector field
def velocity():
    N = 8
    #keeping only every 8th value for better visualizations
    nX = X[::N,::N] ; nY = Y[::N,::N] ; nU = U[::N,::N] ; nV = V[::N,::N] 
    l = np.sqrt(nU**2+nV**2)

    box(169,159)
    box(99,84)
    box(59,49)

    #non-normalized plot
    plt.quiver(nX, nY, nU, nV, l, cmap='jet') 
    plt.plot(xit[0], yit[0], 'k')  #<this line plots the border between water and air
    plt.colorbar()
    plt.title('velocity [mm/s]') 
    plt.xlabel('x-axis [mm]')
    plt.ylabel('y-axis [mm]')
    plt.show()

    box(169,159)
    box(99,84)
    box(59,4)

    #normalized plot
    plt.quiver(nX, nY, nU/l, nV/l, l, cmap='jet')   
    plt.plot(xit[0], yit[0], 'k')  #<this line plots the border between water and air
    plt.colorbar()
    plt.title('velocity normalized') 
    plt.xlabel('x-axis [mm]')
    plt.ylabel('y-axis [mm]')
    plt.show() 



#calculate divergence 
def divergence():
    dudx = np.gradient(U, dx, axis=0) 
    dvdy = np.gradient(V, dy, axis=1) 
    div = dudx + dvdy

    div_max = np.max(div); div_min = np.min(div) 

    plt.figure()

    #overall divergance
    plt.subplot(2,1,1)  
    plt.contourf(X, Y, div, 300, cmap='jet') 
    plt.colorbar()
    plt.title('divergence [s⁻¹]') 
    plt.xlabel('x-axis [mm]')
    plt.ylabel('y-axis [mm]')
    plt.plot(xit[0], yit[0],'k') 

    box(169,159)
    box(99,84)
    box(59,49)

    #Divergance only from 260 s⁻¹ to -260 s⁻¹
    plt.subplot(2,1,2)
    plt.contourf(X, Y, div, levels=np.linspace(-260,260,300), cmap='jet') 
    plt.colorbar()
    plt.title('divergence from -260s⁻¹ to 260s⁻¹ ') 
    plt.xlabel('x-axis [mm]')
    plt.ylabel('y-axis [mm]')
    plt.plot(xit[0], yit[0], 'k') 

    
    box(169,159)
    box(99,84)
    box(59,49)

    plt.show()


#calculate curl
def curl():
    dudy = np.gradient(U, dy, axis=1)
    dvdx = np.gradient(V, dx, axis=0)
    curlz = dvdx - dudy

    plt.contourf(X,Y,curlz, 300, cmap='hsv') 
    plt.colorbar()
    plt.title('Curl [s⁻¹] ') 
    plt.xlabel('x-axis [mm]')
    plt.ylabel('y-axis [mm]')

    plt.plot(xit[0],yit[0],'k') 

    box(169,159)
    box(99,84)
    box(59,49)

    plt.show()


#plot streamlines
def streamline():
    
    l = np.sqrt(u**2+v**2)

    
    plt.streamplot(X[:,0], Y[0,:], u, v, color=l, cmap='plasma') 
    plt.plot(xit[0], yit[0], 'k') 
    plt.colorbar()
    plt.title('Stream') 
    plt.xlabel('x-axis [mm]')
    plt.ylabel('y-axis [mm]')

    plt.show()



#line integral along the lines of the box
def circulation(y1,y2):
    red = np.sum(U[34:69+1, y1]*dx) 
    green = np.sum(V[69, y1:y2+1]*dy)
    blue = np.sum(U[34:69+1, y2]*(-dx)) 
    black =  np.sum(V[34, y1:y2+1]*(-dy))
    return red, green, blue, black 

#fist we test that the data is shaped as expected
sizecheck()

#then we vizualise speed and velocity of the fluids, together whit the stream lines
speed()
velocity()
streamline()

#Next we calculate the divergance and curl over all the field
divergence()
curl()

 
print('Circulation of each box from the top ')

print(np.sum(circulation(159,169)))
print(np.sum(circulation(84,99)))
print(np.sum(circulation(49,59)))

print()
print('Line integral of each side in order of: red, green, blue, black')

print(circulation(159,169))
print(circulation(84,99))
print(circulation(49,59))
print()

#calculation of circulation using Stokes theorem
def stokes(y1,y2):
    dudy = np.gradient(U[34:69+1, y1:y2+1], dy, axis=1)
    dvdx = np.gradient(V[34:69+1, y1:y2+1], dx, axis=0)
    curlz = dvdx - dudy
    res = np.sum(np.sum(curlz, axis=0))*dx*dy
    print(res)

print("Circulation calculated with Stoke's theorem")

stokes(159,169)
stokes(84,99)
stokes(49,59) 
print()

#flux of fluid through the boxes
def flux(y1,y2):
    red = np.sum(-V[34:69+1, y1]*dx) 
    green = np.sum(U[69, y1:y2+1]*dy) 
    blue = np.sum(V[34:69+1, y2]*dx) 
    black =  np.sum(-U[34, y1:y2+1]*dy)
    return red, green, blue, black

print('Flux out of each box from the top')

print(np.sum(flux(159,169)))
print(np.sum(flux(84,99)))
print(np.sum(flux(49,59)))

print()
print('Flux through each line in order of: red, green, blue, black')

print(flux(159,169))
print(flux(84,99))
print(flux(49,59))

"""
Run example:
$ python3 vel_field.py 
Circulation of each box from the top 
2695.5140926958193
-60976.60016211555
9.521016433026077

Line integral of each side in order of: red, green, blue, black
(70100.52387861427, 266.2735761585868, -68332.85609978675, 661.5727377096991)
(198.47559740489237, 300.21661027011714, -61243.464778495945, -231.82759129461652)
(5133.347850903835, 207.91001043390142, -5410.039721925996, 78.30287702128548)

Circulation calculated with Stoke's theorem
2794.1435500042353
-61322.592754967714
2.160413041262494

Flux out of each box from the top
104.8526049082102
-6476.939182097958
-124.5686660449619

Flux through each line in order of: red, green, blue, black
(1556.867943941396, 21664.567474322168, -2059.677184793871, -21056.905628561482)
(-5187.564033067892, 14782.532896182347, -4074.0522144394345, -11997.85583077298)
(-195.57014792583357, 1536.8217966413547, 284.9436464350764, -1750.7639611955597)


"""