# start from Example 8.5 program odesim.py
import numpy as np
import matplotlib.pyplot as plt

# constants
#g     = 9.81	# m s^-2
#m     = 1.0	# kg
#rho   = 1.22	# kg m^-3
#C     = 0.47	# unitless
#R     = 0.08    # m
#h     = 0.001   # seconds
#theta = 30.0*(np.pi/180) # radians
#v0    = 100.0	# m s^-1

def universal_artillery(g=9.81, m=1.0, rho=1.22, C=0.47, R=0.08, h=0.001, theta=30.0*(np.pi/180),\
        v0=100):
    const = (rho*C*np.pi*R**2)/(2.0*m)
    # define the equations of motion
    def f(r,const):
        x   = r[0]
        y   = r[1]
        vx  = r[2]
        vy  = r[3]
        fx  = vx
        fy  = vy
        fvx = -const*vx*np.sqrt(vx**2+vy**2)
        fvy = -g-const*vy*np.sqrt(vx**2+vy**2)
        return np.array([fx,fy,fvx,fvy],float)
    
    # try different values of m
    masses = np.linspace(0.1, 10, num=10)
    # list of ranges. Range is last point in ypoints
    yranges = []
    xranges = []
    for m in masses:
        const = (rho*C*np.pi*R**2)/(2.0*m)
        r = np.array([0.0,0.0,v0*np.cos(theta),v0*np.sin(theta)],float)
        xpoints = []
        ypoints = []

        # use fourth-order Runge-Kutta
        while r[1]>=0:
            k1 = h*f(r,const)
            k2 = h*f(r+0.5*k1,const)
            k3 = h*f(r+0.5*k2,const)
            k4 = h*f(r+k3,const)
            r += (k1+2*k2+2*k3+k4)/6
            xpoints.append(r[0])
            ypoints.append(r[1])

        # append range to ranges
        xranges.append(xpoints[-1])
        yranges.append(max(ypoints))

    fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(8, 4))
    ax[0].plot(masses, xranges)
    ax[0].set_xlabel("mass [kg]")
    ax[0].set_ylabel('range [m]')
    ax[1].plot(masses, yranges)
    ax[1].set_xlabel("mass [kg]")
    ax[1].set_ylabel('maximum height [m]')
    
    #plt.legend()
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    # DONE: part a
    # part b
    universal_artillery(g=3.71, rho=0.20)
    # part c: try 10 times the v0

