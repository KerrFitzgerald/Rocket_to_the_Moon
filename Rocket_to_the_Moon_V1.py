# Program to calculate the positions of the Earth & Moon at set time
# intervals as they orbit around a common centre of mass. The motion
# of the rocket is calculated by solving differential equations of
# motion using the Runge-Kutta method. Different rocket trajectories
# are simulated by applying velocity boosts at different time steps.

import math
import numpy
from matplotlib import animation
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# Earth mass
e_mass = 5.9742e24
# Moon mass
m_mass = 7.3500e22
# Total mass
t_mass = m_mass + e_mass
# Earth Moon seperation
em_sep = 3.8440e8
# Gravitational constant
G = 6.6726e-11
# comd_r from the centre of mass to the Earth
comd_e = em_sep*(m_mass/t_mass)
# comd_r from the centre of mass to the Moon
comd_m = em_sep*(e_mass/t_mass)
# Angular Velocity of both the Moon & Earth
w = (math.sqrt((G*e_mass)/(comd_m*((em_sep)**2))))
# L1 calculated position using timestep of 10 seconds
L1 = 58022627.9
# L2 calculated position using timestep of 10 seconds
L2 = 64521222.9550
# Rocket orbital height above Earth
orbh_r = 6800000


def EMR_Animation(xE, yE, xM, yM, xR, yR):
    '''
    Generate animation of the positions of the Earth, Moon & Rocket.
    '''

    xRa, yRa, xEa, yEa, xMa, yMa = [], [], [], [], [], []
    for i in range(len(xR)):
        if i % 500 == 0:
            xRa.append((xR[i]))
            yRa.append((yR[i]))
            xEa.append((xE[i]))
            yEa.append((yE[i]))
            xMa.append((xM[i]))
            yMa.append((yM[i]))
    fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(15, 15))
    axes.set_ylim(-6E8, 6E8)
    axes.set_xlim(-6E8, 6E8)
    plt.style.use("ggplot")
    xRb, yRb, xEb, yEb, xMb, yMb = [], [], [], [], [], []

    def animate(i):
        xRb.append((xRa[i]))
        yRb.append((yRa[i]))
        xEb.append((xEa[i]))
        yEb.append((yEa[i]))
        xMb.append((xMa[i]))
        yMb.append((yMa[i]))
        axes.plot(xRb, yRb, 'g--')
        axes.plot(xEb, yEb, 'b--')
        axes.plot(xMb, yMb, 'r--')
    anim = FuncAnimation(fig, animate, interval=1)
    plt.show()

    return None


def EMR_Plotter(xE, yE, xM, yM, xR, yR):
    '''
    Plots the positions of the Earth, Moon & Rocket.
    '''

    fig = plt.figure()
    plt.ylabel('y/10^8m')
    plt.xlabel('x/10^8m')
    plt.plot(xR/1e8, yR/1e8, 'g-', label='Rocket')
    plt.plot(xE/1e8, yE/1e8, 'b-', label='Earth')
    plt.plot(xM/1e8, yM/1e8, 'r-', label='Moon')
    plt.show()

    return None


def Orbital_Positions_em(comd_e, comd_m, w, t_step, fract):
    '''
    Calculates the positions of the Earth & Moon over one orbital
    period of the system at sucessive time steps. The positions are
    outputted as arrays.
    '''

    # Time period of the Earth & Moon Orbit
    T = (2*math.pi)/w
    T_array = numpy.arange(0, T*fract, t_step)
    # x & y positions of the Earth
    xE = comd_e*numpy.cos(w*T_array + math.pi)
    yE = comd_e*numpy.sin(w*T_array + math.pi)
    # x & y positions of the Moon
    xM = comd_m*numpy.cos(w*T_array)
    yM = comd_m*numpy.sin(w*T_array)

    return xE, yE, xM, yM


def Initialize_Arrays(w, t_step, fract):
    '''
    Creates empty arrays for the rockets position, velocity,
    acceleration & moon separation.
    '''

    T = (2*math.pi)/w
    T_array = numpy.arange(0, T*fract, t_step)
    xR = numpy.zeros((len(T_array), 1), dtype=float)
    yR = numpy.zeros((len(T_array), 1), dtype=float)
    vxR = numpy.zeros((len(T_array), 1), dtype=float)
    vyR = numpy.zeros((len(T_array), 1), dtype=float)
    axR = numpy.zeros((len(T_array), 1), dtype=float)
    ayR = numpy.zeros((len(T_array), 1), dtype=float)
    # Rocket-Moon seperation distance
    rm_sep = numpy.zeros((len(T_array), 1), dtype=float)

    return xR, yR, vxR, vyR, axR, ayR, rm_sep


def Acceleration(xR, yR, xE, yE, xM, yM):
    '''
    Calculates the Acceleration experienced by the rocket due to the
    gravitational attraction of the Earth & Moon.
    '''

    # Rocket-Earth distance
    dE = math.sqrt((xR-xE)**2 + (yR-yE)**2)
    # Rocket-Moon distance
    dM = math.sqrt((xR-xM)**2 + (yR-yM)**2)
    # Calculates rocket x & y accelerations from gravitational forces
    axR = (-G*e_mass*(xR-xE)/(dE**3)) - (G*m_mass*(xR-xM)/(dM**3))
    ayR = (-G*e_mass*(yR-yE)/(dE**3)) - (G*m_mass*(yR-yM)/(dM**3))

    return axR, ayR


def RungeZ1(X, Y, vX, vY, aX, aY, t_step):
    '''
    Calculates Z1 position & velocity components for use in the
    Runge-kutta Method.
    '''

    Z1x = X + (t_step/2.0)*vX
    Z1y = Y + (t_step/2.0)*vY
    vZ1x = vX + (t_step/2.0)*aX
    vZ1y = vY + (t_step/2.0)*aY

    return Z1x, Z1y, vZ1x, vZ1y


def RungeZ2(X, Y, vX, vY, vZ1X, vZ1Y, aZ1X, aZ1Y, t_step):
    '''
    Calculates Z2 position & velocity components for use in the
    Runge-kutta Method.
    '''

    Z2x = X + (t_step/2.0)*vZ1X
    Z2y = Y + (t_step/2.0)*vZ1Y
    vZ2x = vX + (t_step/2.0)*aZ1X
    vZ2y = vY + (t_step/2.0)*aZ1Y
    return Z2x, Z2y, vZ2x, vZ2y


def RungeZ3(X, Y, vX, vY, vZ2X, vZ2Y, aZ2X, aZ2Y, t_step):
    '''
    Calculates Z3 position & velocity components for use in the
    Runge-kutta Method.
    '''

    Z3x = X + (t_step/1.0)*vZ2X
    Z3y = Y + (t_step/1.0)*vZ2Y
    vZ3x = vX + (t_step/1.0)*aZ2X
    vZ3y = vY + (t_step/1.0)*aZ2Y
    return Z3x, Z3y, vZ3x, vZ3y


def Rocket_Motion(traject, w, t_step, comd_e, comd_m, fract):
    '''
    Function that calculates the rocket trajectory allowing for
    different starting locations & boost sequences.
    '''

    theta = math.pi
    T = (2*math.pi)/w
    T_array = numpy.arange(0, T*fract, t_step)
    x, y, vx, vy, ax, ay, comd_r = Initialize_Arrays(w, t_step, fract)
    rm_sep = numpy.zeros((len(x), 1), dtype=float)
    xE, yE, xM, yM = Orbital_Positions_em(comd_e, comd_m, w, t_step, fract)
    if traject[0] == 'Earth_Orbit':
        print('The rocket is starting in ', traject[0])
        wR = math.sqrt((G*e_mass)/(orbh_r**3))
        # Calculated the angle at which the boost is applied
        x[0] = -comd_e + (orbh_r*math.cos(theta))
        y[0] = orbh_r*math.sin(theta)
        vx[0] = -wR*orbh_r*math.sin(theta) + w*yE[0]
        vy[0] = wR*orbh_r*math.cos(theta) + w*xE[0]
    if traject[0] == 'L1':
        x[0] = comd_m - L1
        vy[0] = w*x[0]
        print('Rocket starts at Lagrangian Point 1')
    if traject[0] == 'L2':
        x[0] = comd_m + L2
        vy[0] = w*x[0]
        print('Rocket starts at Lagrangian Point 2')
    ax[0], ay[0] = Acceleration(x[0], y[0], xE[0], yE[0], xM[0], yM[0])
    tracker = 0
    # Update the rockets position using the Runge-Kutta method
    for i in range(1, len(xM)):
        Z1x, Z1y, vZ1x, vZ1y = RungeZ1(
                                       x[i-1], y[i-1],
                                       vx[i-1], vy[i-1],
                                       ax[i-1], ay[i-1],
                                       t_step
                                       )
        aZ1x, aZ1y = Acceleration(
                                  Z1x, Z1y,
                                  xE[i-1], yE[i-1],
                                  xM[i-1], yM[i-1]
                                  )
        Z2x, Z2y, vZ2x, vZ2y = RungeZ2(x[i-1], y[i-1], vx[i-1], vy[i-1],
                                       vZ1x, vZ1y,
                                       aZ1x, aZ1y,
                                       t_step
                                       )
        aZ2x, aZ2y = Acceleration(
                                  Z2x, Z2y,
                                  xE[i-1], yE[i-1],
                                  xM[i-1], yM[i-1]
                                  )
        Z3x, Z3y, vZ3x, vZ3y = RungeZ3(
                                       x[i-1], y[i-1],
                                       vx[i-1], vy[i-1],
                                       vZ2x, vZ2y,
                                       aZ2x, aZ2y,
                                       t_step)
        aZ3x, aZ3y = Acceleration(Z3x, Z3y,
                                  xE[i-1], yE[i-1],
                                  xM[i-1], yM[i-1]
                                  )
        # Calculates new rocket x & y positions
        x[i] = x[i-1] + ((t_step/6.0)*(vx[i-1] + 2*vZ1x + 2*vZ2x + vZ3x))
        y[i] = y[i-1] + ((t_step/6.0)*(vy[i-1] + 2*vZ1y + 2*vZ2y + vZ3y))
        # Calculates new rocket x & y velocities
        vx[i] = vx[i-1] + ((t_step/6.0)*(ax[i-1] + 2*aZ1x + 2*aZ2x + aZ3x))
        vy[i] = vy[i-1] + ((t_step/6.0)*(ay[i-1] + 2*aZ1y + 2*aZ2y + aZ3y))
        # Calculates new rocket x & y accelerations
        ax[i], ay[i] = Acceleration(x[i], y[i], xE[i], yE[i], xM[i], yM[i])
        comd_r[i] = math.sqrt((x[i])**2 + (y[i])**2)
        rm_sep[i] = math.sqrt((x[i]-xM[i])**2 + (y[i]-yM[i])**2)
        # Check if a boost should be applied at the current timestep
        for j in range(0, len(traject[1])):
            if i == traject[1][j]:
                angle = math.atan(vy[i]/vx[i])
                ang_deg = angle*57.2958
                boost = traject[2][j]
                print('Applying boost', j+1, 'at step:', i)
                print('The boost magnitude is', boost)
                print('Angle of rocket with horizontal is:', ang_deg)
                # Applies the velocity boost to x & y velocities
                vx[i] = vx[i] + boost*math.cos(angle)
                vy[i] = vy[i] + boost*math.sin(angle)

    return x, y, xE, yE, xM, yM, comd_r, rm_sep, T_array


if __name__ == '__main__':
    t_step = 10
    fract = 0.5
    # Langragian Point 1 (10s)
    traject1 = ['L1', [0], [0]]
    # Langragian Point 2 (10s)
    traject2 = ['L2', [0], [0]]
    # Free Return Trajectory (10s)
    traject3 = ['Earth_Orbit', [85], [3086]]
    # Moon Orbit (10s)
    traject4 = ['Earth_Orbit', [85, 31445], [3086.51, -600]]

    xR, yR, xE, yE, xM, yM, comd_r, rm_sep, T_array = \
        Rocket_Motion(traject2, w, t_step, comd_e, comd_m, fract)
    EMR_Animation(xE, yE, xM, yM, xR, yR)
