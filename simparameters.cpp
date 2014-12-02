#include "simparameters.h"

SimParameters::SimParameters()
{
    simRunning = false;
    timeStep = 0.0001;

    activeForces = F_GRAVITY | F_STRETCHING | F_DAMPING | F_BENDING | F_CONTACT;
    gravityG = -3.4;
    penaltyStiffness = 100000;
    pinCorner = false;

    stretchingK = 10000;
    bendingK = 1e-2;
    dampingCoeff = 0.1;
    cor = 0.3;

    gridSize = 30;
    viscosity = .00028;
}
