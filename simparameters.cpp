#include "simparameters.h"

SimParameters::SimParameters()
{
    simRunning = false;
    timeStep = 0.0001;

    activeForces = F_GRAVITY | F_STRETCHING | F_DAMPING | F_BENDING | F_CONTACT;
    gravityG = -9.8;
    penaltyStiffness = 100000;
    pinCorner = false;

    stretchingK = 10000;
    bendingK = 1e-2;
    dampingCoeff = 0.1;
    cor = 0.3;
}
