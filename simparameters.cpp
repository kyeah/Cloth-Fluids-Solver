#include "simparameters.h"

SimParameters::SimParameters()
{
    simRunning = false;
    timeStep = 0.0001;

    activeForces = F_GRAVITY;
    gravityG = -3.4;

    gridSize = 30;
    viscosity = .00028;
    kDiffusion = .0001;
    applyFluidDragForce = false;
    dragForceMag = 10000;
}
