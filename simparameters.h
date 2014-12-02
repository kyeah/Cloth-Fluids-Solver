#ifndef SIMPARAMETERS_H
#define SIMPARAMETERS_H

struct SimParameters
{
    SimParameters();

    const static int F_GRAVITY = 1;

    bool simRunning;
    double timeStep;

    int activeForces;
    double gravityG;

    int gridSize;
    double viscosity;
    double kDiffusion;

    bool applyFluidDragForce;
};

#endif // SIMPARAMETERS_H

