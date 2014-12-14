#ifndef SIMPARAMETERS_H
#define SIMPARAMETERS_H

struct SimParameters
{
    SimParameters();

    const static int F_GRAVITY = 1;
    const static int F_STRETCHING = 2;
    const static int F_BENDING = 4;
    const static int F_CONTACT = 8;
    const static int F_DAMPING = 16;

    bool simRunning;
    double timeStep;
    double penaltyStiffness;

    int activeForces;
    double gravityG;
    double stretchingK;
    double bendingK;
    double dampingCoeff;

    bool pinCorner;

    double cor;

};

#endif // SIMPARAMETERS_H

