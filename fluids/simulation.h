#ifndef SIMULATION_H
#define SIMULATION_H

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <vector>
#include <set>
#include <QMutex>
#include "simparameters.h"
#include <QGLWidget>
#include "mat3d.h"

typedef Eigen::Triplet<double> Tr;

class SimParameters;

class Simulation
{
public:
    Simulation(const SimParameters &params);
    ~Simulation();

    void takeSimulationStep();
    void initializeGL();

    void renderFloor();
    void renderObjects();
    void clearScene();
    void addVelocity(Eigen::Vector2d pos, Eigen::Vector2d vel);

    void stableFluidSolve();
    void diffuse(Mat3D *x, Mat3D *xprev);
    void advect(Mat3D *x, Mat3D *xprev, Mat3D *vx, Mat3D *vy, Mat3D *vz);
    void project(Mat3D *fluidvx, Mat3D *fluidvy, Mat3D *fluidvz);

    void bound_mat(int gridSize, int b, Mat3D *x);

private:
    void loadFloorTexture();

    const SimParameters &params_;
    QMutex renderLock_;

    double time_;
    GLuint floorTex_;

    Mat3D *fluidvx, *fluidvy, *fluidvz,
    *fluidvx_prev, *fluidvy_prev, *fluidvz_prev,
    *fluiddensity, *fluiddensity_prev;
};

#endif // SIMULATION_H

