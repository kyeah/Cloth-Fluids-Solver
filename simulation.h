#ifndef SIMULATION_H
#define SIMULATION_H

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <vector>
#include <set>
#include <QMutex>
#include "simparameters.h"
#include <QGLWidget>
#include "cloth.h"
#include "mat3d.h"

class RigidBodyTemplate;
class RigidBodyInstance;

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
    void accelerateBody(double vx, double vy, double vz, double wx, double wy, double wz);
    Eigen::VectorXd computeGravForce();
    Eigen::VectorXd computeClothForce();
    void stableFluidSolve();
    void set_bnd ( int N, int b, Mat2D *x );

private:
    void loadFloorTexture();
    void loadRigidBodies();

    const SimParameters &params_;
    QMutex renderLock_;

    double time_;
    GLuint floorTex_;

    RigidBodyTemplate * bodyTemplate_;
    RigidBodyInstance * bodyInstance_;
    Cloth *cloth_;
    Mat2D *fluidvx, *fluidvy, *fluidvx_prev, *fluidvy_prev, *fluiddensity, *fluiddensity_prev;
};

#endif // SIMULATION_H

