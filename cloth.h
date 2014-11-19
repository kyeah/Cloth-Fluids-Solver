#ifndef CLOTH_H
#define CLOTH_H

#include <Eigen/Core>
#include "mesh.h"

class Cloth {

public:
    Cloth();
    Eigen::VectorXd getVertNormals();
    void render();
    void computeMassMatrices();
    void resetState();

    Mesh *mesh_;
    Eigen::VectorXd verts_;
    Eigen::VectorXd velocities_;
    Eigen::MatrixXd mass_;
    Eigen::MatrixXd massInv_;
};

#endif // CLOTH_H
