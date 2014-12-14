#ifndef CLOTH_H
#define CLOTH_H

#include <Eigen/Core>
#include "mesh.h"

class Hinge {
public:
    Hinge(int _ep1, int _ep2, int _v1, int _v2, int _face1, int _face2, double _rls, double _area, double _theta) :
        ep1(_ep1), ep2(_ep2), v1(_v1), v2(_v2), face1(_face1), face2(_face2), restLengthSq(_rls), totalArea(_area), restTheta(_theta) {}

    int ep1, ep2;  // Shared edge point indices
    int v1, v2;    // Unshared vertex indices
    int face1, face2;
    double restLengthSq;
    double totalArea;
    double restTheta;
};

class Cloth {

public:
    Cloth();
    Eigen::VectorXd getVertNormals();
    void render();
    void computeMassMatrices();
    void calculateAmats();
    void fillHinges();
    void resetState();

    Eigen::Vector3d getVert(int i) { return verts_.segment<3>(3*i); }

    Mesh *mesh_;
    Eigen::VectorXd verts_;
    Eigen::VectorXd velocities_;
    Eigen::MatrixXd mass_;
    Eigen::MatrixXd massInv_;
    std::vector<Eigen::MatrixXd> bmats_;
    std::vector<Eigen::MatrixXd> amats_;
    std::vector<Eigen::Matrix2d> gmats_;
    std::vector<Hinge> hinges_;
};

#endif // CLOTH_H
