#ifndef SIGNEDDISTANCEFIELD_H
#define SIGNEDDISTANCEFIELD_H

#include <Eigen/Core>

class RigidBodyTemplate;

class SignedDistanceField
{
public:
    SignedDistanceField(const RigidBodyTemplate &rbt, const char *SDFfile, int defaultGridDim=50);
    ~SignedDistanceField();

    void exportSignedDistanceField(const char *filename);
    bool signedDistanceAndGradient(const Eigen::Vector3d &pos, double &dist, Eigen::Vector3d &Ddist) const;

private:
    void computeSignedDistanceField(const RigidBodyTemplate &rbt, int gridDim);
    bool isInside(const RigidBodyTemplate &rbt, const Eigen::Vector3d &pt);
    double unsignedDistance(const RigidBodyTemplate &rbt, const Eigen::Vector3d &pt);

    double valAt(int i, int j, int k) const;

    int gridDim_;
    double cellsize_;
    double *gridVals_;
};

#endif // SIGNEDDISTANCEFIELD_H
