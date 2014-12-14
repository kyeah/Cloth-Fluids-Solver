#include "signeddistancefield.h"
#include "rigidbodytemplate.h"
#include "mesh.h"
#include "collisions/Distance.h"
#include "collisions/CTCD.h"
#include "vectormath.h"
#include <fstream>
#include <iostream>

using namespace Eigen;
using namespace std;

SignedDistanceField::SignedDistanceField(const RigidBodyTemplate &rbt, const char *SDFfile, int defaultGridDim)
{
    ifstream ifs(SDFfile);
    if(!ifs)
    {
        computeSignedDistanceField(rbt, defaultGridDim);
        exportSignedDistanceField(SDFfile);
        return;
    }

    ifs >> gridDim_;

    cellsize_ = 2.0/double(gridDim_-1);
    gridVals_ = new double[gridDim_*gridDim_*gridDim_];
    for(int i=0; i<gridDim_; i++)
    {
        for(int j=0; j<gridDim_; j++)
        {
            for(int k=0; k<gridDim_; k++)
            {
                ifs >> gridVals_[gridDim_*gridDim_*i + gridDim_*j + k];
            }
        }
    }
}

void SignedDistanceField::computeSignedDistanceField(const RigidBodyTemplate &rbt, int gridDim)
{
    gridDim_ = gridDim;
    cellsize_ = 2.0/double(gridDim-1);

    gridVals_ = new double[gridDim*gridDim*gridDim];

    for(int i=0; i<gridDim; i++)
    {
        for(int j=0; j<gridDim; j++)
        {
            for(int k=0; k<gridDim; k++)
            {
                Vector3d pos(-1.0 + i*cellsize_, -1.0+j*cellsize_, -1.0+k*cellsize_);
                double sign = isInside(rbt, pos) ? -1.0 : 1.0;
                double dist = sign*unsignedDistance(rbt, pos);
                gridVals_[gridDim*gridDim*i + gridDim*j + k] = dist;
            }
        }
    }
}

SignedDistanceField::~SignedDistanceField()
{
    delete[] gridVals_;
}

bool SignedDistanceField::isInside(const RigidBodyTemplate &rbt, const Vector3d &pt)
{
    int votesinside = 0;
    int votesoutside = 0;

    for(int trial=0; trial<10; trial++)
    {
        Vector3d randray;
        for(int i=0; i<3; i++)
            randray[i] = 2.0*VectorMath::randomUnitIntervalReal()-1.0;
        randray.normalize();

        double r = 1.0;

        const Mesh &m = rbt.getMesh();

        int numfaces = m.getNumFaces();
        int count = 0;
        for(int i=0; i<numfaces; i++)
        {
            const Vector3i face = m.getFace(i);
            Vector3d endpos = pt+4.0*r*randray;
            double t;
            if(CTCD::vertexFaceCTCD(pt, m.getVert(face[0]), m.getVert(face[1]), m.getVert(face[2]), endpos, m.getVert(face[0]), m.getVert(face[1]), m.getVert(face[2]),1e-8,t))
            {
                count++;
            }
        }
        if(count % 2 == 0)
            votesoutside++;
        else
            votesinside++;
    }
    return votesinside > votesoutside;
}

double SignedDistanceField::unsignedDistance(const RigidBodyTemplate &rbt, const Eigen::Vector3d &pt)
{
    const Mesh &m = rbt.getMesh();
    int numfaces = m.getNumFaces();
    double mindist = std::numeric_limits<double>::max();
    for(int i=0; i<numfaces; i++)
    {
        const Vector3i face = m.getFace(i);
        Vector3d pt0 = m.getVert(face[0]);
        Vector3d pt1 = m.getVert(face[1]);
        Vector3d pt2 = m.getVert(face[2]);

        double d;
        double dist = Distance::vertexFaceDistance(pt, pt0, pt1, pt2, d, d, d).norm();
        mindist = std::min(mindist, dist);
    }
    return mindist;
}

void SignedDistanceField::exportSignedDistanceField(const char *filename)
{
    ofstream ofs(filename);
    ofs << gridDim_;
    ofs << endl;
    for(int i=0; i<gridDim_; i++)
    {
        for(int j=0; j<gridDim_; j++)
        {
            for(int k=0; k<gridDim_; k++)
            {
                ofs << gridVals_[gridDim_*gridDim_*i + gridDim_*j + k] << " ";
            }
        }
    }
}

double SignedDistanceField::valAt(int i, int j, int k) const
{
    return gridVals_[gridDim_*gridDim_*i + gridDim_*j + k];
}

bool SignedDistanceField::signedDistanceAndGradient(const Vector3d &pos, double &dist, Vector3d &Ddist) const
{
    int i = floor((1.0 + pos[0])/cellsize_);
    int j = floor((1.0 + pos[1])/cellsize_);
    int k = floor((1.0 + pos[2])/cellsize_);

    if(i < 0 || i >= gridDim_-1
            || j < 0 || j >= gridDim_-1
            || k < 0 || k >= gridDim_-1)
    {
        dist = std::numeric_limits<double>::infinity();
        Ddist.setZero();
        return false;
    }

    double r = (1.0+pos[0] - i*cellsize_)/cellsize_;
    double s = (1.0+pos[1] - j*cellsize_)/cellsize_;
    double t = (1.0+pos[2] - k*cellsize_)/cellsize_;

    dist=0;
    Ddist.setZero();

    for(int a=0; a<2; a++)
    {
        double term1 = (a==0 ? 1.0-r : r);
        double dterm1 = (a==0 ? -1.0 : 1.0);
        for(int b=0; b<2; b++)
        {
            double term2 = (b==0 ? 1.0-s : s);
            double dterm2 = (b==0 ? -1.0 : 1.0);
            for(int c=0; c<2; c++)
            {
                double term3 = (c==0 ? 1.0-t : t);
                double dterm3 = (c==0 ? -1.0 : 1.0);
                dist += term1*term2*term3*valAt(i+a, j+b, k+c);
                Ddist[0] += dterm1*term2*term3*valAt(i+a, j+b, k+c)/cellsize_;
                Ddist[1] += term1*dterm2*term3*valAt(i+a, j+b, k+c)/cellsize_;
                Ddist[2] += term1*term2*dterm3*valAt(i+a, j+b, k+c)/cellsize_;
            }
        }
    }
    return true;
}
