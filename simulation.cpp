#include "simulation.h"
#include <QGLWidget>
#include "simparameters.h"
#include <iostream>
#include <Eigen/Geometry>
#include <QDebug>
#include "SOIL.h"
#include "rigidbodytemplate.h"
#include "rigidbodyinstance.h"
#include "vectormath.h"
#include <Eigen/Dense>
#include "mesh.h"
#include "signeddistancefield.h"

const double PI = 3.1415926535898;

using namespace Eigen;
using namespace std;

Simulation::Simulation(const SimParameters &params) : params_(params), time_(0), floorTex_(0)
{
    loadRigidBodies();
    cloth_ = new Cloth();
    bodyInstance_ = NULL;
    clearScene();
}

Simulation::~Simulation()
{
    delete bodyInstance_;
    delete bodyTemplate_;
}

void Simulation::initializeGL()
{
    loadFloorTexture();
}

void Simulation::loadRigidBodies()
{
    string objname("resources/sphere.obj");
    bodyTemplate_ = new RigidBodyTemplate(objname);
    string sdfname("resources/sphere.sdf");
    bodyTemplate_->computeSDF(sdfname.c_str());
}

void Simulation::loadFloorTexture()
{
    floorTex_ = SOIL_load_OGL_texture("resources/grid.jpg", SOIL_LOAD_AUTO, SOIL_CREATE_NEW_ID, SOIL_FLAG_INVERT_Y |  SOIL_FLAG_NTSC_SAFE_RGB | SOIL_FLAG_COMPRESS_TO_DXT | SOIL_FLAG_MIPMAPS);
    if(floorTex_ != 0)
    {
        glBindTexture(GL_TEXTURE_2D, floorTex_);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
    }
}


void Simulation::renderFloor()
{
    renderLock_.lock();

    glColor4f(1.0, 1.0, 1.0, 1.0);

    if(floorTex_)
    {
        glBindTexture(GL_TEXTURE_2D, floorTex_);
        glEnable(GL_TEXTURE_2D);
    }
    else
    {
        glColor3f(0.5, 0.5, 0.5);
        glDisable(GL_TEXTURE_2D);
    }

    double texsize = 5.0;
    double gridsize = 1000.0;

    double texmax = gridsize/texsize;

    Vector3d tangent1(1.0,0,0);
    Vector3d tangent2(0,-1.0,0);

    Vector3d corner;

    glBegin(GL_QUADS);
    {
        glTexCoord2f(texmax, texmax);
        glNormal3f(0, 0, 1.0);
        corner = -gridsize*tangent1 + gridsize*tangent2;
        glVertex3f(corner[0], corner[1], corner[2]);

        glTexCoord2f(texmax, -texmax);
        glNormal3f(0, 0, 1.0);
        corner = -gridsize*tangent1 - gridsize*tangent2;
        glVertex3f(corner[0], corner[1], corner[2]);

        glTexCoord2f(-texmax, -texmax);
        glNormal3f(0, 0, 1.0);
        corner = gridsize*tangent1 - gridsize*tangent2;
        glVertex3f(corner[0], corner[1], corner[2]);

        glTexCoord2f(-texmax, texmax);
        glNormal3f(0, 0, 1.0);
        corner = gridsize*tangent1 + gridsize*tangent2;
        glVertex3f(corner[0], corner[1], corner[2]);
    }
    glDisable(GL_TEXTURE_2D);
    glEnd();
    renderLock_.unlock();
}

void Simulation::renderObjects()
{
    renderLock_.lock();
    {
        bodyInstance_->render();
        cloth_->render();
    }
    renderLock_.unlock();
}

void Simulation::takeSimulationStep()
{
    time_ += params_.timeStep;

    bodyInstance_->c += params_.timeStep*bodyInstance_->cvel;
    bodyInstance_->theta = VectorMath::axisAngle(VectorMath::rotationMatrix(params_.timeStep*bodyInstance_->w)*VectorMath::rotationMatrix(bodyInstance_->theta));

    cloth_->verts_ += params_.timeStep * cloth_->velocities_;

    VectorXd gravForce = computeGravForce();
    VectorXd clothForce = computeClothForce();
    VectorXd dampingForce = -params_.dampingCoeff * cloth_->velocities_;

    if (!(params_.activeForces & SimParameters::F_DAMPING)) {
        dampingForce.setZero();
    }

    cloth_->velocities_ += params_.timeStep * (gravForce + dampingForce + cloth_->massInv_ * clothForce);

    if (params_.pinCorner) {
        cloth_->velocities_.segment<3>(0) = Vector3d(0,0,0);
    }
}

void Simulation::clearScene()
{
    renderLock_.lock();
    {
        delete bodyInstance_;
        Vector3d pos(5, 0, 3);
        Vector3d zero(0,0,0);
        bodyInstance_ = new RigidBodyInstance(*bodyTemplate_, pos, zero, 1.0);
        cloth_->resetState();
    }
    renderLock_.unlock();
}

void Simulation::accelerateBody(double vx, double vy, double vz, double wx, double wy, double wz)
{
    bodyInstance_->cvel += Vector3d(vx,vy,vz);
    bodyInstance_->w += Vector3d(wx,wy,wz);
}

VectorXd Simulation::computeGravForce() {
    VectorXd vertPos = cloth_->verts_;
    VectorXd force(vertPos.rows());
    force.setZero();

    if (params_.activeForces & SimParameters::F_GRAVITY) {
        for (int i = 2; i < force.rows(); i+=3) {
            force[i] += params_.gravityG;
        }
    }
    return force;
}

VectorXd Simulation::computeClothForce() {
    VectorXd vertPos = cloth_->verts_;
    VectorXd force(vertPos.rows());
    force.setZero();

    // Stretching Force
    // DqV = 2Ev.transpose() * DqEv
    // DqEv = DdgEvDqd
    // DdgEv = 0.5[A]Dqdgv
    // Dqdgv = Dq(gv - gvtemplate) = Dqgv = DevgvDqev = [C]dqev
    // dqev = [D]dq
    // DqV = 2Ev.transpose() * 0.5[A][C][D]dqdqd
    // F = -(DqV).transpose() = -[D.t][C.t][A.t]Ev where Ev = 9x1

    Matrix3d I = Matrix3d::Identity();

    if (params_.activeForces & SimParameters::F_STRETCHING) {
        MatrixXd D(6,9);
        D.setZero();
        D.block<3,3>(0,0) = -I;
        D.block<3,3>(0,3) = I;
        D.block<3,3>(3,0) = -I;
        D.block<3,3>(3,6) = I;

        for (int i = 0; i < cloth_->mesh_->getNumFaces(); i++) {
            Vector3i fverts = cloth_->mesh_->getFace(i);
            Vector3d pts[3];
            for(int j=0; j<3; j++)
                pts[j] = cloth_->verts_.segment<3>(3*fverts[j]);

            Vector3d e1 = pts[1] - pts[0];
            Vector3d e2 = pts[2] - pts[0];
            MatrixXd bmat = cloth_->bmats_[i];

            Matrix2d g;
            g(0,0) = e1.dot(e1);
            g(0,1) = e1.dot(e2);
            g(1,0) = g(0,1);
            g(1,1) = e2.dot(e2);

            Matrix2d dg = g - cloth_->gmats_[i];

            MatrixXd C(4,6);
            C.setZero();
            C.block<1,3>(0,0) = 2*e1.transpose();
            C.block<1,3>(1,0) = e2.transpose();
            C.block<1,3>(1,3) = e1.transpose();
            C.block<1,3>(2,0) = e2.transpose();
            C.block<1,3>(2,3) = e1.transpose();
            C.block<1,3>(3,3) = 2*e2.transpose();

            Matrix3d epsilon = bmat.transpose() * dg * bmat;
            VectorXd epsilonv(9,1);
            epsilonv.segment<3>(0) = epsilon.block<1,3>(0,0).transpose();
            epsilonv.segment<3>(3) = epsilon.block<1,3>(1,0).transpose();
            epsilonv.segment<3>(6) = epsilon.block<1,3>(2,0).transpose();
            epsilonv /= 2.0;

            MatrixXd A = cloth_->amats_[i];

            VectorXd df = -D.transpose() * C.transpose() * A.transpose() * epsilonv;
            df *= cloth_->mesh_->getFaceArea(i) * params_.stretchingK;

            for (int j=0; j < 3; j++) {
                force.segment<3>(3*fverts[j]) += df.segment<3>(3*j);
            }
        }
    }

    if (params_.activeForces & SimParameters::F_BENDING) {
        for (int i = 0; i < cloth_->hinges_.size(); i++) {
            Hinge hinge = cloth_->hinges_[i];
            Vector3d pi = cloth_->getVert(hinge.ep1);
            Vector3d pj = cloth_->getVert(hinge.ep2);
            Vector3d pk = cloth_->getVert(hinge.v1);
            Vector3d pl = cloth_->getVert(hinge.v2);

            Vector3d n0 = (pj-pi).cross(pk-pi);
            Vector3d n1 = (pl-pi).cross(pj-pi);

            Vector3d n0xn1 = n0.cross(n1);
            double theta = 2*atan2(n0xn1.norm(), n0.norm()*n1.norm() + n0.dot(n1)) - hinge.restTheta;

            if (theta == 0) {
                continue;
            }

            double rest_coeff = params_.bendingK * hinge.restLengthSq * 2 * theta / hinge.totalArea;

            Vector3d Dn1_coeff = (n0xn1/n0xn1.norm()).cross(n1/n1.squaredNorm());
            Vector3d Dn0_coeff = (n0xn1/n0xn1.norm()).cross(n0/n0.squaredNorm());

            Matrix3d Dn0_i = VectorMath::crossProductMatrix(pk-pj);
            Matrix3d Dn0_j = VectorMath::crossProductMatrix(pi-pk);
            Matrix3d Dn0_k = VectorMath::crossProductMatrix(pj-pi);

            Matrix3d Dn1_i = VectorMath::crossProductMatrix(pj-pl);
            Matrix3d Dn1_j = VectorMath::crossProductMatrix(pl-pi);
            Matrix3d Dn1_l = VectorMath::crossProductMatrix(pi-pj);

            Vector3d di = -(Dn1_coeff.transpose()*Dn1_i - Dn0_coeff.transpose()*Dn0_i).transpose();
            Vector3d dj = -(Dn1_coeff.transpose()*Dn1_j - Dn0_coeff.transpose()*Dn0_j).transpose();
            Vector3d dk = -(-Dn0_coeff.transpose()*Dn0_k).transpose();
            Vector3d dl = -(Dn1_coeff.transpose()*Dn1_l).transpose();

            force.segment<3>(3*hinge.ep1) += rest_coeff * di;
            force.segment<3>(3*hinge.ep2) += rest_coeff * dj;
            force.segment<3>(3*hinge.v1) += rest_coeff * dk;
            force.segment<3>(3*hinge.v2) += rest_coeff * dl;
        }
    }

    if (params_.activeForces & SimParameters::F_CONTACT) {
        for (int i=0; i < cloth_->mesh_->getNumVerts(); i++) {
            Vector3d pt = cloth_->getVert(i);

            double dist;
            Vector3d Ddist;
            Matrix3d bodyrot = VectorMath::rotationMatrix(bodyInstance_->theta);
            Vector3d testpt = bodyrot.transpose()*(pt - bodyInstance_->c);

            if(bodyInstance_->getTemplate().getSDF()->signedDistanceAndGradient(testpt, dist, Ddist) && dist < 0)
            {
                Vector3d vel1 = cloth_->velocities_.segment<3>(3*i);
                Vector3d vel2 = bodyInstance_->cvel +bodyrot*(bodyInstance_->w.cross(testpt));
                Vector3d relvel = (vel1-vel2);
                double stiffness = params_.penaltyStiffness * (relvel.dot(bodyrot*Ddist) > 0 ? params_.cor : 1.0);
                force.segment<3>(3*i) += -stiffness*dist*bodyrot*Ddist;
            }
        }
    }
    return force;
}
