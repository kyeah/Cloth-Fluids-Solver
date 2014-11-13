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
    string objname("resources/2by4.obj");
    bodyTemplate_ = new RigidBodyTemplate(objname);
    string sdfname("resources/2by4.sdf");
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
    }
    renderLock_.unlock();
}

void Simulation::takeSimulationStep()
{
    time_ += params_.timeStep;

    bodyInstance_->c += params_.timeStep*bodyInstance_->cvel;
    bodyInstance_->theta = VectorMath::axisAngle(VectorMath::rotationMatrix(params_.timeStep*bodyInstance_->w)*VectorMath::rotationMatrix(bodyInstance_->theta));
}

void Simulation::clearScene()
{
    renderLock_.lock();
    {
        delete bodyInstance_;
        Vector3d pos(5, 0, 3);
        Vector3d zero(0,0,0);
        bodyInstance_ = new RigidBodyInstance(*bodyTemplate_, pos, zero, 1.0);
    }
    renderLock_.unlock();
}

void Simulation::accelerateBody(double vx, double vy, double vz, double wx, double wy, double wz)
{
    bodyInstance_->cvel += Vector3d(vx,vy,vz);
    bodyInstance_->w += Vector3d(wx,wy,wz);
}
