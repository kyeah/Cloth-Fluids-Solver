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
#include <unsupported/Eigen/FFT>

#define SWAP_MAT(A,B) {Mat3D *tmp=A;A=B;B=tmp;}

const double PI = 3.1415926535898;

using namespace Eigen;
using namespace std;

Simulation::Simulation(const SimParameters &params) : params_(params), time_(0), floorTex_(0)
{
    loadRigidBodies();
    cloth_ = new Cloth();
    bodyInstance_ = NULL;
    fluidvx = fluidvy = fluidvz = fluidvx_prev = fluidvy_prev = fluidvz_prev = fluiddensity = fluiddensity_prev = NULL;
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
        cloth_->render();
        glPointSize(10);
        glBegin(GL_POINTS);
        float maxVal = 0;
        for (int i = 0; i < params_.gridSize; i++) {
            for (int j = 0; j < params_.gridSize; j++) {
                for (int k = 0; k < params_.gridSize; k++) {
                    float val = fluiddensity->valAt(i,j,k);
                    if (val > maxVal) {
                        maxVal = val;
                    }
                }
            }
        }

        if (maxVal == 0) maxVal = 1;

        for (int i = 1; i <= params_.gridSize; i++) {
            for (int j = 1; j <= params_.gridSize; j++) {
                float val = (fluiddensity->valAt(i,j) / maxVal) * 255;
                while (val > 255) {
                    val -= 255;
                }
                glColor4f(val, 0, val, val + 0.1);
                glVertex3f(.25*(params_.gridSize - i), .25*k, .25*(j + 5));
            }
            }
        }
        glEnd();
    }
    renderLock_.unlock();
}

void Simulation::takeSimulationStep()
{
    time_ += params_.timeStep;

    /*
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
    }*/


    for (int i = 1; i <= params_.gridSize; i++) {
        for (int j = 1; j <= params_.gridSize; j++) {
            fluidvy_prev->valAt(i, j) -= .98;
            fluidvy->valAt(i, j) -= .98;
            //fluiddensity_prev->valAt(i, j) = fluiddensity->valAt(i, j);
        }
    }
    //fluiddensity->valAt(params_.gridSize / 2, params_.gridSize / 2) += 1;


    //for (int i = 0; i < params_.gridSize; i++) {
    //    fluidvy_prev->valAt(0, i) = -2*fluidvy_prev->valAt(1, i);
   // }

    // Source
    stableFluidSolve();
}

void Simulation::addVelocity(Eigen::Vector2d pos, Eigen::Vector2d vel) {
    int px = floor((pos[0] + 1)*params_.gridSize / 2);
    int py = floor((pos[1] + 1)*params_.gridSize / 2);
    vel = 10 * (vel / vel.norm());

    if (px < 0 || px >= params_.gridSize - 1 || py < 0 || py >= params_.gridSize - 1) {
        return;
    }

    fluiddensity->valAt(px, py, params_.gridSize/2) += .001;

    //fluiddensity_prev->valAt(px, py) += .001;

    //fluidvx->valAt(px, py) += vel[0];
    //fluidvy->valAt(px, py) += vel[1];
}

void Simulation::clearScene()
{
    renderLock_.lock();
    {
        delete bodyInstance_;
        delete fluidvx;
        delete fluidvy;
        delete fluidvz;
        delete fluidvx_prev;
        delete fluidvy_prev;
        delete fluidvz_prev;
        delete fluiddensity;
        delete fluiddensity_prev;
        Vector3d pos(5, 0, 3);
        Vector3d zero(0,0,0);
        bodyInstance_ = new RigidBodyInstance(*bodyTemplate_, pos, zero, 1.0);
        cloth_->resetState();
        fluidvx = new Mat3D(params_.gridSize+2, params_.gridSize+2, params_.gridSize+2);
        fluidvy = new Mat3D(params_.gridSize+2, params_.gridSize+2, params_.gridSize+2);
        fluidvz = new Mat3D(params_.gridSize+2, params_.gridSize+2, params_.gridSize+2);
        fluidvx_prev = new Mat3D(params_.gridSize+2, params_.gridSize+2, params_.gridSize+2);
        fluidvy_prev = new Mat3D(params_.gridSize+2, params_.gridSize+2, params_.gridSize+2);
        fluidvz_prev = new Mat3D(params_.gridSize+2, params_.gridSize+2, params_.gridSize+2);
        fluiddensity = new Mat3D(params_.gridSize+2, params_.gridSize+2, params_.gridSize+2);
        fluiddensity_prev = new Mat3D(params_.gridSize+2, params_.gridSize+2, params_.gridSize+2);
        //fluiddensity->valAt(params_.gridSize/2, params_.gridSize/2, params_.gridSize/2) = 1;
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

void Simulation::set_bnd ( int N, int b, Mat3D *x ) {
    int i, k;

    for ( i=1 ; i<=N ; i++ ) {
        x->valAt(0,i) = b==1 ? -2*x->valAt(1,i) : (b == 0 ? -2*x->valAt(1,i) : -2*x->valAt(1,i));
        x->valAt(N+1,i) = b==1 ? -2*x->valAt(N,i) : (b == 0 ? -2*x->valAt(N,i) : -2 * x->valAt(N,i));
        x->valAt(i,0 ) = b==2 ? -x->valAt(i,1) : (b == 0 ? -2*x->valAt(i,1) : -2*x->valAt(i,1));
        x->valAt(i,N+1) = b==2 ? -x->valAt(i,N) : (b == 0 ? -2*x->valAt(i,N) : -2*x->valAt(i,N));
    }

    x->valAt(0,0,0)    = 0.333*(x->valAt(1,0,0) + x->valAt(0,1,0) + x->valAt(0,0,1));
    x->valAt(0,N+1,0)  = 0.333*(x->valAt(1,N+1,0) + x->valAt(0,N,0)  + x->valAt(0,N+1,1));
    x->valAt(N+1,0,0)  = 0.333*(x->valAt(N,0,0)  + x->valAt(N+1,1,0) + x->valAt(N+1,0,1));
    x->valAt(N+1,N+1,0) = 0.333*(x->valAt(N,N+1,0) + x->valAt(N+1,N,0) + x->valAt(N+1,N+1,1));

    x->valAt(0,0,N+1)    = 0.333*(x->valAt(1,0,N+1) + x->valAt(0,1,N+1) + x->valAt(0,0,N));
    x->valAt(0,N+1,N+1)  = 0.333*(x->valAt(1,N+1,N+1) + x->valAt(0,N,N+1)  + x->valAt(0,N+1,N));
    x->valAt(N+1,0,N+1)  = 0.333*(x->valAt(N,0,N+1)  + x->valAt(N+1,1,N+1) + x->valAt(N+1,0,N));
    x->valAt(N+1,N+1,N+1) = 0.333*(x->valAt(N,N+1,N+1) + x->valAt(N+1,N,N+1) + x->valAt(N+1,N+1,N));
}

void Simulation::diffuse(Mat3D *x, Mat3D *xprev) {
  int iters, i, j, k;

  float dt = params_.timeStep;
  int n = params_.gridSize;
  float kDiffusion = .001;
  float a=dt*kDiffusion*pow(n, 3);

  // Gauss-Seidel Relaxation
  for ( iters=0 ; iters<20 ; iters++ ) {
    for ( i=1 ; i<=n ; i++ ) {
        for ( j=1 ; j<=n ; j++ ) {
            for ( k=1 ; k<=n ; k++ ) {
                x->valAt(i,j,k) = (xprev->valAt(i,j,k) +
                               a*(x->valAt(i-1,j,k) + x->valAt(i+1,j,k) +
                                  x->valAt(i,j-1,k) + x->valAt(i,j+1,k) +
                                  x->valAt(i,j,k-1) + x->valAt(i,j,k+1)))/(1+6*a);
            }
        }
    }
    set_bnd ( n, 0, x );
  }
}

void Simulation::advect(Mat3D *x, Mat3D *xprev, Mat3D *vx, Mat3D *vy, Mat3D *vz) {
    int i, j, k, i0, i1, j0, j1, k0, k1;
    float xi, yi, zi, s0, s1, t0, t1, u0, u1;

    int n = params_.gridSize;
    float dt = params_.timeStep;

    for ( i=1 ; i<=n; i++ ) {
      for ( j=1 ; j<=n; j++ ) {
          for ( k=1 ; k<=n; k++ ) {

          // Trace particle back
          xi = i - dt*n*vx->valAt(i,j,k);
          yi = j - dt*n*vy->valAt(i,j,k);
          zi = k - dt*n*vz->valAt(i,j,k);

          // Clip
          if (xi<0.5) xi=0.5; if (xi>n+0.5) xi=n+ 0.5; i0=(int)xi; i1=i0+1;
          if (yi<0.5) yi=0.5; if (yi>n+0.5) yi=n+ 0.5; j0=(int)yi; j1=j0+1;
          if (zi<0.5) zi=0.5; if (zi>n+0.5) zi=n+ 0.5; k0=(int)zi; k1=k0+1;

          // Interpolate
          s1 = xi-i0; s0 = 1-s1; t1 = yi-j0; t0 = 1-t1; u1 = zi-k0; u0 = 1-u1;
          x->valAt(i,j,k) =
                  u0*(s0*(t0*xprev->valAt(i0,j0,k0) + t1*xprev->valAt(i0,j1,k0))+
                      s1*(t0*xprev->valAt(i1,j0,k0) + t1*xprev->valAt(i1,j1,k0))) +
                  u1*(s0*(t0*xprev->valAt(i0,j0,k1) + t1*xprev->valAt(i0,j1,k1))+
                      s1*(t0*xprev->valAt(i1,j0,k1) + t1*xprev->valAt(i1,j1,k1)));
      }
     }
    }
    set_bnd ( n, 0, x );
}

void Simulation::project(Mat3D *fluidvx, Mat3D *fluidvy, Mat3D *fluidvz, Mat3D *fluidvx_prev, Mat3D *fluidvy_prev, Mat3D *fluidvz_prev) {
    // Conserves Mass within the System
    int iters, i, j, k;
    float h;

    int n = params_.gridSize;

    h = 1.0/n;

    /*
    for ( i=1 ; i<=n ; i++ ) {
     for ( j=1 ; j<=n ; j++ ) {
     fluidvy_prev->valAt(i,j) = -0.5*h*(fluidvx->valAt(i+1,j)-fluidvx->valAt(i-1,j)+
                                        fluidvy->valAt(i,j+1)-fluidvy->valAt(i,j-1));
     fluidvx_prev->valAt(i,j) = 0;
    }
    }
    set_bnd ( n, 0, fluidvy_prev ); set_bnd ( n, 0, fluidvx_prev );

    for ( iters=0 ; iters<20 ; iters++ ) {
     for ( i=1 ; i<=n ; i++ ) {
     for ( j=1 ; j<=n ; j++ ) {
     fluidvx_prev->valAt(i,j) = (fluidvy_prev->valAt(i,j) + fluidvx_prev->valAt(i-1,j) + fluidvx_prev->valAt(i+1,j) +
                                                     fluidvx_prev->valAt(i,j-1) + fluidvx_prev->valAt(i,j+1))/4;
    }
     }
     set_bnd ( n, 0, fluidvx_prev );
     }

     for ( i=1 ; i<=n ; i++ ) {
     for ( j=1 ; j<=n ; j++ ) {
     fluidvx->valAt(i,j) -= 0.5*(fluidvx_prev->valAt(i+1,j) - fluidvx_prev->valAt(i-1,j))/h;
     fluidvy->valAt(i,j) -= 0.5*(fluidvx_prev->valAt(i,j+1) - fluidvx_prev->valAt(i,j-1))/h;
     }
     }*/

    for ( i=1 ; i<=n ; i++ ) {
        for ( j=1 ; j<=n ; j++ ) {
            for ( k=1 ; k<=n ; k++ ) {
                // Find total flow through pt i,j,k
                fluidvy_prev->valAt(i,j,k) = -0.333*h*(fluidvx->valAt(i+1,j,k)-fluidvx->valAt(i-1,j,k)+
                                                       fluidvy->valAt(i,j+1,k)-fluidvy->valAt(i,j-1,k)+
                                                       fluidvz->valAt(i,j,k+1)-fluidvz->valAt(i,j,k-1));
                fluidvx_prev->valAt(i,j,k) = 0;
                fluidvz_prev->valAt(i,j,k) = 0;
            }
        }
    }

    set_bnd ( n, 0, fluidvy_prev ); set_bnd ( n, 0, fluidvx_prev ); set_bnd ( n, 0, fluidvz_prev);

    for ( iters=0 ; iters<20 ; iters++ ) {
        for ( i=1 ; i<=n ; i++ ) {
            for ( j=1 ; j<=n ; j++ ) {
                for ( k=1 ; k<=n ; k++ ) {
                    fluidvx_prev->valAt(i,j,k) = (fluidvy_prev->valAt(i,j,k) +
                                                fluidvx_prev->valAt(i-1,j,k) + fluidvx_prev->valAt(i+1,j,k) +
                                                fluidvx_prev->valAt(i,j-1,k) + fluidvx_prev->valAt(i,j+1,k) +
                                                fluidvx_prev->valAt(i,j,k-1) + fluidvx_prev->valAt(i,j,k+1))/6;
/*
                    fluidvz_prev->valAt(i,j,k) = (fluidvy_prev->valAt(i,j) +
                                                fluidvz_prev->valAt(i-1,j) + fluidvz_prev->valAt(i+1,j) +
                                                fluidvz_prev->valAt(i,j-1) + fluidvz_prev->valAt(i,j+1) +
                                                fluidvz_prev->valAt(i,j,k-1) + fluidvz_prev->valAt(i,j,k+1))/6;*/
                }
            }
        }
        set_bnd ( n, 0, fluidvx_prev );
     }

     for ( i=1 ; i<=n ; i++ ) {
        for ( j=1 ; j<=n ; j++ ) {
            for ( k=1 ; k<=n ; k++ ) {
                fluidvx->valAt(i,j,k) -= 0.5*(fluidvx_prev->valAt(i+1,j,k) - fluidvx_prev->valAt(i-1,j,k))/h;
                fluidvy->valAt(i,j,k) -= 0.5*(fluidvx_prev->valAt(i,j+1,k) - fluidvx_prev->valAt(i,j-1,k))/h;
                fluidvz->valAt(i,j,k) -= 0.5*(fluidvx_prev->valAt(i,j,k+1) - fluidvx_prev->valAt(i,j,k-1))/h;
            }
        }
     }

     set_bnd ( n, 1, fluidvx ); set_bnd ( n, 2, fluidvy ); set_bnd (n, 3, fluidvz);
}

void Simulation::stableFluidSolve() {
  int i, j, k;
  int n = params_.gridSize;

  double dt = params_.timeStep;

  // Density Diffusion
  SWAP_MAT(fluiddensity_prev, fluiddensity)
  diffuse(fluiddensity, fluiddensity_prev);

  // Density Advection
  SWAP_MAT(fluiddensity_prev, fluiddensity)
  advect(fluiddensity, fluiddensity_prev, fluidvx, fluidvy, fluidvz);

  // Add force grid * timestep to velocity field
  /*for (i = 0 ; i < n+2; i++) {
      for (j = 0 ; j < n+2; j++) {
          for (k = 0 ; k < n+2; k++) {
            fluidvx->valAt(i, j, k) += dt*fluidvx_prev->valAt(i, j, k);
            fluidvy->valAt(i, j, k) += dt*fluidvy_prev->valAt(i, j, k);
            fluidvz->valAt(i, j, k) += dt*fluidvz_prev->valAt(i, j, k);
          }
        }
  }*/

  // Velocity Diffusion
  SWAP_MAT(fluidvx_prev, fluidvx)
  SWAP_MAT(fluidvy_prev, fluidvy)
  SWAP_MAT(fluidvy_prev, fluidvy)

  diffuse(fluidvx, fluidvx_prev);
  diffuse(fluidvy, fluidvy_prev);
  diffuse(fluidvz, fluidvz_prev);

  // project
  project(fluidvx, fluidvy, fluidvz, fluidvx_prev, fluidvy_prev, fluidvz_prev);

   // Velocity Self-Advection
  SWAP_MAT(fluidvx_prev, fluidvx)
  SWAP_MAT(fluidvy_prev, fluidvy)
  SWAP_MAT(fluidvz_prev, fluidvz)
  advect(fluidvx, fluidvx_prev, fluidvx_prev, fluidvy_prev, fluidvz_prev);
  advect(fluidvy, fluidvy_prev, fluidvx_prev, fluidvy_prev, fluidvz_prev);
  advect(fluidvz, fluidvz_prev, fluidvx_prev, fluidvy_prev, fluidvz_prev);

  // project
  project(fluidvx, fluidvy, fluidvz, fluidvx_prev, fluidvy_prev, fluidvz_prev);
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

            VectorXd df = -D.transpose() * (C.transpose() * (A.transpose() * epsilonv));
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
            double theta = 2*atan2(n0xn1.norm(), n0.norm()*n1.norm() + n0.dot(n1));

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
