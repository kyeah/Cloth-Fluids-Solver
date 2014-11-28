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

const double PI = 3.1415926535898;

using namespace Eigen;
using namespace std;

Simulation::Simulation(const SimParameters &params) : params_(params), time_(0), floorTex_(0)
{
    loadRigidBodies();
    cloth_ = new Cloth();
    bodyInstance_ = NULL;
    fluidvx = fluidvy = fluidfx = fluidfy = fluiddensity = NULL;
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
        glPointSize(5);
        glBegin(GL_POINTS);
        for (int i = 0; i < params_.gridSize; i++) {
            for (int j = 0; j < params_.gridSize; j++) {
                float val = fluidvx->valAt(i,j);
                float valy = fluidvy->valAt(i,j);
                while (val > 255) {
                    val -= 255;
                }
                while (valy > 255) {
                    valy -= 255;
                }
                glColor3f(valy, 0, val);
                glVertex3f(i, 0, j + 5);
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

    for (int i = 0; i < params_.gridSize; i++) {
        for (int j = 0; j < params_.gridSize; j++) {
            fluidfx->valAt(i, j) = -9.8;
        }
    }
    fluidfx->valAt(1, params_.gridSize - 1) = 20*VectorMath::randomUnitIntervalReal();
    fluidfy->valAt(1, params_.gridSize - 1) = 60*VectorMath::randomUnitIntervalReal();
    fluidfx->valAt(2, params_.gridSize - 1) = 20*VectorMath::randomUnitIntervalReal();
    fluidfy->valAt(2, params_.gridSize - 1) = 60*VectorMath::randomUnitIntervalReal();
    fluidfx->valAt(3, params_.gridSize - 1) = 20*VectorMath::randomUnitIntervalReal();
    fluidfy->valAt(3, params_.gridSize - 1) = 60*VectorMath::randomUnitIntervalReal();
    for (int i = 0; i < params_.gridSize; i++) {
        fluidfx->valAt(i, 0) = -98*3;
        fluidfy->valAt(i, 0) = 5;
    }
    stableFluidSolve(*fluidvx, *fluidvy, *fluidfx, *fluidfy);
}

void Simulation::clearScene()
{
    renderLock_.lock();
    {
        delete bodyInstance_;
        delete fluidvx;
        delete fluidvy;
        delete fluidfx;
        delete fluidfy;
        delete fluiddensity;
        Vector3d pos(5, 0, 3);
        Vector3d zero(0,0,0);
        bodyInstance_ = new RigidBodyInstance(*bodyTemplate_, pos, zero, 1.0);
        cloth_->resetState();
        fluidvx = new Mat2D(params_.gridSize, params_.gridSize);
        fluidvy = new Mat2D(params_.gridSize, params_.gridSize);
        fluidfx = new Mat2D(params_.gridSize, params_.gridSize);
        fluidfy = new Mat2D(params_.gridSize, params_.gridSize);
        fluiddensity = new Mat2D(params_.gridSize, params_.gridSize);
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

void Simulation::stableFluidSolve(Mat2D &u, Mat2D &v, Mat2D &u0, Mat2D &v0) {
  float x, y, x0, y0, f, r, U[2], V[2], s, t;
  int i, j, i0, j0, i1, j1;

  float dt = params_.timeStep;
  float visc = params_.viscosity;
  int n = params_.gridSize;

  // u, v = previous velocities (vx, vy)
  // u0, v0 = forces (fx, fy)

  // Add force grid * timestep to velocity field
  for (i = 0 ; i < pow(n, 2); i++) {
      u[i] += dt*u0[i]; u0[i] = u[i];
      v[i] += dt*v0[i]; v0[i] = v[i];
  }

  // Self-advection
  for ( x=0.5/n, i=0; i < n; i++, x += 1.0/n) {
      for ( y=0.5/n, j=0 ; j < n; j++, y += 1.0/n) {

          x0 = n*(x-dt*u0[i+n*j])-0.5; y0 = n*(y-dt*v0[i+n*j])-0.5;

          i0 = floor(x0); s = x0 -i0; i0 = (n+(i0%n))%n; i1 = (i0+1);
          j0 = floor(y0); t = y0 -j0; j0 = (n+(j0%n))%n; j1 = (j0+1);

          u[i+n*j] = (1-s)*((1-t)*u0[i0+n*j0]+t*u0[i0+n*j1])+
                        s *((1-t)*u0[i1+n*j0]+t*u0[i1+n*j1]);

          v[i+n*j] = (1-s)*((1-t)*v0[i0+n*j0]+t*v0[i0+n*j1])+
                        s *((1-t)*v0[i1+n*j0]+t*v0[i1+n*j1]);
      }
  }


  /*
  // Transform to Fourier Domain
  for ( i=0 ; i<n ; i++ )
    for ( j=0 ; j<n ; j++ )
      { u0[i+(n+2)*j] = u[i+n*j]; v0[i+(n+2)*j] = v[i+n*j]; }
*/

  FFT<double> fft;
  MatrixXcd Uf(n,n), Vf(n,n);
  for (int i = 0; i < n; i++) {
      VectorXd A(n);
      for (int j = 0; j < n; j++) {
          A[j] = u0[i+n*j];
      }
      VectorXcd B(n);
      fft.fwd(B, A);
      Uf.block(0,i,n,1) = B;
  }
  for (int i = 0; i < n; i++) {
      VectorXd A(n);
      for (int j = 0; j < n; j++) {
          A[j] = v0[i+n*j];
      }
      VectorXcd B(n);
      fft.fwd(B, A);
      Vf.block(0,i,n,1) = B;
  }

  // Diffusion and Projection (Mass Conservation)
  for ( i=0 ; i<=n ; i+=2 ) {
      x = 0.5*i;

      for ( j=0 ; j<n ; j++ ) {
          y = j<=n/2 ? j : j-n;

          r = x*x+y*y;
          if ( r==0.0 ) continue;

          f = exp(-r*dt*visc);
          U[0] = u0[i +(n)*j]; V[0] = v0[i +(n)*j];
          U[1] = u0[i+1+(n)*j]; V[1] = v0[i+1+(n)*j];

          u0[i +(n)*j] = f*( (1-x*x/r)*U[0] -x*y/r *V[0] );
          u0[i+1+(n)*j] = f*( (1-x*x/r)*U[1] -x*y/r *V[1] );

          v0[i+ (n)*j] = f*( -y*x/r *U[0] + (1-y*y/r)*V[0] );
          v0[i+1+(n)*j] = f*( -y*x/r *U[1] + (1-y*y/r)*V[1] );
      }
  }

  // Reconvert into Cartesian
  for (int i = 0; i < n; i++) {
      VectorXcd A(n);
      for (int j = 0; j < n; j++) {
          A[j] = u0[i+n*j];
      }
      VectorXd B(n);
      fft.inv(B, A);
      for (int j = 0; j < n; j++) {
          u0[i+n*j] = B[j];
      }
  }
  for (int i = 0; i < n; i++) {
      VectorXcd A(n);
      for (int j = 0; j < n; j++) {
          A[j] = v0[i+n*j];
      }
      VectorXd B(n);
      fft.inv(B, A);
      for (int j = 0; j < n; j++) {
          v0[i+n*j] = B[j];
      }
  }

  //fft(-1,n,u0); FFT(-1,n,v0);
  f = 1.0/(n*n);
/*
  for ( i=0 ; i<n ; i++ )
    for ( j=0 ; j<n ; j++ )
      { u[i+n*j] = f*u0[i+(n+2)*j]; v[i+n*j] = f*v0[i+(n+2)*j ]; }
*/
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
