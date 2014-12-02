#include "simulation.h"
#include <QGLWidget>
#include "simparameters.h"
#include <iostream>
#include <Eigen/Geometry>
#include <QDebug>
#include "SOIL.h"
#include "vectormath.h"
#include <Eigen/Dense>

// Quick macro to swap two matrix pointers
#define SWAP_MAT(A,B) {Mat3D *tmp=A;A=B;B=tmp;}

const double PI = 3.1415926535898;

using namespace Eigen;
using namespace std;

Simulation::Simulation(const SimParameters &params) : params_(params), time_(0), floorTex_(0)
{
    fluidvx = fluidvy = fluidvz = fluidvx_prev = fluidvy_prev = fluidvz_prev = fluiddensity = fluiddensity_prev = NULL;
    clearScene();
}

Simulation::~Simulation()
{
    delete fluidvx;
    delete fluidvy;
    delete fluidvz;
    delete fluidvx_prev;
    delete fluidvy_prev;
    delete fluidvz_prev;
    delete fluiddensity;
    delete fluiddensity_prev;
}

void Simulation::initializeGL()
{
    loadFloorTexture();
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
        float scale = 0.5;
        glPointSize(40*(1/scale));
        int scaledGridSize = params_.gridSize * scale;

        int pts[2] = { 0, scaledGridSize };
        glBegin(GL_LINES);
        glColor3f(0,0,0);
        for (int xi = 0; xi < 2; xi++) {
            for (int yi = 0; yi < 2; yi++) {
                int x = pts[xi];
                int y = pts[yi];
                glVertex3f(0, x, y);
                glVertex3f(scaledGridSize, x, y);
                glVertex3f(x, y, 0);
                glVertex3f(x, y, scaledGridSize);
                glVertex3f(x, 0, y);
                glVertex3f(x, scaledGridSize, y);
            }
        }
        glEnd();

        glPointSize(15);
        glBegin(GL_POINTS);

        for (int i = 0; i < params_.gridSize; i++) {
            for (int j = 0; j < params_.gridSize; j++) {
                for (int k = 0; k < params_.gridSize; k++) {

                float val = (fluiddensity->valAt(i,j,k)) * 2550;
                glColor4f(val, 0, val, val);
                glVertex3f(scale*(params_.gridSize - i), scale*k, scale*j);
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
    stableFluidSolve();
}

void Simulation::addVelocity(Eigen::Vector2d pos, Eigen::Vector2d vel) {
    int px = floor((pos[0] + 1)*params_.gridSize / 2);
    int py = floor((pos[1] + 1)*params_.gridSize / 2);
    vel = 100 * (vel / vel.norm());

    if (px < 0 || px >= params_.gridSize - 1 || py < 0 || py >= params_.gridSize - 1) {
        return;
    }

    if (params_.applyFluidDragForce) {
        fluidvx->valAt(px, py, params_.gridSize/2) += vel[0];
        fluidvy->valAt(px, py, params_.gridSize/2) += vel[1];

    } else {
        fluiddensity->valAt(px, py, params_.gridSize/2) += 100;
        fluidvx->valAt(px, py, params_.gridSize/2) += VectorMath::randomUnitIntervalReal() * 980;
        fluidvy->valAt(px, py, params_.gridSize/2) += VectorMath::randomUnitIntervalReal() * 980;
        fluidvz->valAt(px, py, params_.gridSize/2) += VectorMath::randomUnitIntervalReal() * 980;
        fluidvx_prev->valAt(px, py, params_.gridSize/2) += VectorMath::randomUnitIntervalReal() * 980;
        fluidvy_prev->valAt(px, py, params_.gridSize/2) += VectorMath::randomUnitIntervalReal() * 980;
        fluidvz_prev->valAt(px, py, params_.gridSize/2) += VectorMath::randomUnitIntervalReal() * 980;
    }

    bound_mat ( params_.gridSize, 0, fluidvx );
    bound_mat ( params_.gridSize, 0, fluidvy );
    bound_mat ( params_.gridSize, 0, fluidvz );
}

void Simulation::clearScene()
{
    renderLock_.lock();
    {
        delete fluidvx;
        delete fluidvy;
        delete fluidvz;
        delete fluidvx_prev;
        delete fluidvy_prev;
        delete fluidvz_prev;
        delete fluiddensity;
        delete fluiddensity_prev;
        fluidvx = new Mat3D(params_.gridSize+2, params_.gridSize+2, params_.gridSize+2);
        fluidvy = new Mat3D(params_.gridSize+2, params_.gridSize+2, params_.gridSize+2);
        fluidvz = new Mat3D(params_.gridSize+2, params_.gridSize+2, params_.gridSize+2);
        fluidvx_prev = new Mat3D(params_.gridSize+2, params_.gridSize+2, params_.gridSize+2);
        fluidvy_prev = new Mat3D(params_.gridSize+2, params_.gridSize+2, params_.gridSize+2);
        fluidvz_prev = new Mat3D(params_.gridSize+2, params_.gridSize+2, params_.gridSize+2);
        fluiddensity = new Mat3D(params_.gridSize+2, params_.gridSize+2, params_.gridSize+2);
        fluiddensity_prev = new Mat3D(params_.gridSize+2, params_.gridSize+2, params_.gridSize+2);
    }
    renderLock_.unlock();
}

void Simulation::bound_mat ( int N, int b, Mat3D *x ) {
    /* Sets non-wrapping boundaries along the fluid container.
     *
     * @params:
     *   N: Gridsize, not including the boundary
     *   b: bounding wrap. 0=passthru, 1=no-xwrap, 2=no-ywrap, 3=no-zwrap
     */

    int i, k;

    // Set non-wrapping boundaries
    for ( i=1 ; i<=N ; i++ ) {
        for (k=1; k <= N; k++) {
            // Sides aligned with x axis
            x->valAt(0,i,k) = b==1 ? -x->valAt(1,i,k) : x->valAt(1,i,k);
            x->valAt(N+1,i,k) = b==1 ? -x->valAt(N,i,k) : x->valAt(N,i,k);

            // Sides aligned with y axis
            x->valAt(i,0,k) = b==2 ? -x->valAt(i,1,k) : x->valAt(i,1,k);
            x->valAt(i,N+1,k) = b==2 ? -x->valAt(i,N,k) : x->valAt(i,N,k);

            // Sides aligned with z axis
            x->valAt(i,k,0) = b==3 ? -x->valAt(i,1,k) : x->valAt(i,k,1);
            x->valAt(i,k,N+1) = b==3 ? -x->valAt(i,k,N) : x->valAt(i,k,N);
        }
    }

    // Average the corners with neighboring perimeter values
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

  int n = params_.gridSize;
  float dt = params_.timeStep;
  float kDiffusion = params_.kDiffusion;
  float coeff=dt*kDiffusion*pow(n, 3);

  // Gauss-Seidel Relaxation
  for ( iters=0 ; iters<10 ; iters++ ) {
    for ( i=1 ; i<=n ; i++ ) {
        for ( j=1 ; j<=n ; j++ ) {
            for ( k=1 ; k<=n ; k++ ) {
                x->valAt(i,j,k) = (xprev->valAt(i,j,k) +
                               coeff*(x->valAt(i-1,j,k) + x->valAt(i+1,j,k) +
                                      x->valAt(i,j-1,k) + x->valAt(i,j+1,k) +
                                      x->valAt(i,j,k-1) + x->valAt(i,j,k+1)))/(1+6*coeff);
            }
        }
    }
    bound_mat(n, 0, x);
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

          // Clamp indexes
          if (xi<0.5) xi=0.5; if (xi>n+0.5) xi=n+0.5;
          i0=(int)xi; i1=i0+1;

          if (yi<0.5) yi=0.5; if (yi>n+0.5) yi=n+0.5;
          j0=(int)yi; j1=j0+1;

          if (zi<0.5) zi=0.5; if (zi>n+0.5) zi=n+0.5;
          k0=(int)zi; k1=k0+1;

          // Tri-linear interpolation using previous values of x
          s1 = xi-i0; s0 = 1-s1;
          t1 = yi-j0; t0 = 1-t1;
          u1 = zi-k0; u0 = 1-u1;

          x->valAt(i,j,k) =
                  u0*(s0*(t0*xprev->valAt(i0,j0,k0) + t1*xprev->valAt(i0,j1,k0))+
                      s1*(t0*xprev->valAt(i1,j0,k0) + t1*xprev->valAt(i1,j1,k0))) +
                  u1*(s0*(t0*xprev->valAt(i0,j0,k1) + t1*xprev->valAt(i0,j1,k1))+
                      s1*(t0*xprev->valAt(i1,j0,k1) + t1*xprev->valAt(i1,j1,k1)));
      }
     }
    }

    bound_mat (n, 0, x);
}

void Simulation::project(Mat3D *fluidvx, Mat3D *fluidvy, Mat3D *fluidvz) {
    /*
     * Conserves Mass within the System.
     */

    int iters, i, j, k;

    int n = params_.gridSize;
    float h = 1.0/n;  // Grid cell height

    Mat3D cell_flow(fluidvy->rows, fluidvy->cols, fluidvy->depth);
    Mat3D conservation_mat(fluidvx->rows, fluidvx->cols, fluidvx->depth);

    for ( i=1 ; i<=n ; i++ ) {
        for ( j=1 ; j<=n ; j++ ) {
            for ( k=1 ; k<=n ; k++ ) {
                // Find total flow through pt i,j,k
                cell_flow.valAt(i,j,k) = -0.333*h*(fluidvx->valAt(i+1,j,k)-fluidvx->valAt(i-1,j,k)+
                                                   fluidvy->valAt(i,j+1,k)-fluidvy->valAt(i,j-1,k)+
                                                   fluidvz->valAt(i,j,k+1)-fluidvz->valAt(i,j,k-1));
            }
        }
    }

    bound_mat (n, 0, &cell_flow);

    // Solve Poisson Equation using Gauss-Seidel relaxation
    for ( iters=0 ; iters<10 ; iters++ ) {
        for ( i=1 ; i<=n ; i++ ) {
            for ( j=1 ; j<=n ; j++ ) {
                for ( k=1 ; k<=n ; k++ ) {
                    conservation_mat.valAt(i,j,k) = (cell_flow.valAt(i,j,k) +
                                         conservation_mat.valAt(i-1,j,k) + conservation_mat.valAt(i+1,j,k) +
                                         conservation_mat.valAt(i,j-1,k) + conservation_mat.valAt(i,j+1,k) +
                                         conservation_mat.valAt(i,j,k-1) + conservation_mat.valAt(i,j,k+1))/6;
                }
            }
        }
        bound_mat (n, 0, &conservation_mat);
     }

    // Subtract gradient map from velocity matrices
     for ( i=1 ; i<=n ; i++ ) {
        for ( j=1 ; j<=n ; j++ ) {
            for ( k=1 ; k<=n ; k++ ) {
                fluidvx->valAt(i,j,k) -= 0.5*(conservation_mat.valAt(i+1,j,k) - conservation_mat.valAt(i-1,j,k))/h;
                fluidvy->valAt(i,j,k) -= 0.5*(conservation_mat.valAt(i,j+1,k) - conservation_mat.valAt(i,j-1,k))/h;
                fluidvz->valAt(i,j,k) -= 0.5*(conservation_mat.valAt(i,j,k+1) - conservation_mat.valAt(i,j,k-1))/h;
            }
        }
     }

     bound_mat (n, 1, fluidvx); bound_mat (n, 2, fluidvy); bound_mat (n, 3, fluidvz);
}

void Simulation::stableFluidSolve() {
  int n = params_.gridSize;

  // Density Diffusion
  SWAP_MAT(fluiddensity_prev, fluiddensity)
  diffuse(fluiddensity, fluiddensity_prev);

  // Density Advection
  SWAP_MAT(fluiddensity_prev, fluiddensity)
  advect(fluiddensity, fluiddensity_prev, fluidvx, fluidvy, fluidvz);

  // Add forces to velocity
  for (int i = 1; i <= n; i++) {
      for (int j = 1; j <= n; j++) {
          for (int k = 1; k <= n; k++) {
              fluidvy_prev->valAt(i, j, k) += params_.gravityG;
              fluidvy->valAt(i, j, k) += params_.gravityG;
          }
      }
  }

  bound_mat (n, 1, fluidvx);
  bound_mat (n, 2, fluidvy);
  bound_mat (n, 3, fluidvz);

  // Velocity Diffusion
  SWAP_MAT(fluidvx_prev, fluidvx)
  SWAP_MAT(fluidvy_prev, fluidvy)
  SWAP_MAT(fluidvy_prev, fluidvy)

  diffuse(fluidvx, fluidvx_prev);
  diffuse(fluidvy, fluidvy_prev);
  diffuse(fluidvz, fluidvz_prev);
  bound_mat (n, 1, fluidvx);
  bound_mat (n, 2, fluidvy);
  bound_mat (n, 3, fluidvz);

  // Pre-Projection
  project(fluidvx, fluidvy, fluidvz);

   // Velocity Self-Advection
  SWAP_MAT(fluidvx_prev, fluidvx)
  SWAP_MAT(fluidvy_prev, fluidvy)
  SWAP_MAT(fluidvz_prev, fluidvz)
  advect(fluidvx, fluidvx_prev, fluidvx_prev, fluidvy_prev, fluidvz_prev);
  advect(fluidvy, fluidvy_prev, fluidvx_prev, fluidvy_prev, fluidvz_prev);
  advect(fluidvz, fluidvz_prev, fluidvx_prev, fluidvy_prev, fluidvz_prev);
  bound_mat (n, 1, fluidvx);
  bound_mat (n, 2, fluidvy);
  bound_mat (n, 3, fluidvz);

  // Post-Projection
  project(fluidvx, fluidvy, fluidvz);

}
