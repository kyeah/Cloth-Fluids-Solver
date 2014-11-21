#include "cloth.h"
#include <QGLWidget>
#include <Eigen/Geometry>
#include <iostream>

using namespace std;
using namespace Eigen;

Cloth::Cloth() {
    mesh_ = new Mesh("resources/square.obj");
    mesh_->translate(Vector3d(5,0,4));
    int num_verts = mesh_->getNumVerts();
    verts_.resize(3*num_verts);
    velocities_.resize(3*num_verts);
    resetState();
    computeMassMatrices();
}

void Cloth::computeMassMatrices() {
    int numVerts = mesh_->getNumVerts();
    VectorXd vertMasses(numVerts);

    for(int i=0; i < mesh_->getNumFaces(); i++)
    {
        Vector3i fverts = mesh_->getFace(i);
        Vector3d pts[3];
        for(int j=0; j<3; j++)
            pts[j] = verts_.segment<3>(3*fverts[j]);

        Vector3d normal = (pts[1] - pts[0]).cross(pts[2] - pts[0]);
        double norm = normal.norm();
        double area = 0.5*norm;

        for(int j=0; j<3; j++)
        {
            vertMasses[fverts[j]] += area;
        }
    }

    mass_.resize(3*numVerts, 3*numVerts);
    massInv_.resize(3*numVerts, 3*numVerts);
    for (int i = 0; i < numVerts; i++) {
        vertMasses[i] /= 3.;

        for (int j = 0; j < 3; j++) {
            mass_(i+j, i+j) = vertMasses[i];
            massInv_(i+j, i+j) = 1. / vertMasses[i];
        }
    }
}

VectorXd Cloth::getVertNormals()
{
    VectorXd vertNormals(3*mesh_->getNumVerts());
    vector<double> verttotalarea;

    vertNormals.setZero();
    verttotalarea.clear();

    for(int i=0; i < mesh_->getNumVerts(); i++)
    {
        verttotalarea.push_back(0);
    }

    for(int i=0; i < mesh_->getNumFaces(); i++)
    {
        Vector3i fverts = mesh_->getFace(i);
        Vector3d pts[3];
        for(int j=0; j<3; j++)
            pts[j] = verts_.segment<3>(3*fverts[j]);

        Vector3d normal = (pts[1] - pts[0]).cross(pts[2] - pts[0]);
        double norm = normal.norm();
        for(int j=0; j<3; j++)
        {
            vertNormals.segment<3>(3*fverts[j]) += normal;
            verttotalarea[fverts[j]] += norm;
        }
    }

    for(int i=0; i<(int)verts_.size()/3; i++)
        vertNormals.segment<3>(3*i) /= verttotalarea[i];

    return vertNormals;
}

void Cloth::render()
{
    glShadeModel(GL_SMOOTH);
    glDisable(GL_TEXTURE_2D);
    glEnable(GL_LIGHTING);
    glColorMaterial ( GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE );
    glEnable ( GL_COLOR_MATERIAL );
    Vector3d color(0.9, 0.6, 0.9);
    glColor4d(color[0], color[1], color[2], 1.0);

    VectorXd vertNormals = getVertNormals();
    double *vertNormalsPointer = &vertNormals[0];
    double *vertPointer = &verts_[0];

    glPushMatrix();
    {
        glEnableClientState(GL_VERTEX_ARRAY);
        glEnableClientState(GL_NORMAL_ARRAY);

        glVertexPointer(3, GL_DOUBLE, 0, vertPointer);
        glNormalPointer(GL_DOUBLE, 0, vertNormalsPointer);

        glDrawElements(GL_TRIANGLES, mesh_->getNumFaces()*3, GL_UNSIGNED_INT, mesh_->getFacePointer());

        glDisableClientState(GL_VERTEX_ARRAY);
        glDisableClientState(GL_NORMAL_ARRAY);
    }
    glPopMatrix();
}

void Cloth::resetState() {
    velocities_.setZero();

    for (int i = 0; i < mesh_->getNumVerts(); i++) {
        verts_.segment(3*i, 3) = mesh_->getVert(i);
    }
}
