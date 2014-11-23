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
    calculateAmats();
    fillHinges();
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

void Cloth::calculateAmats() {
    amats_.clear();
    bmats_.clear();
    gmats_.clear();
    for (int i = 0; i < mesh_->getNumFaces(); i++) {
        Vector3i fverts = mesh_->getFace(i);
        Vector3d pts[3];
        for(int j=0; j<3; j++)
            pts[j] = verts_.segment<3>(3*fverts[j]);

        Vector3d e1 = pts[1] - pts[0];
        Vector3d e2 = pts[2] - pts[0];
        MatrixXd e(2,3);
        e.block<1,3>(0,0) = e1.transpose();
        e.block<1,3>(1,0) = e2.transpose();

        Matrix2d g;
        g(0,0) = e1.dot(e1);
        g(0,1) = e1.dot(e2);
        g(1,0) = g(0,1);
        g(1,1) = e2.dot(e2);

        MatrixXd bmat = g.inverse() * e;
        MatrixXd A(9,4);
        A.setZero();
        for (int i = 0; i < 9; i++) {
            int col_idx = (i/3);
            int col_idx2 = (i%3);
            double b1 = bmat(0, col_idx);
            double b2 = bmat(1, col_idx);
            double b3 = bmat(0, col_idx2);
            double b4 = bmat(1, col_idx2);
            A(i, 0) = b1*b3;
            A(i, 1) = b1*b4;
            A(i, 2) = b2*b3;
            A(i, 3) = b2*b4;
        }

        gmats_.push_back(g/1);
        bmats_.push_back(g.inverse() * e);
        amats_.push_back(A/1);
    }
}

void Cloth::fillHinges() {

    // For each face
    cout << mesh_->getNumFaces() << endl;
    for (int faceidx = 0; faceidx < mesh_->getNumFaces(); faceidx++) {
        Vector3i fverts = mesh_->getFace(faceidx);

        // For each other face
        for (int faceidx2 = faceidx+1; faceidx2 < mesh_->getNumFaces(); faceidx2++) {
            Vector3i f2verts = mesh_->getFace(faceidx2);

            // For each edge
            for (int i = 0; i < 3; i++) {
                for (int j = i+1; j < 3; j++) {

                    // For each other edge
                    for (int i2 = 0; i2 < 3; i2++) {
                        for (int j2 = 0; j2 < 3; j2++) {
                            int v1i = fverts[i];
                            int v1j = fverts[j];
                            int v2i = f2verts[i2];
                            int v2j = f2verts[j2];

                            if ((v1i==v2i && v1j==v2j) || (v1i==v2j && v1j==v2i)) {
                                Vector3d pi = getVert(v1i);
                                Vector3d pj = getVert(v1j);
                                double restLengthSq = (pj - pi).squaredNorm();
                                double totalArea = mesh_->getFaceArea(faceidx) + mesh_->getFaceArea(faceidx2);

                                hinges_.push_back(Hinge(v1i, v1j, fverts[3-i-j], f2verts[3-i2-j2], faceidx, faceidx2, restLengthSq, totalArea));
                            }
                        }
                    }

                }
            }
        }
    }
}
