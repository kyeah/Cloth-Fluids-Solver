#ifndef MAT3D_H
#define MAT3D_H

class Mat3D {

public:
    Mat3D(int nx, int ny, int nz) : rows(nx), cols(ny), depth(nz) {
        mat_ = new double[nx*ny*nz]();  // Initialized to zero
    }

    ~Mat3D() {
        delete[] mat_;
    }

    double& operator[](unsigned int i){ return mat_[i];}
    const double& operator[](unsigned int i)const{ return mat_[i];}

    double& valAt(int i, int j, int k) const {
        return mat_[cols*depth*i + depth*j + k];
    }

    double *mat_;
    int rows, cols, depth;
};

class Mat2D {

public:
    Mat2D(int nx, int ny) : rows(nx), cols(ny) {
        mat_ = new double[nx*ny]();  // Initialized to zero
    }

    ~Mat2D() {
        delete[] mat_;
    }

    double& operator[](unsigned int i){ return mat_[i];}
    const double& operator[](unsigned int i)const{ return mat_[i];}

    double& valAt(int i, int j) const {
        return mat_[cols*i + j];
    }

    double *mat_;
    int rows, cols;
};


#endif // MAT3D_H
