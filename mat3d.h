#ifndef MAT3D_H
#define MAT3D_H

class Mat3D {

public:
    Mat3D(int nx, int ny, int nz) {
        mat_ = new double[nx*ny*nz]();  // Initialized to zero
    }

    ~Mat3D() {
        delete[] mat_;
    }

    double& operator[](unsigned int i){ return mat_[i];}
    const double& operator[](unsigned int i)const{ return mat_[i];}

    double& valAt(int i, int j, int k) const {
        return mat_[rows*i + cols*j + k];
    }

    double *mat_;
    int rows, cols, depth;
};

#endif // MAT3D_H
