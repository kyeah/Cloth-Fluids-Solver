#include "controller.h"
#include "mainwindow.h"
#include "simulation.h"
#include <QDebug>
#include <Eigen/Core>

using namespace Eigen;

Controller::Controller(int fps) : QThread(), mw_(NULL), fps_(fps)
{
}

Controller::~Controller()
{
    delete sim_;
}

void Controller::initialize(MainWindow *mw)
{
    mw_ = mw;
    sim_ = new Simulation(params_);
}

void Controller::initializeGL()
{
    sim_->initializeGL();
}

void Controller::run()
{
    reset();
    connect(&simtimer_, SIGNAL(timeout()), this, SLOT(simTick()));
    simtimer_.start(1000/fps_);
    exec();
}

void Controller::reset()
{
    params_ = SimParameters();
    QMetaObject::invokeMethod(mw_, "setUIFromParameters", Q_ARG(SimParameters, params_));
    clearScene();
}

void Controller::clearScene()
{
    sim_->clearScene();
}

void Controller::updateParameters(SimParameters params)
{
    params_ = params;
}

void Controller::renderFloor()
{
    sim_->renderFloor();
}

void Controller::renderObjects()
{
    sim_->renderObjects();
}

void Controller::mouseClicked(double, double, double, double , double , double )
{
}

void Controller::accelerateBody(double vx, double vy, double vz, double wx, double wy, double wz)
{
    sim_->accelerateBody(vx,vy,vz, wx, wy, wz);
}

void Controller::simTick()
{
    if(params_.simRunning)
    {
        for(int i=0; i<10; i++)
            sim_->takeSimulationStep();
    }
}
