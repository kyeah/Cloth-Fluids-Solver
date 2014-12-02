#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "simparameters.h"
#include "controller.h"

MainWindow::MainWindow(Controller &cont, int fps, QWidget *parent) :
    QMainWindow(parent),
    cont_(cont),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    ui->GLWidget->setController(&cont);
    simRunning_ = false;
    connect(&renderTimer_, SIGNAL(timeout()), this, SLOT(updateGL()));
    renderTimer_.start(1000/fps);
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::on_actionExit_triggered()
{
    close();
}

void MainWindow::setParametersFromUI()
{
    SimParameters params;

    params.simRunning = simRunning_;

    params.timeStep = ui->timeStepEdit->text().toDouble();
    params.applyFluidDragForce = ui->fluidDragForceCheckbox->isChecked();

    params.activeForces = 0;
    if(ui->gravityCheckBox->isChecked())
        params.activeForces |= SimParameters::F_GRAVITY;

    params.gravityG = ui->gravityGEdit->text().toDouble();
    params.kDiffusion = ui->kDiffEdit->text().toDouble();
    params.dragForceMag = ui->kDragMagEdit->text().toDouble();

    setUIFromParameters(params);
    QMetaObject::invokeMethod(&cont_, "updateParameters", Q_ARG(SimParameters, params));
}

void MainWindow::setUIFromParameters(const SimParameters &params)
{
    if(params.simRunning)
    {
        ui->startSimulationButton->setText(QString("Pause Simulation"));
        simRunning_ = true;
    }
    else
    {
        ui->startSimulationButton->setText(QString("Start Simulation"));
        simRunning_ = false;
    }

    ui->timeStepEdit->setText(QString::number(params.timeStep));
    ui->fluidDragForceCheckbox->setChecked(params.applyFluidDragForce);

    ui->gravityCheckBox->setChecked(params.activeForces & SimParameters::F_GRAVITY);

    ui->gravityGEdit->setText(QString::number(params.gravityG));
    ui->kDiffEdit->setText(QString::number(params.kDiffusion));
    ui->kDragMagEdit->setText(QString::number(params.dragForceMag));
}

void MainWindow::updateGL()
{
    ui->GLWidget->tick();
    ui->GLWidget->update();
}

void MainWindow::on_actionReset_Everything_triggered()
{
    QMetaObject::invokeMethod(&cont_, "reset");
}

void MainWindow::on_actionReset_triggered()
{
    QMetaObject::invokeMethod(&cont_, "clearScene");
}

void MainWindow::on_startSimulationButton_clicked()
{
    simRunning_ = !simRunning_;
    setParametersFromUI();
}

void MainWindow::on_timeStepEdit_editingFinished()
{
    setParametersFromUI();
}

void MainWindow::on_gravityCheckBox_clicked()
{
    setParametersFromUI();
}

void MainWindow::on_gravityGEdit_editingFinished()
{
    setParametersFromUI();
}

void MainWindow::on_fluidDragForceCheckbox_clicked()
{
    setParametersFromUI();
}

void MainWindow::on_kDragMagEdit_editingFinished()
{
    setParametersFromUI();
}
