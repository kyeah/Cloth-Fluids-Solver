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
    params.penaltyStiffness = ui->penaltyStiffnessEdit->text().toDouble();
    params.pinCorner = ui->pinClothCornersCheckbox->isChecked();

    params.activeForces = 0;
    if(ui->gravityCheckBox->isChecked())
        params.activeForces |= SimParameters::F_GRAVITY;
    if(ui->clothStretchingCheckBox->isChecked())
        params.activeForces |= SimParameters::F_STRETCHING;
    if(ui->clothBendingCheckbox->isChecked())
        params.activeForces |= SimParameters::F_BENDING;
    if(ui->bodyClothContactCheckbox->isChecked())
        params.activeForces |= SimParameters::F_CONTACT;
    if(ui->massDampingCheckBox->isChecked())
        params.activeForces |= SimParameters::F_DAMPING;

    params.gravityG = ui->gravityGEdit->text().toDouble();
    params.stretchingK = ui->stretchingStiffnessEdit->text().toDouble();
    params.bendingK = ui->bendingStiffnessEdit->text().toDouble();
    params.dampingCoeff = ui->dampingCoeffEdit->text().toDouble();

    params.cor = ui->corEdit->text().toDouble();

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
    ui->penaltyStiffnessEdit->setText(QString::number(params.penaltyStiffness));
    ui->pinClothCornersCheckbox->setChecked(params.pinCorner);

    ui->gravityCheckBox->setChecked(params.activeForces & SimParameters::F_GRAVITY);
    ui->clothStretchingCheckBox->setChecked(params.activeForces & SimParameters::F_STRETCHING);
    ui->clothBendingCheckbox->setChecked(params.activeForces & SimParameters::F_BENDING);
    ui->bodyClothContactCheckbox->setChecked(params.activeForces & SimParameters::F_CONTACT);
    ui->massDampingCheckBox->setChecked(params.activeForces & SimParameters::F_DAMPING);

    ui->gravityGEdit->setText(QString::number(params.gravityG));
    ui->stretchingStiffnessEdit->setText(QString::number(params.stretchingK));
    ui->bendingStiffnessEdit->setText(QString::number(params.bendingK));
    ui->corEdit->setText(QString::number(params.cor));   
    ui->dampingCoeffEdit->setText(QString::number(params.dampingCoeff));
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

void MainWindow::on_penaltyStiffnessEdit_editingFinished()
{
    setParametersFromUI();
}

void MainWindow::on_corEdit_editingFinished()
{
    setParametersFromUI();
}

void MainWindow::on_pinClothCornersCheckbox_clicked()
{
    setParametersFromUI();
}

void MainWindow::on_clothStretchingCheckBox_clicked()
{
    setParametersFromUI();
}

void MainWindow::on_clothBendingCheckbox_clicked()
{
    setParametersFromUI();
}

void MainWindow::on_stretchingStiffnessEdit_editingFinished()
{
    setParametersFromUI();
}

void MainWindow::on_bendingStiffnessEdit_editingFinished()
{
    setParametersFromUI();
}

void MainWindow::on_bodyClothContactCheckbox_clicked()
{
    setParametersFromUI();
}

void MainWindow::on_massDampingCheckBox_clicked()
{
    setParametersFromUI();
}

void MainWindow::on_dampingCoeffEdit_editingFinished()
{
    setParametersFromUI();
}
