#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QTimer>

class Controller;
struct SimParameters;

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(Controller &cont, int fps, QWidget *parent = 0);
    ~MainWindow();   

public slots:
    void setUIFromParameters(const SimParameters &params);

private slots:
    void updateGL();

    void on_actionExit_triggered();

    void on_actionReset_Everything_triggered();

    void on_actionReset_triggered();

    void on_startSimulationButton_clicked();

    void on_timeStepEdit_editingFinished();

    void on_gravityCheckBox_clicked();

    void on_gravityGEdit_editingFinished();

    void on_penaltyStiffnessEdit_editingFinished();

    void on_corEdit_editingFinished();

    void on_pinClothCornersCheckbox_clicked();

    void on_clothStretchingCheckBox_clicked();

    void on_clothBendingCheckbox_clicked();

    void on_stretchingStiffnessEdit_editingFinished();

    void on_bendingStiffnessEdit_editingFinished();

    void on_bodyClothContactCheckbox_clicked();

    void on_massDampingCheckBox_clicked();

    void on_dampingCoeffEdit_editingFinished();

private:
    Controller &cont_;
    Ui::MainWindow *ui;
    bool simRunning_;
    QTimer renderTimer_;

    void setParametersFromUI();
};

#endif // MAINWINDOW_H
