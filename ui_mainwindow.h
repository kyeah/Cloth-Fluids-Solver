/********************************************************************************
** Form generated from reading UI file 'mainwindow.ui'
**
** Created by: Qt User Interface Compiler version 4.8.6
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_MAINWINDOW_H
#define UI_MAINWINDOW_H

#include <QtCore/QVariant>
#include <QtGui/QAction>
#include <QtGui/QApplication>
#include <QtGui/QButtonGroup>
#include <QtGui/QCheckBox>
#include <QtGui/QFrame>
#include <QtGui/QGroupBox>
#include <QtGui/QHBoxLayout>
#include <QtGui/QHeaderView>
#include <QtGui/QLabel>
#include <QtGui/QLineEdit>
#include <QtGui/QMainWindow>
#include <QtGui/QMenu>
#include <QtGui/QMenuBar>
#include <QtGui/QPushButton>
#include <QtGui/QStatusBar>
#include <QtGui/QToolBar>
#include <QtGui/QVBoxLayout>
#include <QtGui/QWidget>
#include "glpanel.h"

QT_BEGIN_NAMESPACE

class Ui_MainWindow
{
public:
    QAction *actionExit;
    QAction *actionReset;
    QAction *actionReset_Everything;
    QWidget *centralWidget;
    GLPanel *GLWidget;
    QFrame *parameterFrame;
    QWidget *verticalLayoutWidget;
    QVBoxLayout *verticalLayout;
    QGroupBox *simOptionsBox;
    QWidget *horizontalLayoutWidget;
    QHBoxLayout *horizontalLayout;
    QGroupBox *SimulationBox;
    QPushButton *startSimulationButton;
    QGroupBox *SimParametersBox;
    QLabel *timeStepLabel;
    QLineEdit *timeStepEdit;
    QLabel *penaltyStiffnessLabel;
    QLineEdit *penaltyStiffnessEdit;
    QCheckBox *pinClothCornersCheckbox;
    QGroupBox *activeForcesBox;
    QCheckBox *gravityCheckBox;
    QLabel *gravityGLabel;
    QLineEdit *gravityGEdit;
    QLabel *kDiffLabel;
    QLineEdit *kDiffEdit;
    QCheckBox *fluidDragForceCheckbox;
    QMenuBar *menuBar;
    QMenu *menuFile;
    QMenu *menuScene;
    QToolBar *mainToolBar;
    QStatusBar *statusBar;

    void setupUi(QMainWindow *MainWindow)
    {
        if (MainWindow->objectName().isEmpty())
            MainWindow->setObjectName(QString::fromUtf8("MainWindow"));
        MainWindow->resize(1200, 800);
        actionExit = new QAction(MainWindow);
        actionExit->setObjectName(QString::fromUtf8("actionExit"));
        actionReset = new QAction(MainWindow);
        actionReset->setObjectName(QString::fromUtf8("actionReset"));
        actionReset_Everything = new QAction(MainWindow);
        actionReset_Everything->setObjectName(QString::fromUtf8("actionReset_Everything"));
        centralWidget = new QWidget(MainWindow);
        centralWidget->setObjectName(QString::fromUtf8("centralWidget"));
        GLWidget = new GLPanel(centralWidget);
        GLWidget->setObjectName(QString::fromUtf8("GLWidget"));
        GLWidget->setGeometry(QRect(10, 0, 731, 731));
        QSizePolicy sizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
        sizePolicy.setHorizontalStretch(0);
        sizePolicy.setVerticalStretch(0);
        sizePolicy.setHeightForWidth(GLWidget->sizePolicy().hasHeightForWidth());
        GLWidget->setSizePolicy(sizePolicy);
        GLWidget->setFocusPolicy(Qt::StrongFocus);
        parameterFrame = new QFrame(centralWidget);
        parameterFrame->setObjectName(QString::fromUtf8("parameterFrame"));
        parameterFrame->setGeometry(QRect(749, -1, 441, 731));
        parameterFrame->setFrameShape(QFrame::StyledPanel);
        parameterFrame->setFrameShadow(QFrame::Raised);
        verticalLayoutWidget = new QWidget(parameterFrame);
        verticalLayoutWidget->setObjectName(QString::fromUtf8("verticalLayoutWidget"));
        verticalLayoutWidget->setGeometry(QRect(9, -1, 431, 731));
        verticalLayout = new QVBoxLayout(verticalLayoutWidget);
        verticalLayout->setSpacing(6);
        verticalLayout->setContentsMargins(11, 11, 11, 11);
        verticalLayout->setObjectName(QString::fromUtf8("verticalLayout"));
        verticalLayout->setContentsMargins(0, 0, 0, 0);
        simOptionsBox = new QGroupBox(verticalLayoutWidget);
        simOptionsBox->setObjectName(QString::fromUtf8("simOptionsBox"));
        simOptionsBox->setMaximumSize(QSize(16777215, 220));
        horizontalLayoutWidget = new QWidget(simOptionsBox);
        horizontalLayoutWidget->setObjectName(QString::fromUtf8("horizontalLayoutWidget"));
        horizontalLayoutWidget->setGeometry(QRect(9, 19, 421, 181));
        horizontalLayout = new QHBoxLayout(horizontalLayoutWidget);
        horizontalLayout->setSpacing(6);
        horizontalLayout->setContentsMargins(11, 11, 11, 11);
        horizontalLayout->setObjectName(QString::fromUtf8("horizontalLayout"));
        horizontalLayout->setContentsMargins(0, 0, 0, 0);
        SimulationBox = new QGroupBox(horizontalLayoutWidget);
        SimulationBox->setObjectName(QString::fromUtf8("SimulationBox"));
        startSimulationButton = new QPushButton(SimulationBox);
        startSimulationButton->setObjectName(QString::fromUtf8("startSimulationButton"));
        startSimulationButton->setGeometry(QRect(10, 40, 181, 27));

        horizontalLayout->addWidget(SimulationBox);

        SimParametersBox = new QGroupBox(horizontalLayoutWidget);
        SimParametersBox->setObjectName(QString::fromUtf8("SimParametersBox"));
        timeStepLabel = new QLabel(SimParametersBox);
        timeStepLabel->setObjectName(QString::fromUtf8("timeStepLabel"));
        timeStepLabel->setGeometry(QRect(10, 30, 81, 21));
        timeStepEdit = new QLineEdit(SimParametersBox);
        timeStepEdit->setObjectName(QString::fromUtf8("timeStepEdit"));
        timeStepEdit->setGeometry(QRect(140, 30, 61, 21));
        penaltyStiffnessLabel = new QLabel(SimParametersBox);
        penaltyStiffnessLabel->setObjectName(QString::fromUtf8("penaltyStiffnessLabel"));
        penaltyStiffnessLabel->setGeometry(QRect(10, 50, 131, 21));
        penaltyStiffnessEdit = new QLineEdit(SimParametersBox);
        penaltyStiffnessEdit->setObjectName(QString::fromUtf8("penaltyStiffnessEdit"));
        penaltyStiffnessEdit->setGeometry(QRect(140, 50, 61, 21));
        pinClothCornersCheckbox = new QCheckBox(SimParametersBox);
        pinClothCornersCheckbox->setObjectName(QString::fromUtf8("pinClothCornersCheckbox"));
        pinClothCornersCheckbox->setGeometry(QRect(10, 80, 151, 21));

        horizontalLayout->addWidget(SimParametersBox);


        verticalLayout->addWidget(simOptionsBox);

        activeForcesBox = new QGroupBox(verticalLayoutWidget);
        activeForcesBox->setObjectName(QString::fromUtf8("activeForcesBox"));
        activeForcesBox->setMaximumSize(QSize(16777215, 170));
        gravityCheckBox = new QCheckBox(activeForcesBox);
        gravityCheckBox->setObjectName(QString::fromUtf8("gravityCheckBox"));
        gravityCheckBox->setGeometry(QRect(30, 30, 97, 21));
        gravityGLabel = new QLabel(activeForcesBox);
        gravityGLabel->setObjectName(QString::fromUtf8("gravityGLabel"));
        gravityGLabel->setGeometry(QRect(230, 30, 121, 21));
        gravityGEdit = new QLineEdit(activeForcesBox);
        gravityGEdit->setObjectName(QString::fromUtf8("gravityGEdit"));
        gravityGEdit->setGeometry(QRect(370, 30, 51, 21));
        kDiffLabel = new QLabel(activeForcesBox);
        kDiffLabel->setObjectName(QString::fromUtf8("kDiffLabel"));
        kDiffLabel->setGeometry(QRect(230, 50, 121, 21));
        kDiffEdit = new QLineEdit(activeForcesBox);
        kDiffEdit->setObjectName(QString::fromUtf8("kDiffEdit"));
        kDiffEdit->setGeometry(QRect(370, 50, 51, 21));
        fluidDragForceCheckbox = new QCheckBox(activeForcesBox);
        fluidDragForceCheckbox->setObjectName(QString::fromUtf8("fluidDragForceCheckbox"));
        fluidDragForceCheckbox->setGeometry(QRect(30, 50, 191, 21));

        verticalLayout->addWidget(activeForcesBox);

        MainWindow->setCentralWidget(centralWidget);
        menuBar = new QMenuBar(MainWindow);
        menuBar->setObjectName(QString::fromUtf8("menuBar"));
        menuBar->setGeometry(QRect(0, 0, 1200, 27));
        menuFile = new QMenu(menuBar);
        menuFile->setObjectName(QString::fromUtf8("menuFile"));
        menuScene = new QMenu(menuBar);
        menuScene->setObjectName(QString::fromUtf8("menuScene"));
        MainWindow->setMenuBar(menuBar);
        mainToolBar = new QToolBar(MainWindow);
        mainToolBar->setObjectName(QString::fromUtf8("mainToolBar"));
        MainWindow->addToolBar(Qt::TopToolBarArea, mainToolBar);
        statusBar = new QStatusBar(MainWindow);
        statusBar->setObjectName(QString::fromUtf8("statusBar"));
        MainWindow->setStatusBar(statusBar);

        menuBar->addAction(menuFile->menuAction());
        menuBar->addAction(menuScene->menuAction());
        menuFile->addAction(actionExit);
        menuScene->addAction(actionReset);
        menuScene->addAction(actionReset_Everything);

        retranslateUi(MainWindow);

        QMetaObject::connectSlotsByName(MainWindow);
    } // setupUi

    void retranslateUi(QMainWindow *MainWindow)
    {
        MainWindow->setWindowTitle(QApplication::translate("MainWindow", "Virtual Catwalk", 0, QApplication::UnicodeUTF8));
        actionExit->setText(QApplication::translate("MainWindow", "Exit", 0, QApplication::UnicodeUTF8));
        actionReset->setText(QApplication::translate("MainWindow", "Clear Scene", 0, QApplication::UnicodeUTF8));
        actionReset_Everything->setText(QApplication::translate("MainWindow", "Reset Everything", 0, QApplication::UnicodeUTF8));
        simOptionsBox->setTitle(QApplication::translate("MainWindow", "Simulation Options", 0, QApplication::UnicodeUTF8));
        SimulationBox->setTitle(QApplication::translate("MainWindow", "Simulation Controls", 0, QApplication::UnicodeUTF8));
        startSimulationButton->setText(QApplication::translate("MainWindow", "Start Simulation", 0, QApplication::UnicodeUTF8));
        SimParametersBox->setTitle(QApplication::translate("MainWindow", "Parameters", 0, QApplication::UnicodeUTF8));
        timeStepLabel->setText(QApplication::translate("MainWindow", "Time Step:", 0, QApplication::UnicodeUTF8));
        penaltyStiffnessLabel->setText(QApplication::translate("MainWindow", "Penalty Stiffness:", 0, QApplication::UnicodeUTF8));
        pinClothCornersCheckbox->setText(QApplication::translate("MainWindow", "Pin Cloth Corner", 0, QApplication::UnicodeUTF8));
        activeForcesBox->setTitle(QApplication::translate("MainWindow", "Active Forces", 0, QApplication::UnicodeUTF8));
        gravityCheckBox->setText(QApplication::translate("MainWindow", "Gravity", 0, QApplication::UnicodeUTF8));
        gravityGLabel->setText(QApplication::translate("MainWindow", "Acceleration:", 0, QApplication::UnicodeUTF8));
        kDiffLabel->setText(QApplication::translate("MainWindow", "Diffusion K:", 0, QApplication::UnicodeUTF8));
        fluidDragForceCheckbox->setText(QApplication::translate("MainWindow", "Apply Fluid Drag Force", 0, QApplication::UnicodeUTF8));
        menuFile->setTitle(QApplication::translate("MainWindow", "File", 0, QApplication::UnicodeUTF8));
        menuScene->setTitle(QApplication::translate("MainWindow", "Scene", 0, QApplication::UnicodeUTF8));
    } // retranslateUi

};

namespace Ui {
    class MainWindow: public Ui_MainWindow {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_MAINWINDOW_H
