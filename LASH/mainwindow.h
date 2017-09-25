#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QTextStream>

#include "secondwindow.h"
#include "helpdialog.h"

using namespace std;

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

private slots:

    void on_first_ok_clicked();

    void on_second_ok_clicked();

    void on_third_ok_clicked();

    void on_calculate_clicked();

    void on_actionAjuda_triggered();

    void sceua();

    float cceua(vector<vector<float>> s, vector<float> sf, float bl, float bu, int icall, int maxn);

    float hydrological_routine(vector<float> x, float tot_dias);

private:
    Ui::MainWindow *ui;
    SecondWindow *second;
    HelpDialog *help;
};

#endif // MAINWINDOW_H
