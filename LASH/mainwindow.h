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
    int icall, size_m, total_count, better, worse, ref, con, random;
    vector<float> vazao_calculada;
    vector<float> vazao_observada;
    vector<float> media_vec;
    vector<vector<vector<float>>> matrices;

    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

private slots:

    void on_first_ok_clicked();

    void on_second_ok_clicked();

    void on_third_ok_clicked();

    void on_calculate_clicked();

    void on_actionAjuda_triggered();

    void sceua();

    pair<int, int> best_parameters(int npt, int nopt);

    vector<float> cceua(vector<vector<float>> s, vector<float> sf, vector<float> bl, vector<float> bu, int tot_dias);

    float hydrological_routine(vector<float> x, float tot_dias, bool best_p);

private:
    Ui::MainWindow *ui;
    SecondWindow *second;
    HelpDialog *help;
};

#endif // MAINWINDOW_H
