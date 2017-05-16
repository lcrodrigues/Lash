#ifndef RESULTS_H
#define RESULTS_H

#include <QWidget>

namespace Ui {
class results;
}

class results : public QWidget
{
    Q_OBJECT

public:
    explicit results(QWidget *parent = 0);
    ~results();

private slots:


    void on_show_results_clicked();

private:
    Ui::results *ui;
};

#endif // RESULTS_H
