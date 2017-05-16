#include "results.h"
#include "ui_results.h"
#include <iostream>
results::results(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::results)
{
    ui->setupUi(this);
}

results::~results()
{
    delete ui;
}

void results::on_show_results_clicked()
{
    ui->textBrowser->setText("RESULTADOS AQUI!");

}
