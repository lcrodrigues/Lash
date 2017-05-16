#include "secondwindow.h"
#include "ui_secondwindow.h"
#include <QFile>
#include <QTextStream>
#include <QMessageBox>
SecondWindow::SecondWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::SecondWindow)
{
    ui->setupUi(this);
}

SecondWindow::~SecondWindow()
{
    delete ui;
}

void SecondWindow::on_show_clicked()
{
    QFile file("Resultados_lash.txt");
    if(!file.open(QIODevice::ReadOnly)){
        QMessageBox::information(0, "info", file.errorString());
    }
    QTextStream in(&file);

    ui->result_browser->setText(in.readAll());
}
