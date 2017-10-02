#include <iostream>
#include <fstream>
#include <cmath>
#include <QMessageBox>
#include <QFileDialog>
#include <fstream>
#include <iomanip>
#include <map>
#include <QtWidgets/QApplication>
#include <QtWidgets/QMainWindow>
#include <QtCharts/QChartView>
#include <QtCharts/QLineSeries>
#include <QtCharts/QValueAxis>
#include <QtCharts/QCategoryAxis>

#include "results.h"
#include "secondwindow.h"
#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "varclasses.h"

QT_CHARTS_USE_NAMESPACE
using namespace std;
using namespace QtCharts;
    
    MainWindow::MainWindow(QWidget *parent) : QMainWindow(parent), ui(new Ui::MainWindow)
    {
        ui->setupUi(this);


    }
    
    MainWindow::~MainWindow()
    {
        delete ui;
    }
    
    
    void MainWindow::on_first_ok_clicked()
    {
        //QString filename1=QFileDialog::getOpenFileName(this, tr("Choose 1st file"), "\\", "All files (*.*);; Text Files (*.txt)");
        //QString filename1 = "/home/hidrica/Lash/Planilhas/Dados_Pelotas.txt";

        QString filename1 = "C:/Users/Cliente/Desktop/Lash/Planilhas/Dados_Pelotas.txt";

        //QMessageBox::information(this, tr("Done!"), filename);
        ui->first_file->setText(filename1);
        /*VarMeteorologicas objmet;
        QString arquivo1;
        arquivo1 = ui->first_file->text();
        objmet.file1_name=arquivo1;
        QMessageBox::information(this, "OK!", objmet.file1_name);*/

    }
    
    void MainWindow::on_second_ok_clicked()
    {
        //QString filename2=QFileDialog::getOpenFileName(this, tr("Choose 1st file"), "\\", "All files (*.*);; Text Files (*.txt)");
        //QString filename2 = "/home/hidrica/Lash/Planilhas/Uso_Solo_Pelotas.txt";

        QString filename2 = "C:/Users/Cliente/Desktop/Lash/Planilhas/Uso_Solo_Pelotas.txt";

        //QMessageBox::information(this, tr("Done!"), filename);
        ui->second_file->setText(filename2);
        /*VarUsoDoSolo objusosolo;
        QString arquivo2;
        arquivo2=ui->second_file->text();
        objusosolo.file2_name=arquivo2;
        QMessageBox::information(this, "OK!2", objusosolo.file2_name);*/
    }
    
    void MainWindow::on_third_ok_clicked()
    {
        //QString filename3=QFileDialog::getOpenFileName(this, tr("Choose 1st file"), "\\", "All files (*.*);; Text Files (*.txt)");
        //QString filename3 = "/home/hidrica/Lash/Planilhas/Mapas_Pelotas.txt";

        QString filename3 = "C:/Users/Cliente/Desktop/Lash/Planilhas/Mapas_Pelotas.txt";

        //QMessageBox::information(this, tr("Done!"), filename);
        ui->third_file->setText(filename3);
        /*VarEntrada objentrada;
       QString arquivo3;
       arquivo3=ui->third_file->text();
       objentrada.file3_name=arquivo3;
       QMessageBox::information(this, "OK!3", objentrada.file3_name);*/
    }

    void MainWindow::sceua() {

        icall = 0;
        int ngs = 5;
        int iniflg = 0;

        vector<float> x0;

        x0.push_back(2.9265); //kcr
        x0.push_back(1.8663); //kb
        x0.push_back(137.4891); //kss
        x0.push_back(28.3582); //cs
        x0.push_back(241.1142); //css
        x0.push_back(44.3362); //cb
        x0.push_back(0.0237); //lambda

        // upper / lower
        vector<float> bu;
        bu.push_back(5);
        bu.push_back(6);
        bu.push_back(250);
        bu.push_back(100);
        bu.push_back(500);
        bu.push_back(100);
        bu.push_back(0.5);

        vector<float> bl;
        for(int i = 0; i < 6; i++)
            bl.push_back(0.1);
        bl.push_back(0.01);

        // ---

        int maxn = 1000;

        /*if(ui->checkBox_parametros->isChecked())    //kcr
            x0.push_back(2.0213);

        else
            x0.push_back(ui->lineEdit_kcr->text().toFloat(NULL));


        if(ui->checkBox_parametros->isChecked())    //kb
            x0.push_back(0.6338);
        else
            x0.push_back(ui->lineEdit_kb->text().toFloat(NULL));


        if(ui->checkBox_parametros->isChecked())    //kss
            x0.push_back(5.2609);
        else
            x0.push_back(ui->lineEdit_kss->text().toFloat(NULL));


        if(ui->checkBox_parametros->isChecked())    //cs
            x0.push_back(46.5581);
        else
            x0.push_back(ui->lineEdit_cs->text().toFloat(NULL));


        if(ui->checkBox_parametros->isChecked())    //css
            x0.push_back(165.4442);
        else
            x0.push_back(ui->lineEdit_css->text().toFloat(NULL));


        if(ui->checkBox_parametros->isChecked())    //cb
            x0.push_back(57.5433);
        else
            x0.push_back(ui->lineEdit_cb->text().toFloat(NULL));


        if(ui->checkBox_parametros->isChecked())    //lambda
            x0.push_back(0.1145);
        else
            x0.push_back(ui->lineEdit_coef->text().toFloat(NULL));*/

        //Initialize SCE parameters:

        int nopt = x0.size();
        int npg = 2 * nopt + 1;
        int nps = nopt + 1;
        int nspl = npg;
        int npt = npg * ngs;
        float bestf, worstf, tot_dias;
        tot_dias = ui->num_dias->text().toInt(NULL);

        vector<float> bound, bestx, worstx, xf;

        for(uint i = 0; i < bu.size(); i++)
            bound.push_back(bu.at(i) - bl.at(i));

        vector<vector<float>> x;

        for(int i = 0; i < npt; i++) {
            vector<float> row;
            for(int j = 0; j < nopt; j++)
                row.push_back(0);

            x.push_back(row);
        }

        for(int i = 0; i < npt; i++) {  //inicia x0
            for(int j = 0; j < nopt; j++) {
                float r = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
                x.at(i).at(j) = bl.at(j) + r * bound.at(j);
            }
        }

        if (iniflg == 1) {
            for(uint i = 0; i < x.at(0).size(); i++)
                x.at(0).at(i) = x0.at(i);
        }

        // código real
        for(int i = 0; i < npt; i++) {
            xf.push_back(hydrological_routine(x.at(i), tot_dias, false));
            icall++;
            cout << "RMSE: " << xf.at(i) << endl;
        }

    /*
        //código para fins acadêmicos
        xf.push_back(1.48814);
        xf.push_back(1.48627);
        xf.push_back(1.70392);
        xf.push_back(1.56368);
        xf.push_back(1.47106);
        xf.push_back(2.78047);
        xf.push_back(2.34897);
        xf.push_back(2.57968);
        xf.push_back(1.71421);
        xf.push_back(1.52463);
        xf.push_back(1.84637);
        xf.push_back(1.58905);
        xf.push_back(1.94826);
        xf.push_back(2.47317);
        xf.push_back(1.9118);
        xf.push_back(1.80122);
        xf.push_back(2.90578);
        xf.push_back(2.13002);
        xf.push_back(1.71444);
        xf.push_back(1.79674);
        xf.push_back(3.28654);
        xf.push_back(1.4228);
        xf.push_back(2.69704);
        xf.push_back(2.59931);
        xf.push_back(14.9418);
        xf.push_back(1.70445);
        xf.push_back(1.51432);
        xf.push_back(2.34553);
        xf.push_back(2.57833);
        xf.push_back(2.00606);
        xf.push_back(1.4754);
        xf.push_back(1.80335);
        xf.push_back(2.27863);
        xf.push_back(2.13251);
        xf.push_back(1.6062);
        xf.push_back(3.13469);
        xf.push_back(1.49792);
        xf.push_back(1.70897);
        xf.push_back(1.8492);
        xf.push_back(2.61157);
        xf.push_back(2.55223);
        xf.push_back(2.17108);
        xf.push_back(1.55858);
        xf.push_back(1.48856);
        xf.push_back(2.13879);
        xf.push_back(2.32624);
        xf.push_back(1.63237);
        xf.push_back(1.53213);
        xf.push_back(1.90212);
        xf.push_back(1.71208);
        xf.push_back(1.76119);
        xf.push_back(2.29582);
        xf.push_back(1.81799);
        xf.push_back(3.75195);
        xf.push_back(3.18964);
        xf.push_back(2.84476);
        xf.push_back(2.91534);
        xf.push_back(1.52061);
        xf.push_back(1.49778);
        xf.push_back(1.47806);
        xf.push_back(1.65241);
        xf.push_back(3.17681);
        xf.push_back(2.10404);
        xf.push_back(1.80703);
        xf.push_back(1.54959);
        xf.push_back(2.0237);
        xf.push_back(1.42918);
        xf.push_back(2.99572);
        xf.push_back(3.30334);
        xf.push_back(1.63166);
        xf.push_back(1.83934);
        xf.push_back(1.49267);
        xf.push_back(3.9123); //fake
        xf.push_back(2.8526);
        xf.push_back(1.71648);
        //fim do código para fins acadêmicos    */

        vector<pair<float, vector<float>>> xf_to_f;

        //liga x com xf
        for(int i = 0; i < npt; i++)
            xf_to_f.push_back(make_pair(xf.at(i), x.at(i)));

        sort(xf_to_f.begin(), xf_to_f.end());

        vector<vector<float>> x_ordered;

        //cria matriz x ordenada
        for(int i = 0; i < npt; i++)
            x_ordered.push_back(xf_to_f.at(i).second);

        bestx = x_ordered.at(0);    //bestx recebe melhores parâmetros
        bestf = xf_to_f.at(0).first;   //bestf recebe melhor xf

        worstx = x_ordered.at(npt - 1); //worstx recebe piores parâmetros
        worstf = xf_to_f.at(npt - 1).first; //worstf recebe pior xf

        //gerar gráfico
        /*
        hydrological_routine(bestx, tot_dias, true);

        QChart *chart = new QChart();
        chart->legend()->hide();
        chart->setTitle("Valores observados e calculados");

        QValueAxis *axisX = new QValueAxis;
        chart->addAxis(axisX, Qt::AlignBottom);

        QLineSeries *series = new QLineSeries;
        for(uint i = 0; i < vazao_calculada.size(); i++)
                    series->append(i + 1, vazao_calculada.at(i));
        chart->addSeries(series);

        QValueAxis *axisY = new QValueAxis;
        axisY->setLinePenColor(series->pen().color());

        chart->addAxis(axisY, Qt::AlignLeft);
        series->attachAxis(axisX);
        series->attachAxis(axisY);

        series = new QLineSeries;
        for(uint i = 0; i < vazao_observada.size(); i++) {
            if(vazao_observada.at(i) != -9999)
                series->append(i + 1, vazao_observada.at(i));
        }
        chart->addSeries(series);

        series->attachAxis(axisX);

        QChartView *chartView = new QChartView(chart);
        chartView->setRenderHint(QPainter::Antialiasing);

        second = new SecondWindow(this);
        second->setCentralWidget(chartView);
        second->show(); */


        //exibir dados
        cout << "The initial loop: 0" << endl;
        cout << "BESTF: " << bestf << endl;
        cout << "BESTX: ";
        for(uint i = 0; i < bestx.size(); i++)
            cout << bestx.at(i) << " ";
        cout << "\nWORSTF: " << worstf << endl;
        cout << "WORSTX: ";
        for(uint i = 0; i < worstx.size(); i++)
            cout << worstx.at(i) << " ";
        cout << endl;

        int nloop = 0;
        int size_m = 0;

        cout << "Entra em while(maxn)\n";
        while(icall < maxn) { //loop contando iterações
            cout << "ICALL: " << icall << endl;
            nloop++;
            cout << "For de 0 a 5\n";
            for(int igs = 0; igs < ngs; igs++) {    //0 a 5
                vector<vector<float>> cx(npg, vector<float>(nopt));
                vector<float> cf(npg);

                for(int k1 = 0; k1 < npg; k1++) {   //0 a 15
                    int k2 = k1 * ngs + igs;

                    cf.at(k1) = xf.at(k2);    //cf recebe xf embaralhado
                    cx.at(k1) = x_ordered.at(k2); //cx recebe x embaralhado
                }
                cout << "Embaralhou intervalo " << igs + 1 << endl;
            
                cout << "Entra em loop de 0 a 15\n";
                for(int i = 0; i < nspl; i++) { //0 a 15
                    vector<int> lcs(8);
                    lcs.at(0) = 0;

                    int size = 14;
                    vector<int> aux_p(size);
                    for(uint l = 0; l < aux_p.size(); l++)
                        aux_p.at(l) = l + 1;

                    for(int j = 1; j < nps; j++) {

                        int r = rand() % size;
                        size--;

                        lcs.at(j) = aux_p.at(r);
                        aux_p.erase(aux_p.begin() + r);

                        /*for(int k = 0; k < 1000; k++) {

                            float r = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
                            int p = 1 + floor(npg + 0.5 - sqrt(pow((npg + 0.5), 2) - npg * (npg + 1) * r));

                            for(int l = 0; l < lcs.size(); l++) {
                                if(p == lcs.at(l)) {
                                    auxb = false; //se já possui a posição, sai pra gerar uma nova
                                    break;
                                }
                            }

                            if(auxb) {  //se não possui a posição, adiciona em lcs
                                lcs.at(j) = p;
                                break;
                            }
                        } */
                    }

                    sort(lcs.begin(), lcs.end());   //ordena lcs

                    vector<vector<float>> s(nps, vector<float>(nopt));
                    vector<float> snew;
                    vector<float> sf(nps);

                    for(int j = 0; j < nps; j++) {
                        s.at(j) = cx.at(lcs.at(j));
                        sf.at(j) = cf.at(lcs.at(j));
                     }

                    //chamar cceua
                    cout << "chamando cceua\n";
                    float cce = cceua(&snew, s, sf, bl, bu, tot_dias);
                    cout << "fim do cceua\n";

                    s.at(nps - 1) = snew;
                    sf.at(nps - 1) = cce;

                    //põe o simplex no complex
                    for(int j = 0; j < nps; j++) {
                        cx.at(lcs.at(j)) = s.at(j);
                        cf.at(lcs.at(j)) = sf.at(j); 
                    }
                }

                //põe o complex na população
                for(int k1 = 0; k1 < npg; k1++) {   //0 a 15
                    int k2 = k1 * ngs + igs;

                    x_ordered.at(k2) = cx.at(k1);
                    xf.at(k2) = cf.at(k1);
                }
            }

            //ordena x de acordo com xf
            vector<pair<float, vector<float>>> new_xf_to_f;
            for(uint j = 0; j < xf.size(); j++)
                new_xf_to_f.push_back(make_pair(xf.at(j), x_ordered.at(j)));

            sort(new_xf_to_f.begin(), new_xf_to_f.end());
            for(int j = 0; j < npt; j++)
                x_ordered.at(j) = new_xf_to_f.at(j).second;

            //melhor x e xf
            bestx = x_ordered.at(0);
            bestf = new_xf_to_f.at(0).first;

            //pior x e xf
            worstx = x_ordered.at(npt - 1);
            worstf = new_xf_to_f.at(npt - 1).first;

            cout << "Evolution loop: " << nloop << " - Trial - " << icall << endl;
            cout << "BESTF: " << bestf << endl;
            cout << "BESTX: ";
            for(uint i = 0; i < bestx.size(); i++)
                cout << bestx.at(i) << " ";
            cout << "\nWORSTF: " << worstf << endl;
            cout << "WORSTX: ";
            for(uint i = 0; i < worstx.size(); i++)
                cout << worstx.at(i) << " ";
            cout << endl;

            vector<vector<float>> f_matrix(npt, vector<float>(nopt + 1));
            for(int j = 0; j < npt; j++) {
                f_matrix.at(j).push_back(j + 1);

                for(int k = 0; k < nopt; k++)
                    f_matrix.at(j).push_back(new_xf_to_f.at(j).second.at(k));

                f_matrix.at(j).push_back(new_xf_to_f.at(j).first);
            }

            matrices.push_back(f_matrix);
            size_m++;
        }

        cout << "\n---------------------------------------------------\n\n";

        for(int i = 0; i < size_m; i++) {
            for(int j = 0; j < npt; j++) {
                for(int k = 0; k < nopt; k++) {
                    cout << matrices.at(i).at(j).at(k) << " ";
                }
                cout << endl;
            }
            cout << "\n---------------------------------------------------\n\n";
        }
    }

    float MainWindow::cceua(vector<float> *snew, vector<vector<float>> s, vector<float> sf, vector<float> bl, vector<float> bu, int tot_dias) {
        vector<float> sb = s.at(0); //sb melhor s
        vector<float> sw = s.at(6); //sw pior s

        float fw = sf.at(sf.size() - 1);    //fw pior sf
        float alpha = 1.0, beta = 0.5;

        vector<float> ce;
        for(uint j = 0; j < sf.size() - 1; j++) { //até 8 - 1
            float sum = 0;

            for(uint k = 0; k < s.at(0).size(); k++) //até 7
                sum += s.at(j).at(k);

            float media = sum / (7 - 1); //media da linha de s, excluindo o pior ponto
                                         //7 é o equivalente a nopt
            ce.push_back(media);

            float svalue = ce.at(j) + alpha * (ce.at(j) - sw.at(j));  //calcula valor pro snew
            snew->push_back(svalue);
        }

        bool ibound = false;
        for(uint j = 0; j < snew->size(); j++) { //verifica os limites
            float s1, s2;
            s1 = snew->at(j) - bl.at(j);
            s2 = bu.at(j) - snew->at(j);

            if(s1 < 0 || s2 < 0) {
                ibound = true;
                break;
            }
        }

        if(ibound) {   //recalcula caso fora de limites
            for(uint j = 0; j < snew->size(); j++)
                snew->at(j) = bl.at(j) + ((1 + (rand() % 7 + 1)) * (bu.at(j) - bl.at(j)));  //7 -> nopt
        }

        float fnew;
        cout << "HR reflection" << endl;
        fnew = hydrological_routine(*snew, tot_dias, false);
        icall++;

        if(fnew > fw) { //reflexão falhou, tentando ponto de contração
            for(uint j = 0; j < snew->size(); j++)
                snew->at(j) = sw.at(j) + beta * (ce.at(j) - sw.at(j));

            cout << "HR contraction (reflection failed)" << endl;
            fnew = hydrological_routine(*snew, tot_dias, false);
            icall++;

            if(fnew > fw) { //reflexão e contração falharam, tentando ponto aleatório
                for(uint j = 0; j < snew->size(); j++)
                    snew->at(j) = bl.at(j) + ((rand() % 7 + 1) * (bu.at(j) - bl.at(j)));  //7 -> nopt

                cout << "HR random (reflection and contraction failed)" << endl;
                fnew = hydrological_routine(*snew, tot_dias, false);
                icall++;
            }
        }

        return fnew;
    }


    float MainWindow::hydrological_routine(vector<float> x, float tot_dias, bool best_p) {
        int i, j, dia, sub_b=0;
        float kcr, Kb, Kss, Cs, Css, Cb, Coef_Ia, evap_lamina, rmse = 0;
        bool flag = false;

        kcr = x.at(0);
        Kb = x.at(1);
        Kss = x.at(2);
        Cs = x.at(3);
        Css = x.at(4);
        Cb = x.at(5);
        Coef_Ia = x.at(6);

        // float p5tst;
        float h_sistrad2;

        VarDiaAnterior anterior;
        SubBacia subw;

        VarMeteorologicas met;
        VarUsoDoSolo usosolo;
        VarEntrada entrada;

        VarEvapotranspiracao evapo;
        VarInterceptacao intercep;
        VarSolo solo;
        VarESD ESD;
        VarESS ESS;
        VarEB EB;

        //Método que seta o número de sub-bacias.
        subw.setNumSub_b(ui->third_file->text());

        //Armazena arquivos de entrada em memória para facilitar acesso aos dados.
        met.criaVetor1(subw, tot_dias);
        met.loadData1(subw, ui->first_file->text(), tot_dias);
        usosolo.criaVetor2(subw, tot_dias);
        usosolo.loadData2(subw, ui->second_file->text(), tot_dias);
        entrada.criaVetor3(subw);
        entrada.loadData3(subw, ui->third_file->text());
        anterior.iniciaVetores(subw);

    //======================ABERTURA DE ARQUIVOS DE SAÍDA=============================

                //Abre arquivo de saída "Resultados" (Que será apresentado ao final da execução).
       /* std::ofstream outFile("Resultados_lash.txt", std::ios_base::out);
        if(!outFile ) { // Abertura falhou ...
            std::cout << "outFile não pode ser aberto para saída\n";
            QMessageBox::information(this, tr("info"), "Nao pode abrir arquivo!");
            exit(-1);
        }
        outFile << "RESULTADOS LASH PARA " << tot_dias << " DIAS." << std::endl << std::endl;

                //Abre arquivo de saída "Lash_ETc".
        //if(ui->checkBox_ETc->isChecked()){
        std::ofstream fileETc("Lash_ETc.txt", std::ios_base::out);
        if(!fileETc ) { // Abertura falhou ...
            std::cout << "Arquivo de ETc não pode ser aberto para saída\n";
            QMessageBox::information(this, tr("info"), "Nao pode abrir arquivo!");
            exit(-1);
        }
        fileETc << "RESULTADOS LASH PARA ETc, DURANTE " << tot_dias << " DIAS." << std::endl << std::endl;
        //}

                //Abre arquivo de saída "Lash_ETr".
        //if(ui->checkBox_ETr->isChecked()){
        std::ofstream fileETr("Lash_ETr.txt", std::ios_base::out);
        if(!fileETr ) { // Abertura falhou ...
            std::cout << "Arquivo de ETr não pode ser aberto para saída\n";
            QMessageBox::information(this, tr("info"), "Nao pode abrir arquivo!");
            exit(-1);
        }
        //fileETr << "RESULTADOS LASH PARA ETr, DURANTE " << tot_dias << " DIAS." << std::endl << std::endl;
        //}

                //Abre arquivo de saída "Lash_Vazao".
        //if(ui->checkBox_Vazao->isChecked()){
        std::ofstream fileVazao("Lash_Vazão.txt", std::ios_base::out);
        if(!fileVazao ) { // Abertura falhou ...
            std::cout << "Arquivo de Vazao não pode ser aberto para saída\n";
            QMessageBox::information(this, tr("info"), "Nao pode abrir arquivo!");
            exit(-1);
        }
        //fileVazao << "RESULTADOS LASH PARA VAZÃO, DURANTE " << tot_dias << " DIAS." << std::endl << std::endl;
        //}

                //Abre arquivo de saída "PE".
        //if(ui->checkBox_X->isChecked()){
        std::ofstream filePe("Lash_Pe.txt", std::ios_base::out);
        if(!filePe) { // Abertura falhou ...
            std::cout << "Arquivo de P.E não pode ser aberto para saída\n";
            QMessageBox::information(this, tr("info"), "Nao pode abrir arquivo!");
            exit(-1);
        }
        filePe << "RESULTADOS LASH PARA PREC. EFETIVA, DURANTE " << tot_dias << " DIAS." << std::endl << std::endl;
        //}

                //Abre arquivo de saída "PTS".
        //if(ui->checkBox_X->isChecked()){
        std::ofstream filePts("Lash_Pts.txt", std::ios_base::out);
        if(!filePts) { // Abertura falhou ...
            std::cout << "Arquivo de PTS não pode ser aberto para saída\n";
            QMessageBox::information(this, tr("info"), "Nao pode abrir arquivo!");
            exit(-1);
        }
        filePts << "RESULTADOS LASH PARA PREC. TOTAL, DURANTE " << tot_dias << " DIAS." << std::endl << std::endl;
        //}

                //Abre arquivo de saída "DCR".
        //if(ui->checkBox_X->isChecked()){
        std::ofstream fileDcr("Lash_DCR.txt", std::ios_base::out);
        if(!fileDcr) { // Abertura falhou ...
            std::cout << "Arquivo de DCR não pode ser aberto para saída\n";
            QMessageBox::information(this, tr("info"), "Nao pode abrir arquivo!");
            exit(-1);
        }
        fileDcr << "RESULTADOS LASH PARA ASCENSÃO CAPILAR, DURANTE " << tot_dias << " DIAS." << std::endl << std::endl;
        //}

                //Abre arquivo de saída "Ia".
        //if(ui->checkBox_X->isChecked()){
        std::ofstream fileIa("Lash_Ia.txt", std::ios_base::out);
        if(!fileIa) { // Abertura falhou ...
            std::cout << "Arquivo de IA não pode ser aberto para saída\n";
            QMessageBox::information(this, tr("info"), "Nao pode abrir arquivo!");
            exit(-1);
        }
        fileIa << "RESULTADOS LASH PARA ABSTRAÇÃO INICIAL, DURANTE " << tot_dias << " DIAS." << std::endl << std::endl;
        //}

                //Abre arquivo de saída "At".
        //if(ui->checkBox_X->isChecked()){
        std::ofstream fileAt("Lash_At.txt", std::ios_base::out);
        if(!fileAt) { // Abertura falhou ...
            std::cout << "Arquivo de AT não pode ser aberto para saída\n";
            QMessageBox::information(this, tr("info"), "Nao pode abrir arquivo!");
            exit(-1);
        }
        //fileAt << "RESULTADOS LASH PARA ARMAZENAMENTO TOTAL, DURANTE " << tot_dias << " DIAS." << std::endl << std::endl;
        //}

                //Abre arquivo de saída "Vazao_ESD".
        //if(ui->checkBox_X->isChecked()){
        std::ofstream fileVESD("Lash_VazaoESD.txt", std::ios_base::out);
        if(!fileVESD) { // Abertura falhou ...
            std::cout << "Arquivo de VAZAO ESD não pode ser aberto para saída\n";
            QMessageBox::information(this, tr("info"), "Nao pode abrir arquivo!");
            exit(-1);
        }
        fileVESD << "RESULTADOS LASH PARA VAZAO ESD, DURANTE " << tot_dias << " DIAS." << std::endl << std::endl;
        //}

                //Abre arquivo de saída "Vazao_ESS".
        //if(ui->checkBox_X->isChecked()){
        std::ofstream fileVESS("Lash_VazaoESS.txt", std::ios_base::out);
        if(!fileVESS){ // Abertura falhou ...
            std::cout << "Arquivo de VAZAO ESS não pode ser aberto para saída\n";
            QMessageBox::information(this, tr("info"), "Nao pode abrir arquivo!");
            exit(-1);
        }
        fileVESS << "RESULTADOS LASH PARA VAZAO ESS, DURANTE " << tot_dias << " DIAS." << std::endl << std::endl;
        //}

                //Abre arquivo de saída "Vazao_EB".
        //if(ui->checkBox_X->isChecked()){
        std::ofstream fileVEB("Lash_VazaoEB.txt", std::ios_base::out);
        if(!fileVEB) { // Abertura falhou ...
            std::cout << "Arquivo de VAZAO EB não pode ser aberto para saída\n";
            QMessageBox::information(this, tr("info"), "Nao pode abrir arquivo!");
            exit(-1);
        }
        fileVEB << "RESULTADOS LASH PARA VAZAO EB, DURANTE " << tot_dias << " DIAS." << std::endl << std::endl;
        //}

               //Abre arquivo de saida teste "Lash_TESTE"
        std::ofstream fileTest("Lash_TESTE.txt", std::ios_base::out);
        if(!fileTest) { // Abertura falhou ...
            std::cout << "Arquivo de TESTE não pode ser aberto para saída\n";
            QMessageBox::information(this, tr("info"), "Nao pode abrir arquivo!");
            exit(-1);
        }   */

    //========================================================================================================

        vector<float> vazao_total;
        float total_dias_observados = tot_dias;
        float dif_valores = 0;

        for(i=0; i<tot_dias; i++){
            float vazao_diaria = 0;
            dia=i+1;

            met.setVarMet(dia, sub_b, subw); //função chamada para ter a informação da data.

 /*           fileETc << std::defaultfloat << met.dia << "/" << met.mes << "/" << met.ano;
            //fileETr << std::defaultfloat << met.dia << "/" << met.mes << "/" << met.ano;
            //fileVazao << std::defaultfloat <<met.dia << "/" << met.mes << "/" << met.ano;
            filePe << std::defaultfloat << met.dia << "/" << met.mes << "/" << met.ano;
            filePts << std::defaultfloat << met.dia << "/" << met.mes << "/" << met.ano;
            fileDcr << std::defaultfloat << met.dia << "/" << met.mes << "/" << met.ano;
            fileIa << std::defaultfloat << met.dia << "/" << met.mes << "/" << met.ano;
            //fileAt << std::defaultfloat << met.dia << "/" << met.mes << "/" << met.ano;
            fileVESD << std::defaultfloat << met.dia << "/" << met.mes << "/" << met.ano;
            fileVESS << std::defaultfloat << met.dia << "/" << met.mes << "/" << met.ano;
            fileVEB << std::defaultfloat << met.dia << "/" << met.mes << "/" << met.ano;    */

            for(j=0; j<(subw.numSub_b); j++){
                float vazao_sub = 0;
                sub_b=j;
                //Busca na matriz as váriaveis do arquivo de entrada.
                met.setVarMet(dia, sub_b, subw);
                usosolo.setVarUS(dia, sub_b, subw);
                h_sistrad2=usosolo.h_sistrad;
                entrada.setVarEntrada(sub_b, h_sistrad2);

                //Layout do arquivo de saída "Resultados".
                //outFile << "Dia " << dia << ", " << "Sub-bacia " << sub_b+1 << ":" << std::endl;
                //std::cout << "Dia " << dia << "," <<" Sub_b " << sub_b+1 << ":" <<std::endl;

                //SETA TODOS OS ATRIBUTOS DA CLASSE EVAPOTRANSPIRAÇÃO
                evapo.setAll(dia, &evap_lamina, met, entrada, usosolo);

                /*
                evapo.setTmedia(met);
                evapo.setEs(met);
                evapo.setEa(met);
                evapo.setSsvpc();
                evapo.setAp(entrada);
                evapo.setPc();
                evapo.setRns(usosolo, met);
                evapo.setSbc();
                evapo.setJday(dia);
                evapo.setDr();
                evapo.setSolar_declination();
                evapo.setWs(entrada);
                evapo.setRa(entrada);
                evapo.setRso(entrada);
                evapo.setRnl(met);
                evapo.setRn();
                evapo.setTkv();
                evapo.setU10(usosolo, met);
                evapo.setAeroResist(usosolo);
                evapo.setETc(usosolo);
                evapo.setOthers(&evap_lamina);
                */


                //SETA TODOS OS ATRIBUTOS DA CLASSE INTERCEPTAÇÃO
                intercep.setAll(dia, sub_b, evap_lamina, anterior, met, usosolo);
                /*intercep.setCri(usosolo);
                intercep.setLam1(dia, sub_b, anterior);
                intercep.setLam2(met);
                intercep.setEvapLam(evap_lamina);
                intercep.setLam3();
                */


                solo.setAll(dia, sub_b, Coef_Ia, kcr, anterior, entrada, intercep, met);
                evapo.setETr(solo);
                /*solo.setAm(entrada); //BLOCO COM ERROS!!!!(CHAMADAS ERRADAS: LAM ess, eb)
                solo.setAt(dia, sub_b, anterior);
                solo.setPts(intercep, met);
                solo.setS();
                solo.setM(sub_b, anterior);
                if(ui->checkBox_parametros->isChecked()){
                    Coef_Ia=0.1145;
                }
                else{
                    Coef_Ia=ui->lineEdit_coef->text().toFloat(NULL);
                }
                solo.setIa(Coef_Ia);0.552709 0.873929 1.10771 1.29059 1.4442 1.57715 1.69497
                solo.setPe();
                solo.setAcc();
                 if(ui->checkBox_parametros->isChecked()){
                    Kss=5.2609;
                }
                else{
                    Kss=ui->lineEdit_kss->text().toFloat(NULL);
                }
                ESS.setLamESS(solo, Kss);
                solo.setAc();
                if(ui->checkBox_parametros->isChecked()){
                    Kb=0.6338;
                }
                else{
                    Kb=ui->lineEdit_kb->text().toFloat(NULL);
                }
                EB.setLamEB(solo, Kb);
                solo.setAcr();
                if(ui->checkBox_parametros->isChecked()){
                    kcr=2.0213;
                }
                else{
                    kcr=ui->lineEdit_kcr->text().toFloat(NULL);
                }
                solo.setDcr(kcr);
                solo.setAl();
                solo.setApmp();
                solo.setKs();
                evapo.setETr(solo);
                */

                //SETA TODOS OS ATRIBUTOS DA CLASSE INTERCEPTAÇÃO
                intercep.setAll(dia, sub_b, evap_lamina, anterior, met, usosolo);
                /*
                intercep.setCri(usosolo);
                intercep.setLam1(dia, sub_b, anterior);
                intercep.setLam2(met);
                intercep.setEvapLam(evap_lamina);
                intercep.setLam3();*/

                //SETA TODOS OS ATRIBUTOS DA CLASSE ESD

                ESD.setAll(dia, sub_b, Cs, anterior, solo, entrada);

                //SETA TODOS OS ATRIBUTOS DA CLASSE ESS

                ESS.setAll(dia, sub_b, Kss, Css, solo, anterior, entrada);

                //SETA TODOS OS ATRIBUTOS DA CLASSE EB

                EB.setAll(dia, sub_b, Kb, Cb, solo, anterior, entrada);


                //Seta a variável "vazao_total".
                vazao_sub = ESD.vazao_ESD + ESS.vazao_ESS + EB.vazao_EB;
                vazao_diaria += vazao_sub;

                //cout << "sub bacia " << sub_b + 1 << ", vazão sub bacia: " << vazao_sub << endl;

    //===================Escreve os resultados desejados em seus respectivos arquivos de saída.=================

                /*if(ui->checkBox_ETc->isChecked()){
                    std::cout << "Etc: " << evapo.ETc << std::endl;
                    fileETc << "\t" << std::fixed << std::setprecision(4) << evapo.ETc;
                }
                if(ui->checkBox_ETr->isChecked()){
                    std::cout << "ETr: " << evapo.ETr << std::endl;
                    //fileETr << "\t" << std::fixed << std::setprecision(4) << evapo.ETr;
                    fileETr << std::fixed << std::setprecision(2) << evapo.ETr << "\t";
                }
                if(ui->checkBox_Vazao->isChecked()){
                    std::cout << "Vazao: " << vazao_total <<std::endl;
                    //fileVazao << "\t" << std::fixed << std::setprecision(4) << vazao_total;
                    fileVazao << std::fixed << std::setprecision(4) << vazao_total << "\t";
                }
                if(ui->checkBox_Pe->isChecked()){
                    filePe << "\t" << std::fixed << std::setprecision(4) << solo.pe;
                }
                if(ui->checkBox_Pts->isChecked()){
                    filePts << "\t" << std::fixed << std::setprecision(4) << solo.pts;
                }
                if(ui->checkBox_Dcr->isChecked()){
                    fileDcr << "\t" << std::fixed << std::setprecision(4) << solo.dcr;
                }
                if(ui->checkBox_Ia->isChecked()){
                    fileIa << "\t" << std::fixed << std::setprecision(4) << solo.ia;
                }
                if(ui->checkBox_At->isChecked()){
                    //fileAt << "\t" << std::fixed << std::setprecision(4) << solo.at;
                    fileAt << std::fixed << std::setprecision(2) << solo.at << "\t";
                }
                if(ui->checkBox_vESD->isChecked()){
                    fileVESD << "\t" << std::fixed << std::setprecision(4) << ESD.vazao_ESD;
                }
                if(ui->checkBox_vESS->isChecked()){
                    fileVESS << "\t" << std::fixed << std::setprecision(4) << ESS.vazao_ESS;
                }
                if(ui->checkBox_vEB->isChecked()){
                    fileVEB << "\t" << std::fixed << std::setprecision(4) << EB.vazao_EB;
                }

                outFile << std::endl;
                if(ui->checkBox_ETc->isChecked()){
                    std::cout << "Etc: " << evapo.ETc << std::endl;
                    outFile<< "ETc: " << evapo.ETc << std::endl;
                }
                if(ui->checkBox_ETr->isChecked()){
                    std::cout << "ETr: " << evapo.ETr << std::endl;
                    outFile << "ETr: " << evapo.ETr << std::endl;
                }
                if(ui->checkBox_Vazao->isChecked()){
                    std::cout << "Vazao: " << vazao_total <<std::endl;
                    outFile << "Vazao: " << vazao_total << std::endl;
                }
                outFile << std::endl;
                if(sub_b == 22){
                    fileTest << std::fixed << std::setprecision(2) << solo.at << std::endl;
                }*/


                //AJUSTA AS VARIÁVEIS REFERENTES AO DIA ANTERIOR:
                anterior.SetPrec(sub_b, met.precMedia);
                anterior.SetL3(sub_b, intercep.lam3);
                anterior.SetAm(sub_b, solo.am);
                anterior.SetAt(sub_b, solo.at);
                anterior.SetPts(sub_b, solo.pts);
                anterior.SetPe(sub_b, solo.pe);
                anterior.SetLam_ESS(sub_b, ESS.lam_ESS);
                anterior.SetLam_EB(sub_b, EB.lam_EB);
                anterior.SetDCR(sub_b, solo.dcr);
                anterior.SetETr(sub_b, evapo.ETr);
                anterior.SetVfinal_ESD(sub_b, ESD.vol_final_ESD);
                anterior.SetVfinal_ESS(sub_b, ESS.vol_final_ESS);
                anterior.SetVfinal_EB(sub_b, EB.vol_final_EB);

            }
            //fileTest << std::endl;
         /*   fileETc << std::endl;
            fileETr << std::endl;
            fileVazao << std::endl;
            filePe << std::endl;
            filePts << std::endl;
            fileDcr << std::endl;
            fileIa << std::endl;
            fileAt << std::endl;
            fileVESD << std::endl;
            fileVESS << std::endl;
            fileVEB << std::endl;   */

            if(met.dado_observado == -9999){
                total_dias_observados--;
                flag = false;
            }
            else{
                dif_valores += pow((met.dado_observado - vazao_diaria), 2);
                flag = true;
            }

            if(best_p) {
                vazao_observada.push_back(met.dado_observado);
                vazao_calculada.push_back(vazao_diaria);
            }

            vazao_total.push_back(vazao_diaria);

        }

        if(flag)
             rmse = sqrt(dif_valores/total_dias_observados);

        //Depois de terminar a execução ele abre a segunda tela para mostrar os resultados.
        //second = new SecondWindow(this); //secondwindow é a classe; second é o ponteiro declarado antes;
        //second->show();

    //======================FECHA PONTEIROS USADOS============================================================
        //if(ui->checkBox_ETc->isChecked()){
     /*       fileETc.close();
        //}
        //if(ui->checkBox_ETr->isChecked()){
            fileETr.close();
        //}
        //if(ui->checkBox_Vazao->isChecked()){
            fileVazao.close();
        //}
        filePe.close();
        filePts.close();
        fileDcr.close();
        fileIa.close();
        fileAt.close();
        fileVESD.close();
        fileVESS.close();
        fileVEB.close();
        outFile.close();
        fileTest.close();   */

        return rmse;
    //=========================================================================================================
    }

    void MainWindow::on_calculate_clicked() {
        sceua();
    }
    
    void MainWindow::on_actionAjuda_triggered()
    {
        help = new HelpDialog(this); //secondwindow é a classe; second é o ponteiro declaradoantes;
        help->show();
    }
    
