#include <iostream>
#include <fstream>
#include <ctime>
#include <cmath>
#include <chrono>
#include <QMessageBox>
#include <QFileDialog>
#include <iomanip>
#include <map>
#include <algorithm>
#include <bits/stdc++.h>

#include <omp.h>

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
using namespace std::chrono; // nanoseconds, system_clock, seconds
    
    MainWindow::MainWindow(QWidget *parent) : QMainWindow(parent), ui(new Ui::MainWindow)
    {
        ui->setupUi(this);
    }
    
    MainWindow::~MainWindow()
    {
        delete ui;
    }
    
    
    void MainWindow::on_first_ok_clicked() {
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
    
    void MainWindow::on_second_ok_clicked() {
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
    
    void MainWindow::on_third_ok_clicked() {
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

        //clock_t begin_sceua;
        //begin_sceua = clock();

        total_count = 0;
        better = 0;
        worse = 0;
        ref = 0;
        con = 0;
        random = 0;

        int tot_dias = ui->num_dias->text().toInt(NULL);
        icall = 0;
        int ngs = 5;
        int iniflg = 1;

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
        float bestf, worstf;


        vector<float> bound, bestx, worstx, xf(npt);

        for(uint i = 0; i < bu.size(); i++)
            bound.push_back(bu.at(i) - bl.at(i));

        vector<vector<float>> x(npt, vector<float>(nopt));

        for(int i = 0; i < npt; i++) {
            for(int j = 0; j < nopt; j++) {
                float r = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
                x.at(i).at(j) = bl.at(j) + r * bound.at(j);
            }
        }

        if (iniflg == 0) {
            for(uint i = 0; i < x.at(0).size(); i++)
                x.at(0).at(i) = x0.at(i);
        }

        bool obj = ui->rmse_radio->isChecked();
        int i;

        omp_set_num_threads(omp_get_max_threads());

        time_t t_st;
        time_t t_end;

        SubBacia subw;

        VarMeteorologicas met;
        VarUsoDoSolo usosolo;
        VarEntrada entrada;

        //Método que seta o número de sub-bacias.
        subw.setNumSub_b(ui->third_file->text());

        //Armazena arquivos de entrada em memória para facilitar acesso aos dados.
        met.criaVetor1(subw, tot_dias);
        met.loadData1(subw, ui->first_file->text(), tot_dias);

        usosolo.criaVetor2(subw, tot_dias);
        usosolo.loadData2(subw, ui->second_file->text(), tot_dias);

        entrada.criaVetor3(subw);
        entrada.loadData3(subw, ui->third_file->text());


        time(&t_st);
        #pragma omp parallel for private(i) schedule(static)
        for(i = 0; i < npt; i++) {
            float r = hydrological_routine(x.at(i), tot_dias, false, subw, met, usosolo, entrada);

            #pragma omp critical
            {
                if(obj)
                    cout << "RMSE " << i << ": " << r << endl;
                else
                    cout << "CNS " << i << ": " << r << endl;

                icall++;
            }
            xf.at(i) = r;
        }
        #pragma omp barrier
        time(&t_end);

        //cout << "Tempo de calibracao para " << tot_dias << " dias com " << omp_get_max_threads() << " threads: " << t_end - t_st << "s.\n";
        //cout << "Tempo de calibracao para " << tot_dias << " dias sequencialmente: " << t_end - t_st << "s.\n";

        cout << "Fim do cálculo de RMSE, iniciando pt2\n";

        vector<pair<float, vector<float>>> link_xxf;

        //liga x com xf
        for(int i = 0; i < npt; i++)
            link_xxf.push_back(make_pair(xf.at(i), x.at(i)));

        sort(link_xxf.begin(), link_xxf.end());
        sort(xf.begin(), xf.end());

        //põe primeira matriz no conjunto de todas as matrizes
        vector<vector<float>> first_matrix(npt, vector<float>(nopt + 2));
        for(int i = 0; i < npt; i++) {
            first_matrix.at(i).at(0) = i +1;

            for(int j = 0; j < nopt; j++)
                first_matrix.at(i).at(j + 1) = link_xxf.at(i).second.at(j);

            first_matrix.at(i).at(nopt + 1) = link_xxf.at(i).first;
        }

        matrices.push_back(first_matrix);
        size_m++;


        //cria matriz x ordenada
        vector<vector<float>> x_ordered;

        for(int i = 0; i < npt; i++)
            x_ordered.push_back(link_xxf.at(i).second);

        bestx = link_xxf.at(0).second;    //bestx recebe melhores parâmetros
        bestf = link_xxf.at(0).first;   //bestf recebe melhor xf

        worstx = link_xxf.at(npt - 1).second; //worstx recebe piores parâmetros
        worstf = link_xxf.at(npt - 1).first; //worstf recebe pior xf

        //calcula média
        float media_rmse = 0;

        #pragma omp parallel for private(i) reduction(+ : media_rmse)
        for(uint i = 0; i < xf.size(); i++)  media_rmse += xf.at(i);

        media_rmse /= xf.size();
        media_vec.push_back(media_rmse);

        cout << "---------------------------------------------------\n";
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

        cout << "AVERAGE (ERRORS): " << media_rmse << endl << endl;

        int nloop = 0;
        size_m = 0;

        cout << "While\n";
        int icount = 1;

        while(icall < maxn) { //loop contando iterações
            srand(time(NULL));
            cout << "Icall: " << icall << endl;
            nloop++;

            int igs;

            #pragma omp parallel for private(igs)
            for(igs = 0; igs < ngs; igs++) {    //0 a 5

                vector<vector<float>> cx(npg, vector<float>(nopt));
                vector<float> cf(npg);

                for(int k1 = 0; k1 < npg; k1++) {   //0 a 15

                    int k2 = k1 * ngs + igs;

                    cf.at(k1) = xf.at(k2);    //cf recebe xf embaralhado
                    cx.at(k1) = x_ordered.at(k2); //cx recebe x embaralhado

                }

                for(int i = 0; i < nspl; i++) { //0 a 15

                    vector<int> lcs;
                    vector<int>::iterator itr;
                    lcs.push_back(0);

                    for(int k3 = 1; k3 < nps; k3++) {
                        int pos;

                        for(int q = 0; q < 1000; q++) {
                            int n_aux = npg - 1; //para gerar random de 1 a 14
                            float r = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);

                            pos = 1 + floor(n_aux + 0.5 - sqrt(pow(n_aux + 0.5, 2.0) - n_aux * (n_aux + 1) * r));

                            itr = std::find(lcs.begin(), lcs.end(), pos);
                            if(itr == lcs.end())
                                break;
                        }

                        lcs.push_back(pos);
                    }

                    sort(lcs.begin(), lcs.end());   //ordena lcs
                    vector<vector<float>> s(nps, vector<float>(nopt));
                    vector<float> sf(nps);

                    for(int j = 0; j < nps; j++) {
                        s.at(j) = cx.at(lcs.at(j));
                        sf.at(j) = cf.at(lcs.at(j));
                    }

                    vector<float> cce = cceua(s, sf, bl, bu, tot_dias, subw, met, usosolo, entrada);

                    vector<float>:: const_iterator first = cce.begin();
                    vector<float>:: const_iterator last = cce.begin() + nopt;
                    vector<float> snew(first, last);

                  /*  #pragma omp critical
                    {
                        if(omp_get_thread_num() == 0) {
                            cout << "IN: ";
                            for(uint ii = 0; ii < cce.size(); ii++) {
                                cout << cce.at(ii) << " ";
                            }
                            cout << endl << "OUT: ";
                            for(int ii = 0; ii < nopt; ii++) {
                                cout << s.at(nps-1).at(ii) << " ";
                            }
                            cout << sf.at(nps-1) << endl;
                        }
                    }   */

                    s.at(nps - 1) = snew;
                    sf.at(nps - 1) = cce.back();

                    //põe o simplex no complex
                    for(int ij = 0; ij < nps; ij++) {
                        cx.at(lcs.at(ij)) = s.at(ij);
                        cf.at(lcs.at(ij)) = sf.at(ij);
                    }

                    //ordenar complex
                    vector<pair<float, vector<float>>> order_c;
                    for(uint kj = 0; kj < cf.size(); kj++)
                        order_c.push_back(make_pair(cf.at(kj), cx.at(kj)));

                    sort(order_c.begin(), order_c.end());

                    for(uint lj = 0; lj < cf.size(); lj++) {
                        cf.at(lj) = order_c.at(lj).first;
                        cx.at(lj) = order_c.at(lj).second;
                    }
                }

                //põe o complex na população
                for(int k1 = 0; k1 < npg; k1++) {   //0 a 15

                    int k2 = k1 * ngs + igs;

                    x_ordered.at(k2) = cx.at(k1);
                    xf.at(k2) = cf.at(k1);
                }

            }
            #pragma omp barrier

            //cout << "Fim: " <<  float(clock() - bt) /  CLOCKS_PER_SEC  << "s" << endl;

            //ordena x de acordo com xf
            vector<pair<float, vector<float>>> new_xf_to_f;
            //first -> xf (rmse)
            //second -> x (vetor de parâmetros)

            for(uint j = 0; j < xf.size(); j++)
                new_xf_to_f.push_back(make_pair(xf.at(j), x_ordered.at(j)));

            sort(new_xf_to_f.begin(), new_xf_to_f.end());

            //melhor x e xf
            bestx = new_xf_to_f.at(0).second;
            bestf = new_xf_to_f.at(0).first;

            //pior x e xf
            worstx = new_xf_to_f.at(npt - 1).second;
            worstf = new_xf_to_f.at(npt - 1).first;

            //calcula média
            media_rmse = 0;
            #pragma omp parallel for private(i) reduction(+ : media_rmse)
            for(uint i = 0; i < xf.size(); i++)  media_rmse += xf.at(i);

            media_rmse /= xf.size();
            media_vec.push_back(media_rmse);

            cout << "\nEvolution loop: " << nloop << " - Trial - " << icall << endl;
            cout << "BESTF: " << bestf << endl;

            cout << "BESTX: ";
            for(uint i = 0; i < bestx.size(); i++)
                cout << bestx.at(i) << " ";

            cout << "\nWORSTF: " << worstf << endl;

            cout << "WORSTX: ";
            for(uint i = 0; i < worstx.size(); i++)
                cout << worstx.at(i) << " ";
            cout << endl;

            cout << "AVERAGE (ERRORS): " << media_rmse << endl << endl;

            cout << "----------------------------------------------\n\n";

            vector<vector<float>> f_matrix(npt, vector<float>(nopt + 2));
            for(int j = 0; j < npt; j++) {

                f_matrix.at(j).at(0) = icount + npt;

                for(int k = 0; k < nopt; k++)
                    f_matrix.at(j).at(k + 1) = new_xf_to_f.at(j).second.at(k);

                f_matrix.at(j).at(nopt + 1) = new_xf_to_f.at(j).first;
                icount++;
            }

            matrices.push_back(f_matrix);
            size_m++;
        }

        cout << "MATRIZES FINAIS:\n---------------------------------------------------\n\n";

        for(int i = 0; i < size_m; i++) {
            for(int j = 0; j < npt; j++) {
                for(int k = 0; k < nopt + 2; k++)
                    cout << matrices.at(i).at(j).at(k) << " ";
                cout << endl;
            }
            cout << "\n---------------------------------------------------\n\n";
        }

        cout << "MEDIA DE ERROS DE CADA ITERACAO:\n";
        for(uint i = 0; i < media_vec.size(); i++)
            cout << i + 1 << ": " << media_vec.at(i) << endl;

        cout << endl;
        cout << "Numero de tentativas: " << total_count << endl;
        cout << "Melhora nos parametros: " << better << endl;
        cout << "Piora nos parametros: " << worse << endl;
        cout << "Reflexao: " << ref << endl;
        cout << "Contracao: " << con << endl;
        cout << "Aleatorio: " << random << endl;

        pair<int, int> bp = best_parameters(npt, nopt);

        vector<float> real_bp;
        for(int i = 0; i < 7; i++)
            real_bp.push_back(matrices.at(bp.first).at(bp.second).at(i + 1));

        cout << "\n\nFIM\n\n";
        for(uint i = 0; i < real_bp.size(); i++)
            cout << real_bp.at(i) << " ";
        cout << endl;

        //cout << "\nTempo total: " << float( clock() - begin_sceua ) /  CLOCKS_PER_SEC << "s\n";

        //GRÁFICO
        hydrological_routine(bestx, tot_dias, true, subw, met, usosolo, entrada);
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
        second->show();

    }

    pair<int, int> MainWindow::best_parameters(int npt, int nopt) {
        float menor = 999999;
        int indice_i, indice_j;

        for(int i = 0; i < size_m; i++) {
            for(int j = 0; j < npt; j++) {
                if(matrices.at(i).at(j).at(nopt + 1) < menor) {
                    menor = matrices.at(i).at(j).at(nopt + 1);
                    indice_i = i;
                    indice_j = j;
                }
            }
        }

        return make_pair(indice_i, indice_j);
    }

    vector<float> MainWindow::cceua(vector<vector<float>> s, vector<float> sf, vector<float> bl, vector<float> bu, int tot_dias,
                                    SubBacia subw, VarMeteorologicas met, VarUsoDoSolo usosolo, VarEntrada entrada) {
        vector<float> sb = s.at(0); //sb melhor s
        vector<float> sw = s.at(6); //sw pior s

        vector<float> snew;

        srand(time(NULL));

        float fw = sf.at(sf.size() - 1);    //fw pior sf
        float alpha = 1.0, beta = 0.5;

        vector<float> ce;

        /*
            ce[0] -> kcr
            ce[1] -> kb
            ce[2] -> kss
            ce[3] -> cs
            ce[4] -> css
            ce[5] -> cb
            ce[6] -> lambda
        */

        for(uint k = 0; k < s.at(0).size(); k++) { //7
            float sum = 0;

            for(uint j =0; j < s.at(0).size(); j++)
                sum += s.at(j).at(k);

            float media = sum / 7; //nps - 1
            ce.push_back(media);

            float svalue = ce.at(k) + alpha * (ce.at(k) - sw.at(k));
            snew.push_back(svalue);
        }

        for(uint j = 0; j < snew.size(); j++) { //verifica os limites
            float s1 = snew.at(j) - bl.at(j);
            float s2 = bu.at(j) - snew.at(j);

            if(s1 < 0 || s2 < 0)
                snew.at(j) = bl.at(j) + ((1 + (rand() % 7 + 1)) * (bu.at(j) - bl.at(j))); //7 nopt
        }

        float fnew;
        bool reflection = true, contraction = false;
        fnew = hydrological_routine(snew, tot_dias, false, subw, met, usosolo, entrada);
        icall++;

        if(fnew > fw) { //reflexão falhou, tentando ponto de contração

            reflection = false;
            contraction = true;

            for(uint j = 0; j < snew.size(); j++)
                snew.at(j) = sw.at(j) + beta * (ce.at(j) - sw.at(j));

            fnew = hydrological_routine(snew, tot_dias, false, subw, met, usosolo, entrada);
            icall++;

            if(fnew > fw) { //reflexão e contração falharam, tentando ponto aleatório

                contraction = false;

                for(uint j = 0; j < snew.size(); j++)
                    snew.at(j) = bl.at(j) + ((rand() % 7 + 1) * (bu.at(j) - bl.at(j)));  //7 -> nopt

                fnew = hydrological_routine(snew, tot_dias, false, subw, met, usosolo, entrada);
                icall++;
            }
        }


        #pragma omp critical
        {
            total_count++;
            if(fnew < fw)
                better++;
            else {	//caso não melhore, mantém o erro antigo e seus parâmetros
               /* fnew = fw;
                for(uint j = 0; j < snew->size(); j++)
                    snew->at(j) = sw.at(j);     */
                worse++;
            }

            if(reflection)
                ref++;
            else if(contraction)
                con++;
            else
                random++;
        }

        snew.push_back(fnew);
        return snew;
    }

    float MainWindow::hydrological_routine(vector<float> x, float tot_dias, bool best_p, SubBacia subw,
                                           VarMeteorologicas met, VarUsoDoSolo usosolo, VarEntrada entrada) {
        int i, j, dia, sub_b = 0;
        float kcr, Kb, Kss, Cs, Css, Cb, Coef_Ia, obj_f = 0, vmedia_obs = 0;
        vector<float> vobs_aux;
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

        VarEvapotranspiracao evapo;
        VarInterceptacao intercep;
        VarSolo solo;
        VarESD ESD;
        VarESS ESS;
        VarEB EB;

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
        int total_dias_observados = tot_dias;
        float dif_valores = 0;

        //cara... faz favor
        int aux = tot_dias;

        for(i = 0; i < aux; i++) {
            float vazao_diaria = 0;
            dia = i + 1;

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


            for(j = 0; j < subw.numSub_b; j++) {
                float vazao_sub = 0, evap_lamina;
                sub_b = j;

                met.setVarMet(dia, sub_b, subw);

                usosolo.setVarUS(dia, sub_b, subw);
                h_sistrad2 = usosolo.h_sistrad;
                entrada.setVarEntrada(sub_b, h_sistrad2);

                evapo.setAll(dia, &evap_lamina, met, entrada, usosolo);
                intercep.setAll(dia, sub_b, evap_lamina, anterior, met, usosolo);

                solo.setAll(dia, sub_b, Coef_Ia, kcr, anterior, entrada, intercep, met);
                evapo.setETr(solo);

                intercep.setAll(dia, sub_b, evap_lamina, anterior, met, usosolo);

                ESD.setAll(dia, sub_b, Cs, anterior, solo, entrada);
                ESS.setAll(dia, sub_b, Kss, Css, solo, anterior, entrada);
                EB.setAll(dia, sub_b, Kb, Cb, solo, anterior, entrada);

                vazao_sub = ESD.vazao_ESD + ESS.vazao_ESS + EB.vazao_EB;
                vazao_diaria += vazao_sub;

                //AJUSTA AS VARIÁVEIS REFERENTES AO DIA ANTERIOR:
                anterior.SetPrec(sub_b, met.precMedia);

                anterior.SetL3(sub_b, intercep.lam3);

                anterior.SetAm(sub_b, solo.am);
                anterior.SetAt(sub_b, solo.at);
                anterior.SetPts(sub_b, solo.pts);
                anterior.SetPe(sub_b, solo.pe);
                anterior.SetDCR(sub_b, solo.dcr);

                anterior.SetETr(sub_b, evapo.ETr);

                anterior.SetVfinal_ESD(sub_b, ESD.vol_final_ESD);

                anterior.SetLam_ESS(sub_b, ESS.lam_ESS);
                anterior.SetVfinal_ESS(sub_b, ESS.vol_final_ESS);

                anterior.SetLam_EB(sub_b, EB.lam_EB);
                anterior.SetVfinal_EB(sub_b, EB.vol_final_EB);

            }

            if(met.dado_observado == -9999){
                total_dias_observados--;
                flag = false;
            }
            else{
                dif_valores += pow((met.dado_observado - vazao_diaria), 2);

                vobs_aux.push_back(met.dado_observado);
                vmedia_obs += met.dado_observado;
                flag = true;
            }

            if(best_p) {
                vazao_observada.push_back(met.dado_observado);
                vazao_calculada.push_back(vazao_diaria);
            }

            vazao_total.push_back(vazao_diaria);
        }

        if(flag) {

            if(ui->rmse_radio->isChecked()) {
                //RMSE
                obj_f = sqrt(dif_valores/total_dias_observados);
            }
            else if(ui->cns_radio->isChecked()) {
                //CNS
                float sum_obs = 0;
                vmedia_obs /= total_dias_observados;

                for(int k = 0; k < total_dias_observados; k++)
                    sum_obs += pow(vobs_aux.at(k) - vmedia_obs, 2);

                obj_f = 1 - (dif_valores/sum_obs);
            }
        }

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

        return obj_f;
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
    
