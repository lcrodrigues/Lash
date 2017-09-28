#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include "varclasses.h"
#include<QTextStream>
#include <stdio.h>
#include <iterator>
#include <vector>

using namespace std;

//===MÉTODOS METEOROLOGICOS===========================
void VarMeteorologicas::criaVetor1(SubBacia subw, int tot_dias)
{
    int col=(7*subw.numSub_b)+4;
    int lin=tot_dias+1;
    dadosMet.resize(lin, std::vector<float>(col, 0));
}
void VarMeteorologicas::loadData1(SubBacia subw, QString filename1, int tot_dias){
    int i, j;
    int col=(7*subw.numSub_b)+4;
    int lin=tot_dias;

    std::string nome_arquivo1;
    nome_arquivo1=filename1.toStdString();

    std::fstream fp;
    fp.open(nome_arquivo1, std::ios::in);

    if(fp.is_open() && fp.good())
    {
        for(i=0; i<lin; i++){
            for(j=0; j<col; j++){
                fp >> this->dadosMet[i][j];
            }
        }
    }
    fp.close();

}
void VarMeteorologicas::setVarMet(int dia, int sub_b, SubBacia subw){

    float U_height = 2;

    this->dia = this->dadosMet[dia-1][0];
    this->mes = this->dadosMet[dia-1][1];
    this->ano = this->dadosMet[dia-1][2];
    this->tempMin = this->dadosMet[dia-1][(4+(sub_b-1))];
    this->tempMax = this->dadosMet[dia-1][4+(subw.numSub_b*1)+(sub_b-1)]; //[(4+(numSubw*1)+(sub_b-1))]
    this->Rh = this->dadosMet[dia-1][(4+(subw.numSub_b*2)+(sub_b-1))];
    this->Rs = this->dadosMet[dia-1][(4+(subw.numSub_b*4)+(sub_b-1))];
    this->U = this->dadosMet[dia-1][(4+(subw.numSub_b*3)+(sub_b-1))];
    this->precMedia = this->dadosMet[dia-1][(4+(subw.numSub_b*5)+(sub_b-1))];
    this->p5 = this->dadosMet[dia-1][(4+(subw.numSub_b*6)+(sub_b-1))];
    this->dado_observado = this->dadosMet[dia - 1][304];

    if(U_height != 2){
        this->U2 = this->U * (4.87 / (log((67.8*U_height)-5.42)));
    }
    else{
        this->U2 = this->U;
    }

}
void VarMeteorologicas::setU2()
{
    float U_height=2;
    if(U_height!=2){
                this->U2=this->U*(4.87 / (log((67.8*U_height)-5.42)));
            }
            else{
                this->U2=this->U;
            }
}
//====================================================


//===MÉTODOS USO DO SOLO==============================
void VarUsoDoSolo::criaVetor2(SubBacia subw, int tot_dias)
{
    int col=(5*subw.numSub_b)+3;
    int lin=tot_dias+1;
    dadosUsoSolo.resize(lin, std::vector<float>(col, 0));
}
void VarUsoDoSolo::loadData2(SubBacia subw, QString filename2, int tot_dias)
{
    int i, j;
    int col=(5*subw.numSub_b)+3;
    int lin=tot_dias;

    std::string nome_arquivo2;
    nome_arquivo2=filename2.toStdString();

    std::fstream fp;
    fp.open(nome_arquivo2, std::ios::in);

    if(fp.is_open() && fp.good())
    {
        for(i=0; i<lin; i++){
            for(j=0; j<col; j++){
                fp >> this->dadosUsoSolo[i][j];
            }
        }
    }
    fp.close();

}
void VarUsoDoSolo::setVarUS(int dia, int sub_b, SubBacia subw)
{
    this->iaf=this->dadosUsoSolo[dia-1][(4+(sub_b-1))];
    this->h_veg=this->dadosUsoSolo[dia-1][4+(subw.numSub_b*1)+(sub_b-1)];
    this->albedo=this->dadosUsoSolo[dia-1][(4+(subw.numSub_b*2)+(sub_b-1))];
    this->stomatal_resistence=this->dadosUsoSolo[dia-1][(4+(subw.numSub_b*3)+(sub_b-1))];
    this->h_sistrad=this->dadosUsoSolo[dia-1][(4+(subw.numSub_b*4)+(sub_b-1))];
}
//====================================================


//===MÉTODOS ENTRADA==================================
void VarEntrada::criaVetor3(SubBacia subw)
{
    int lin=subw.numSub_b;
    int col=8;
    dadosEntrada.resize(lin, std::vector<float>(col, 0));
}
void VarEntrada::loadData3(SubBacia subw, QString filename3)
{
    int i,j;
    int lin=subw.numSub_b;
    int col=8;

    string nome_arquivo3;
    nome_arquivo3=filename3.toStdString();

    fstream fp;
    fp.open(nome_arquivo3, std::ios::in);

    if(fp.is_open() && fp.good()) {
        for(i=0; i<lin; i++){//percorre as linhas da matriz
            for(j=0; j<col; j++){//percorre as colunas da matriz
                fp >> this->dadosEntrada[i][j];
            }
        }
    }
    fp.close();

}
void VarEntrada::setVarEntrada(int sub_b, float h_sistrad)
{

    this->h_balhid=this->dadosEntrada[sub_b][1];
    this->altitude=this->dadosEntrada[sub_b][2];
    this->Upmp=this->dadosEntrada[sub_b][3];
    this->Usaturacao=this->dadosEntrada[sub_b][4];
    this->Tc=this->dadosEntrada[sub_b][5];
    this->latitude=this->dadosEntrada[sub_b][6];
    this->Asubw=this->dadosEntrada[sub_b][7];

    if(h_sistrad > this->h_balhid){
         this->h_balhid = this->dadosEntrada[sub_b][1];
    }
    else{
        this->h_balhid = h_sistrad;
    }
/*

    cout << this->h_balhid << ", " << this->altitude << ", " << this->Upmp << ", " << this->Usaturacao << ", " << this->Tc << ", " << this->latitude << ", " << this->Asubw << endl;
    cout << "===============================================================\n";*/

}
//====================================================


//===MÉTODOS EVAPOTRANSPIRAÇÃO========================
void VarEvapotranspiracao::setTmedia(VarMeteorologicas objmet)
{
    this->Tmedia=(objmet.tempMax + objmet.tempMin)/2;
}
void VarEvapotranspiracao::setEs(VarMeteorologicas objmet)
{
    this->es=((0.6108*(exp(((17.27*objmet.tempMin)/(objmet.tempMin+237.3))))) + (0.6108*(exp(((17.27*objmet.tempMax)/(objmet.tempMax+237.3))))))/2;
}
void VarEvapotranspiracao::setEa(VarMeteorologicas objmet)
{
    this->ea = this->es*(objmet.Rh/100.0);
}
void VarEvapotranspiracao::setSsvpc()
{
    this->ssvpc=(4098*this->es)/pow((this->Tmedia+237.3), 2);
}
void VarEvapotranspiracao::setAp(VarEntrada objentrada)
{
    this->ap=101.3*pow(((293-(0.0065*objentrada.altitude))/293), 5.26);
}
void VarEvapotranspiracao::setPc()
{
   this->pc=0.665*0.001*this->ap;
}
void VarEvapotranspiracao::setRns(VarUsoDoSolo objusodosolo, VarMeteorologicas objmet)
{
    this->rns=(1-objusodosolo.albedo)*objmet.Rs;
}
void VarEvapotranspiracao::setSbc()
{
    this->sbc=(4.903*0.000000001);
}
void VarEvapotranspiracao::setJday(int dia)
{
    this->Jday = dia;
}
void VarEvapotranspiracao::setDr()
{
    this->dr = 1+ (0.033*cos(((2*M_PI)/365.0)*this->Jday));
}
void VarEvapotranspiracao::setSolar_declination()
{
   this->solar_declination=0.409*sin((((2*M_PI)/365.0)*this->Jday)-1.39);
}
void VarEvapotranspiracao::setWs(VarEntrada objentrada)
{
    this->ws=acos(-tan(objentrada.latitude) * tan(this->solar_declination));
}
void VarEvapotranspiracao::setRa(VarEntrada objentrada)
{
    const float Gsc=0.0820; //mj/m2/min
    this->ra = ((24*60)/M_PI)*Gsc*this->dr*((this->ws*sin(objentrada.latitude)*sin(this->solar_declination))+(cos(objentrada.latitude)*cos(this->solar_declination)*sin(this->ws)));
}
void VarEvapotranspiracao::setRso(VarEntrada objentrada)
{
    this->rso=(0.75+(0.00002*objentrada.altitude))*this->ra;
}
void VarEvapotranspiracao::setRnl(VarMeteorologicas objmet)
{
    if(this->rso==0){
        this->rnl=0;
    }
    else if(objmet.Rs<=this->rso){
        this->rnl= this->sbc * ((pow((objmet.tempMax+273.15), 4) + pow((objmet.tempMin+273.15), 4))/2.0) * (0.34-(0.14*sqrt(this->ea))) * ((1.35*(objmet.Rs/this->rso))-0.35);
    }
    else if(objmet.Rs>this->rso){
        this->rnl= this->sbc * (pow(objmet.tempMax+273.15, 4) + pow(objmet.tempMin+273.15, 4))/2 * (0.34-0.14*sqrt(this->ea));
    }
}
void VarEvapotranspiracao::setRn()
{
    this->rn = this->rns-this->rnl;
}
void VarEvapotranspiracao::setTkv()
{
    this->tkv=1.01*(this->Tmedia+273.15);
}
void VarEvapotranspiracao::setU10(VarUsoDoSolo objusodosolo, VarMeteorologicas objmet)
{
    float zo;
    if(objusodosolo.h_veg < 10){
        zo=(objusodosolo.h_veg/10.0);
        this->U10=objmet.U2*((log(10/zo))/(log(2/zo)));
    }
    else{
        zo=(objusodosolo.h_veg/10.0);
        this->U10=objmet.U2*((log(10/zo))/(log(2/zo)));
    }
}
void VarEvapotranspiracao::setAeroResist(VarUsoDoSolo objusodosolo)
{
    float zo;
    if(objusodosolo.h_veg < 10){
        zo=(objusodosolo.h_veg/10.0);
        this->aero_resist=((6.25/this->U10)*(pow(log(10/zo), 2)));
    }
    else{
        zo=(objusodosolo.h_veg/10.0);
        this->aero_resist=94/this->U10;
    }
}
void VarEvapotranspiracao::setETc(VarUsoDoSolo objusodosolo)
{
    const float e=0.622; //ratio molecular weight of water vapour/dry air
    const float R=0.287; //specific gas constant (kj/kg/k)
    int G=0;

    this->ETc= (0.408*this->ssvpc*(this->rn-G)+((86400*this->pc*e)/(this->tkv*R*this->aero_resist))*(this->es-this->ea))/(this->ssvpc+this->pc*(1+(objusodosolo.stomatal_resistence/this->aero_resist)));
}
/*void VarEvapotranspiracao::setKs(VarSolo solo)
{
    if(solo.at<solo.al){
        if((solo.at-solo.Apmp)<0){
            this->ks=0;
        }
        else{
            this->ks=(log(solo.at-solo.Apmp))/(log(solo.al-solo.Apmp));
        }
    }
    else{
        this->ks=1;
    }
    if(this->ks<0){
        this->ks=0;
    }
}*/
void VarEvapotranspiracao::setETr(VarSolo solo)
{
    this->ETr=solo.ks*this->ETc;
}
void VarEvapotranspiracao::setOthers(float *evap_l)
{
    const float e=0.622; //ratio molecular weight of water vapour/dry air
    const float R=0.287; //specific gas constant (kj/kg/k)
    *evap_l = (0.408*this->ssvpc*(this->rn-0)+((86400*this->pc*e)/(this->tkv*R*this->aero_resist))*(this->es-this->ea))/(this->ssvpc+this->pc*(1+(0/this->aero_resist)));
}
void VarEvapotranspiracao::setAll(int dia, float *evap_l, VarMeteorologicas objmet, VarEntrada objentrada, VarUsoDoSolo objusodosolo)
{
    const float Gsc=0.0820; //mj/m2/min
    float zo;
    const float e=0.622; //ratio molecular weight of water vapour/dry air
    const float R=0.287; //specific gas constant (kj/kg/k)
    int G=0;

    this->Tmedia=(objmet.tempMax + objmet.tempMin)/2;
    this->es=((0.6108*(exp(((17.27*objmet.tempMin)/(objmet.tempMin+237.3))))) + (0.6108*(exp(((17.27*objmet.tempMax)/(objmet.tempMax+237.3))))))/2;
    this->ea = this->es*(objmet.Rh/100.0);
    this->ssvpc=(4098*this->es)/pow((this->Tmedia+237.3), 2);
    this->ap=101.3*pow(((293-(0.0065*objentrada.altitude))/293), 5.26);
    this->pc=0.665*0.001*this->ap;
    this->rns=(1-objusodosolo.albedo)*objmet.Rs;
    this->sbc=(4.903*0.000000001);
    this->Jday = dia;
    this->dr = 1+ (0.033*cos(((2*M_PI)/365.0)*this->Jday));
    this->solar_declination=0.409*sin((((2*M_PI)/365.0)*this->Jday)-1.39);
    this->ws=acos(-tan(objentrada.latitude) * tan(this->solar_declination));

    this->ra = ((24*60)/M_PI)*Gsc*this->dr*((this->ws*sin(objentrada.latitude)*sin(this->solar_declination))+(cos(objentrada.latitude)*cos(this->solar_declination)*sin(this->ws)));
    this->rso=(0.75+(0.00002*objentrada.altitude))*this->ra;
    if(this->rso==0){
        this->rnl=0;
    }
    else if(objmet.Rs<=this->rso){
        this->rnl= this->sbc * ((pow((objmet.tempMax+273.15), 4) + pow((objmet.tempMin+273.15), 4))/2.0) * (0.34-(0.14*sqrt(this->ea))) * ((1.35*(objmet.Rs/this->rso))-0.35);
    }
    else if(objmet.Rs>this->rso){
        this->rnl= this->sbc * (pow(objmet.tempMax+273.15, 4) + pow(objmet.tempMin+273.15, 4))/2 * (0.34-0.14*sqrt(this->ea));
    }
    this->rn = this->rns-this->rnl;
    this->tkv=1.01*(this->Tmedia+273.15);

    if(objusodosolo.h_veg < 10){
        zo=(objusodosolo.h_veg/10.0);
        this->U10=objmet.U2*((log(10/zo))/(log(2/zo)));
    }
    else{
        zo=(objusodosolo.h_veg/10.0);
        this->U10=objmet.U2*((log(10/zo))/(log(2/zo)));
    }
    if(objusodosolo.h_veg < 10){
        zo=(objusodosolo.h_veg/10.0);
        this->aero_resist=((6.25/this->U10)*(pow(log(10/zo), 2)));
    }
    else{
        zo=(objusodosolo.h_veg/10.0);
        this->aero_resist=94/this->U10;
    }
    this->ETc= (0.408*this->ssvpc*(this->rn-G)+((86400*this->pc*e)/(this->tkv*R*this->aero_resist))*(this->es-this->ea))/(this->ssvpc+this->pc*(1+(objusodosolo.stomatal_resistence/this->aero_resist)));
    *evap_l = (0.408*this->ssvpc*(this->rn-0)+((86400*this->pc*e)/(this->tkv*R*this->aero_resist))*(this->es-this->ea))/(this->ssvpc+this->pc*(1+(0/this->aero_resist)));
//break
}

//=======================================================

//===MÉTODOS SOLO========================================
void VarSolo::setAm(VarEntrada entrada)
{
    this->am=(entrada.Usaturacao-entrada.Upmp)*entrada.h_balhid;
}
void VarSolo::setAt(int dia, int sub_b, VarDiaAnterior anterior)
{
    if(dia==1){
        this->at=(50/100.0)*this->am;

    }
    if(dia!=1){
        this->at= anterior.At[sub_b] + (anterior.Pts[sub_b][0] - anterior.Pe[sub_b] - anterior.Lam_ESS[sub_b] - anterior.Lam_EB[sub_b] - anterior.ETr[sub_b] + anterior.DCR[sub_b]);
        if(this->at<0){
            this->at=0;

        }
    }

}
        //Variaveis necessárias para calcular At, quando dia!=1.
void VarSolo::setPts(VarInterceptacao intercep, VarMeteorologicas met) //depende de lam1 e lam2!!!
{
    this->pts= met.precMedia-intercep.lam2-intercep.lam1;

    if(this->pts<0){
        this->pts=0;
    }
}
void VarSolo::setS()
{
    this->S=this->am-this->at;
    if(this->S<0){
        this->S=0;
    }
}
void VarSolo::setM(int sub_b, VarDiaAnterior anterior)
{
    float const coef_Ia=0.1145; //BACIA E10;
    float p5=0;
    int i;

        //calcula p5:
    for(i=0; i<5; i++){
        p5=p5+anterior.Pts[sub_b][i];
    }

    if( (0.5) * ((-(1+coef_Ia)*this->S)) + (sqrt((pow((1-coef_Ia), 2)) * (pow(coef_Ia, 2)) + (4*p5*coef_Ia))) ){
        this->M=0;
    }
    else{
        this->M = (0.5 * (((-(1+coef_Ia)*this->S)) + (sqrt((pow((1-coef_Ia), 2)) * (pow(this->S, 2)) + (4*p5*coef_Ia)))));
    }
}
void VarSolo::setIa(float coef_Ia)
{
    //float const coef_Ia=0.1145; //BACIA E10;

    if((this->S==0) && (this->M==0)){
        this->ia=0;
    }
    else{
        //this->ia=coef_Ia*((pow(this->S, 2))/(this->S+this->M));
        this->ia= ((coef_Ia)*(pow(this->S, 2)))/(this->S+this->M);
    }
}
void VarSolo::setPe()
{
    if(this->pts>this->ia){
        this->pe=(((this->pts-this->ia)*(this->pts-this->ia+this->M))/(this->pts-this->ia+this->M+this->S));
    }
    else{
        this->pe=0;
    }
}
void VarSolo::setAcc()
{
    this->acc=(10/100.0)*this->am;
}
void VarSolo::setAc()
{
    this->ac=(1/100.0)*this->am;
}
void VarSolo::setAcr()
{
    this->acr=(10/100.0)*this->am;
}
void VarSolo::setDcr(float Kcr)
{
    //float const Kcr=2.0213;//BACIA E4;

    if(this->acr>this->at){
        this->dcr= Kcr*((this->acr-this->at)/(this->acr));
    }
    else{
        this->dcr=0;
    }
}
void VarSolo::setAl()
{
    this->al=(50/100.0)*this->am;
}
void VarSolo::setApmp()
{
    this->Apmp=(10/100.0)*this->am;
}
void VarSolo::setKs()
{
    if(this->at < this->al){
        if((this->at-this->Apmp)<0){
            this->ks=0;
        }
        else{
            this->ks=(log(this->at-this->Apmp))/(log(this->al-this->Apmp));
        }
    }
    else{
        this->ks=1;
    }
    if(this->ks<0){
        this->ks=0;
    }
}
void VarSolo::setAll(int dia, int sub_b, float coef_Ia, float Kcr, VarDiaAnterior anterior, VarEntrada entrada, VarInterceptacao intercep, VarMeteorologicas met)
{
    float p5=0;
    int i;

    this->am=(entrada.Usaturacao-entrada.Upmp)*entrada.h_balhid;
    if(dia==1){
        this->at=(50/100.0)*this->am;

    }
    if(dia!=1){
        this->at= anterior.At[sub_b] + (anterior.Pts[sub_b][0] - anterior.Pe[sub_b] - anterior.Lam_ESS[sub_b] - anterior.Lam_EB[sub_b] - anterior.ETr[sub_b] + anterior.DCR[sub_b]);
        if(this->at<0){
            this->at=0;

        }
    }
    this->pts= met.precMedia-intercep.lam2-intercep.lam1;

    if(this->pts<0){
        this->pts=0;
    }
    this->S=this->am-this->at;
    if(this->S<0){
        this->S=0;
    }

    //calcula p5:
    for(i=0; i<5; i++){
        p5=p5+anterior.Pts[sub_b][i];
    }
    //*p5tst=p5;

    if( 0.5 * ((-(1+coef_Ia)*this->S) + sqrt( (pow(1-coef_Ia,2)) * (pow(this->S,2)) + (4 * met.p5 * this->S) )) <0){
        this->M=0;
    }
    else{
        this->M=0.5 * ((-(1+coef_Ia)*this->S) + sqrt( (pow(1-coef_Ia,2)) * (pow(this->S,2)) + (4 * met.p5 * this->S) ));
    }

    if((this->S==0) && (this->M==0)){
        this->ia=0;
    }
    else{
        //this->ia=coef_Ia*((pow(this->S, 2))/(this->S+this->M));
        this->ia= ((coef_Ia)*(pow(this->S, 2)))/(this->S+this->M);
    }

    if(this->pts>this->ia){
        this->pe=(((this->pts-this->ia)*(this->pts-this->ia+this->M))/(this->pts-this->ia+this->M+this->S));
    }
    else{
        this->pe=0;
    }
    this->acc=(10/100.0)*this->am;
    this->ac=(1/100.0)*this->am;
    this->acr=(10/100.0)*this->am;

    if(this->acr>this->at){
        this->dcr= Kcr*((this->acr-this->at)/(this->acr));
    }
    else{
        this->dcr=0;
    }

    this->al=(50/100.0)*this->am;
    this->Apmp=(10/100.0)*this->am;

    if(this->at < this->al){
        if((this->at-this->Apmp)<0){
            this->ks=0;
        }
        else{
            this->ks=(log(this->at-this->Apmp))/(log(this->al-this->Apmp));
        }
    }
    else{
        this->ks=1;
    }
    if(this->ks<0){
        this->ks=0;
    }
}

//========================================================


//===MÉTODOS INTERCEPTAÇÃO================================
void VarInterceptacao::setCri(VarUsoDoSolo usosolo)
{
    this->cri=usosolo.iaf*0.2;
}
void VarInterceptacao::setLam1(int dia, int sub_b, VarDiaAnterior anterior)
{
    if(dia==1){
        this->lam1=0;
    }
    else{
        this->lam1=anterior.L3[sub_b]; //dado da matriz = (L3-)cod=207;
    }
}
void VarInterceptacao::setLam2(VarMeteorologicas objmet)
{
    if((objmet.precMedia+this->lam1)>this->cri){
        this->lam2=this->cri;
    }
    else{
        this->lam2=objmet.precMedia+this->lam1;
    }
}
void VarInterceptacao::setEvapLam(float evap_lamina)
{
    this->evap_lam=evap_lamina;
    if(this->evap_lam > this->lam2){
        this->evap_lam = this->lam2;
    }
}
void VarInterceptacao::setLam3()
{
    this->lam3=this->lam2-this->evap_lam;
    //lamina
}
void VarInterceptacao::setAll(int dia, int sub_b, float evap_lamina, VarDiaAnterior anterior, VarMeteorologicas objmet, VarUsoDoSolo usosolo)
{
    this->cri=usosolo.iaf*0.2;
    if(dia==1){
        this->lam1=0;
    }
    else{
        this->lam1=anterior.L3[sub_b]; //dado da matriz = (L3-)cod=207;
    }
    if((objmet.precMedia+this->lam1)>this->cri){
        this->lam2=this->cri;
    }
    else{
        this->lam2=objmet.precMedia+this->lam1;
    }
    this->evap_lam=evap_lamina;
    if(this->evap_lam > this->lam2){
        this->evap_lam = this->lam2;
    }
    this->lam3=this->lam2-this->evap_lam;
}

//=========================================================

//===MÉTODOS ESD===========================================
void VarESD::setVolInicESD(int dia, int sub_b, VarDiaAnterior anterior)
{
    if(dia==1){
        this->vol_inic_ESD=0;
    }
    else{
        this->vol_inic_ESD = anterior.Vfinal_ESD[sub_b];
    }
}
void VarESD::setVolGerESD(VarSolo objsolo, VarEntrada objentrada)
{
    this->vol_ger_ESD = this->vol_inic_ESD+(objsolo.pe*objentrada.Asubw*1000);
}
void VarESD::setVazaoESD(VarEntrada objentrada, float Cs)
{
    //const float Cs=46.5581; //E7

    this->vazao_ESD = (this->vol_ger_ESD)/(objentrada.Tc*Cs*60);
}
void VarESD::setVolFinalESD()
{
    this->vol_final_ESD = this->vol_ger_ESD-(this->vazao_ESD*1440*60);
    if(this->vol_final_ESD<0){
        this->vol_final_ESD=0;
    }
}
void VarESD::setAll(int dia, int sub_b, float Cs, VarDiaAnterior anterior, VarSolo objsolo, VarEntrada objentrada)
{
    if(dia==1){
        this->vol_inic_ESD=0;
    }
    else{
        this->vol_inic_ESD = anterior.Vfinal_ESD[sub_b];
    }

    this->vol_ger_ESD = this->vol_inic_ESD+(objsolo.pe*objentrada.Asubw*1000);
    this->vazao_ESD = (this->vol_ger_ESD)/(objentrada.Tc*Cs*60);

    this->vol_final_ESD = this->vol_ger_ESD-(this->vazao_ESD*1440*60);
    if(this->vol_final_ESD<0){
        this->vol_final_ESD=0;
    }


}

//==========================================================

//===MÉTODOS ESS============================================
void VarESS::setLamESS(VarSolo solo, float Kss)
{
    //float const Kss=5.2609; //BACIA E6;
    float Pr=0.4;

    if(solo.at>=solo.acc){
        this->lam_ESS=Kss*(pow(((solo.at-solo.acc)/(solo.am-solo.acc)),(3+(2/Pr))));
//        anterior.Lam_ESS[sub_b]=this->lam_ESS;

    }
    else{
        this->lam_ESS=0;
//      anterior.Lam_ESS[sub_b]=this->lam_ESS;

    }
}
void VarESS::setVolInicESS(int dia, int sub_b, VarDiaAnterior anterior)
{
    if(dia==1){
        this->vol_inic_ESS=0;
    }
    else{
        this->vol_inic_ESS=anterior.Vfinal_ESS[sub_b];
    }
}
void VarESS::setVolGerESS(VarEntrada objentrada)
{
    this->vol_ger_ESS = this->vol_inic_ESS + (this->lam_ESS*objentrada.Asubw*1000);
}
void VarESS::setVazaoESS(VarEntrada objentrada, float Css)
{
    //const float Css=165.4442;//E8

    this->vazao_ESS = this->vol_ger_ESS / (Css*objentrada.Tc*60);
}
void VarESS::setVolFinalESS()
{
    this->vol_final_ESS = this->vol_ger_ESS-(this->vazao_ESS*1440*60);
    if(this->vol_final_ESS<0){
        this->vol_final_ESS=0;
    }
//    anterior.Vfinal_ESS[sub_b]=this->vol_final_ESS;

}
void VarESS::setAll(int dia, int sub_b, float Kss, float Css, VarSolo solo, VarDiaAnterior anterior, VarEntrada objentrada)
{
    float Pr=0.4;

    if(solo.at>=solo.acc){
        this->lam_ESS=Kss*(pow(((solo.at-solo.acc)/(solo.am-solo.acc)),(3+(2/Pr))));
    }
    else{
        this->lam_ESS=0;
    }
    if(dia==1){
        this->vol_inic_ESS=0;
    }
    else{
        this->vol_inic_ESS=anterior.Vfinal_ESS[sub_b];
    }

    this->vol_ger_ESS = this->vol_inic_ESS + (this->lam_ESS*objentrada.Asubw*1000);
    this->vazao_ESS = this->vol_ger_ESS / (Css * objentrada.Tc * 60);
    this->vol_final_ESS = this->vol_ger_ESS-(this->vazao_ESS*1440*60);

    if(this->vol_final_ESS<0){
        this->vol_final_ESS=0;
    }

    cout << "AT " << solo.at << endl;
    cout << "ACC " << solo.acc << endl;
    cout << "AM " << solo.am << endl;
    cout << "lam ESS " << this->lam_ESS<< endl;
    cout << "vol ini " << this->vol_inic_ESS << endl;
    cout << "vol ger " << this->vol_ger_ESS << endl;
    cout << "vazao " << this->vazao_ESS << endl;
    cout << "==============================================\n";
}

//==========================================================

//===MÉTODOS EB=============================================
void VarEB::setLamEB(VarSolo solo, float Kb)
{
    //float const Kb=0.6338; //BACIA E5;

    if(solo.at>=solo.ac){
        this->lam_EB=Kb*((solo.at-solo.ac)/(solo.am-solo.ac));
//        anterior.Lam_EB[sub_b]=this->lam_EB;
    }
    else{
        this->lam_EB=0;
//        anterior.Lam_EB[sub_b]=this->lam_EB;
    }

}
void VarEB::setVolInicEB(int dia, int sub_b, VarDiaAnterior anterior)
{
    if(dia==1){
        this->vol_inic_EB=0;
    }
    else{
        this->vol_inic_EB = anterior.Vfinal_EB[sub_b];
    }
}
void VarEB::setVolGerEB(VarEntrada objentrada)
{
    this->vol_ger_EB= this->vol_inic_EB+(this->lam_EB*objentrada.Asubw*1000);
}
void VarEB::setVazaoEB(float Cb)
{
    //const float Cb=57.5433;//E9
    this->vazao_EB = this->vol_ger_EB/(Cb*24*60*60);
}
void VarEB::setVolFinalEB()
{
    this->vol_final_EB = this->vol_ger_EB-(this->vazao_EB*1440*60);
    if(this->vol_final_EB<0){
        this->vol_final_EB=0;
    }
}
void VarEB::setAll(int dia, int sub_b, float Kb, float Cb, VarSolo solo, VarDiaAnterior anterior, VarEntrada objentrada)
{
    if(solo.at>=solo.ac){
        this->lam_EB=Kb*((solo.at-solo.ac)/(solo.am-solo.ac));
    }
    else{
        this->lam_EB=0;
    }
    if(dia==1){
        this->vol_inic_EB=0;
    }
    else{
        this->vol_inic_EB = anterior.Vfinal_EB[sub_b];
    }

    this->vol_ger_EB= this->vol_inic_EB+(this->lam_EB*objentrada.Asubw*1000);
    this->vazao_EB = this->vol_ger_EB/(Cb*24*60*60);
    this->vol_final_EB = this->vol_ger_EB-(this->vazao_EB*1440*60);
    if(this->vol_final_EB<0){
        this->vol_final_EB=0;
    }

}

//==========================================================

//===MÉTODOS DIAS_ANTERIORES================================
void VarDiaAnterior::iniciaVetores(SubBacia subw)
{
    Prec.resize(subw.numSub_b);
    Pts.resize(subw.numSub_b, std::vector<float>(5, 0));
    //Pts.resize(subw.numSub_b);
    Am.resize(subw.numSub_b);
    At.resize(subw.numSub_b);
    Pe.resize(subw.numSub_b);
    Lam_ESS.resize(subw.numSub_b);
    Lam_EB.resize(subw.numSub_b);
    ETr.resize(subw.numSub_b);
    DCR.resize(subw.numSub_b);
    L3.resize(subw.numSub_b);
    Vfinal_ESD.resize(subw.numSub_b);
    Vfinal_ESS.resize(subw.numSub_b);
    Vfinal_EB.resize(subw.numSub_b);
}

void VarDiaAnterior::SetPrec(int pos, float valor){
    this->Prec[pos]=valor;
}
void VarDiaAnterior::SetPts(int pos, float valor){
    int i;

    for(i=4; i>=0; i--){
        this->Pts[pos][i]=this->Pts[pos][i-1];
    }
    this->Pts[pos][0]=valor;
}
void VarDiaAnterior::SetAm(int pos, float valor){
    this->Am[pos]=valor;
}
void VarDiaAnterior::SetAt(int pos, float valor){
    this->At[pos]=valor;
}
void VarDiaAnterior::SetPe(int pos, float valor){
    this->Pe[pos]=valor;
}
void VarDiaAnterior::SetLam_ESS(int pos, float valor){
    this->Lam_ESS[pos]=valor;
}
void VarDiaAnterior::SetLam_EB(int pos, float valor){
    this->Lam_EB[pos]=valor;
}
void VarDiaAnterior::SetETr(int pos, float valor){
    this->ETr[pos]=valor;
}
void VarDiaAnterior::SetDCR(int pos, float valor){
    this->DCR[pos]=valor;
}
void VarDiaAnterior::SetL3(int pos, float valor){
    this->L3[pos]=valor;
}
void VarDiaAnterior::SetVfinal_ESD(int pos, float valor){
    this->Vfinal_ESD[pos]=valor;
}
void VarDiaAnterior::SetVfinal_ESS(int pos, float valor){
    this->Vfinal_ESS[pos]=valor;
}
void VarDiaAnterior::SetVfinal_EB(int pos, float valor){
    this->Vfinal_EB[pos]=valor;
}


void SubBacia::setNumSub_b(QString filename)
{
    std::fstream fp;
    std::string nome_arquivo;
    char c;
    int count=0;

    nome_arquivo=filename.toStdString();

    fp.open(nome_arquivo, std::ios::in);

    if (fp.is_open())
    {
        while(fp.get(c)){
            if(c == '\n'){
                count++;
            }
        }
    }
    else std::cout << "Unable to open file";

    this->numSub_b=count;
    fp.close();
}

/*
void SCE_UA::sceua_routine(vector<float> x0, vector<float> bu, vector<float> bl, int maxn,
                                              int peps, int ngs, float iseed, int iniflg) {









}   */

void CCE_UA::cceua_routine(vector< vector<float> > s, vector<float> sf, vector<float> bu, vector<float> bl, int icall, int maxn) {

    vector< vector<float> > result;
    this->result = result;
}


