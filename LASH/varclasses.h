#ifndef VARCLASSES_H
#define VARCLASSES_H

#include <QTextStream>
#include <vector>

using namespace std;
//==VARIAVEIS GLOBAIS======

//=========================


//=======CLASSES============
class SubBacia
{
public:
    int numSub_b;
public:
    void setNumSub_b(QString filename);
};
class VarDiaAnterior
{
public:
    std::vector <float> Prec;
    //    float Prec[43];
    std::vector< std::vector<float> > Pts;
    //    float Pts[43];
    std::vector <float> Am;
    //    float Am[43];
    std::vector <float> At;
    //    float At[43];
    std::vector <float> Pe;
    //    float Pe[43];
    std::vector <float> Lam_ESS;
    //    float Lam_ESS[43];
    std::vector <float> Lam_EB;
    //    float Lam_EB[43];
    std::vector <float> ETr;
    //    float ETr[43];
    std::vector <float> DCR;
    //    float DCR[43];
    std::vector <float> L3;
    //    float L3[43];
    std::vector <float> Vfinal_ESD;
    //    float Vfinal_ESD[43];
    std::vector <float> Vfinal_ESS;
    //    float Vfinal_ESS[43];
    std::vector <float> Vfinal_EB;
    //    float Vfinal_EB[43];
public:
    void iniciaVetores(SubBacia subw);
    void SetPrec(int pos, float valor);
    void SetPts(int pos, float valor);
    void SetAm(int pos, float valor);
    void SetAt(int pos, float valor);
    void SetPe(int pos, float valor);
    void SetLam_ESS(int pos, float valor);
    void SetLam_EB(int pos, float valor);
    void SetETr(int pos, float valor);
    void SetDCR(int pos, float valor);
    void SetL3(int pos, float valor);
    void SetVfinal_ESD(int pos, float valor);
    void SetVfinal_ESS(int pos, float valor);
    void SetVfinal_EB(int pos, float valor);
};

class VarEntrada
{
//Váriveis referentes aos dados do arquivo "mapaslash".
public:
    std::vector< std::vector<float> > dadosEntrada;
    float h_balhid;
    float altitude;
    float latitude;
    float Upmp;
    float Usaturacao;
    float Asubw;
    float Tc;
    //float U_height;
public:
    void criaVetor3(SubBacia subw);
    void loadData3(SubBacia subw, QString filename3);
    void setVarEntrada(int sub_b, float h_sistrad);
    //declarar metodos.
};
class VarMeteorologicas
{
//Váriaveis referentes aos dados do arquivo "Dados".
public:
    std::vector< std::vector<float> > dadosMet;
    float dia;
    float mes;
    float ano;
    float precMedia;
    float tempMin;
    float tempMax;
    float Rh;
    float Rs;
    float U;
    float p5;
    float U2;
    float dado_observado;
public:
    void criaVetor1(SubBacia subw, int tot_dias);
    void loadData1(SubBacia subw, QString filename1, int tot_dias);
    void setVarMet(int dia, int sub_b, SubBacia subw);
    void setU2();
    //declarar metodos.
};
class VarUsoDoSolo
{
//Variáveis referentes aos dados do arquivo "Uso_Solo_Fragata".
public:
    std::vector< std::vector<float> > dadosUsoSolo;
    float iaf;
    float h_veg;
    float albedo;
    float stomatal_resistence;
    float h_sistrad;
public:
    void criaVetor2(SubBacia subw, int tot_dias);
    void loadData2(SubBacia subw, QString filename2, int tot_dias);
    void setVarUS(int dia, int sub_b, SubBacia subw);
    //definir metodos.
};
class VarInterceptacao
{
public:
    float cri;
    float lam1;
    float lam2;
    float lam3;
    float evap_lam;
public:
    void setCri(VarUsoDoSolo usosolo);
    void setLam1(int dia, int sub_b, VarDiaAnterior anterior);
    void setLam2(VarMeteorologicas objmet);
    void setLam3();
    void setEvapLam(float evap_lamina);
    void setAll(int dia, int sub_b, float evap_lamina, VarDiaAnterior anterior, VarMeteorologicas objmet, VarUsoDoSolo usosolo);
    //definir metodos.
};
class VarSolo
{
public:
    float at;
    float am;
    float al;
    float Apmp;
    float acc;
    float ac;
    float acr;
    float pts;
    float dcr;
    float S;
    float M;
    float ia;
    float pe;
    float ks;
public:
    void setAm(VarEntrada entrada);
    void setAt(int dia, int sub_b, VarDiaAnterior anterior);
    void setPts(VarInterceptacao interceptacao, VarMeteorologicas met);
    void setS();
    void setM(int sub_b, VarDiaAnterior anterior);
    void setIa(float coef_Ia);
    void setPe();
    void setAcc();
    void setAc();
    void setAcr();
    void setDcr(float Kcr);
    void setAl();
    void setApmp();
    void setKs();
    void setAll(int dia, int sub_b, float coef_Ia, float Kcr, VarDiaAnterior anterior, VarEntrada entrada, VarInterceptacao intercep, VarMeteorologicas met);
    //definir metodos.
};
class VarEvapotranspiracao
{
public:
    float Tmedia;
    float Jday;
    float dr;
    float solar_declination;
    float ws;
    float ra;
    float rso;
    float es;
    float ea;
    float rnl;
    float rns;
    float tkv;
    float ssvpc;
    float U10;
    float aero_resist;
    float ETc;
    float ETr;
    //float ks;
    float ap;
    float pc;
    float sbc;
    float rn;
public:
    void setTmedia(VarMeteorologicas objmet);
    void setEs(VarMeteorologicas objmet);
    void setEa(VarMeteorologicas objmet);
    void setSsvpc();
    void setAp(VarEntrada objentrada);
    void setPc();
    void setRns(VarUsoDoSolo objusodosolo, VarMeteorologicas objmet);
    void setSbc();
    void setJday(int dia);
    void setDr();
    void setSolar_declination();
    void setWs(VarEntrada objentrada);
    void setRa(VarEntrada objentrada);
    void setRso(VarEntrada objentrada);
    void setRnl(VarMeteorologicas objmet);
    void setRn();
    void setTkv();
    void setU10(VarUsoDoSolo objusodosolo, VarMeteorologicas objmet);
    void setAeroResist(VarUsoDoSolo objusodosolo);
    void setETc(VarUsoDoSolo objusodosolo);
    //void setKs(VarSolo solo);
    void setETr(VarSolo solo);
    void setOthers(float *evap_l);
    void setAll(int dia, float *evap_l, VarMeteorologicas objmet, VarEntrada objentrada, VarUsoDoSolo objusodosolo);
    //definir metodos.
};
class VarESD
{
public:
    float vol_inic_ESD;
    float vol_ger_ESD;
    float vazao_ESD;
    float vol_final_ESD;
public:
    void setVolInicESD(int dia, int sub_b, VarDiaAnterior anterior);
    void setVolGerESD(VarSolo objsolo, VarEntrada objentrada);
    void setVolFinalESD();
    void setVazaoESD(VarEntrada objentrada, float Cs);
    void setAll(int dia, int sub_b, float Cs, VarDiaAnterior anterior, VarSolo objsolo, VarEntrada objentrada);
    //definir metodos.
};
class VarESS
{
public:
    float vol_inic_ESS;
    float vol_ger_ESS;
    float vazao_ESS;
    float vol_final_ESS;
    float lam_ESS;
public:
    void setLamESS(VarSolo solo, float Kss);
    void setVolInicESS(int dia, int sub_b, VarDiaAnterior anterior);
    void setVolGerESS(VarEntrada objentrada);
    void setVolFinalESS();
    void setVazaoESS(VarEntrada objentrada, float Css);
    void setAll(int dia, int sub_b, float Kss, float Css, VarSolo solo, VarDiaAnterior anterior, VarEntrada objentrada);
    //definir metodos.
};
class VarEB
{
public:
    float vol_inic_EB;
    float vol_ger_EB;
    float vazao_EB;
    float vol_final_EB;
    float lam_EB;
public:
    void setLamEB(VarSolo solo, float Kb);
    void setVolInicEB(int dia, int sub_b, VarDiaAnterior anterior);
    void setVolGerEB(VarEntrada objentrada);
    void setVolFinalEB();
    void setVazaoEB(float Cb);
    void setAll(int dia, int sub_b, float Kb, float Cb, VarSolo solo, VarDiaAnterior anterior, VarEntrada objentrada);
    //definir metodos.
};

class SCE_UA {

public:
    vector< vector<float> > result;

public:
    void sceua_routine(vector<float> x0, vector<float> bu, vector<float> bl, int maxn,
                                          int peps, int ngs, float iseed, int iniflg);
};

class CCE_UA {

public:
    vector< vector<float> > result;

public:
    void cceua_routine(vector< vector<float> > s, vector<float> sf, vector<float> bu, vector<float> bl, int icall, int maxn);
};

//==========================

#endif // VARCLASSES_H
