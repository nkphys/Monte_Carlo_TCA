#include <math.h>
#include "tensor_type.h"
#include "ParametersEngine.h"
#include "Coordinates.h"
#include "random"
#include <stdlib.h>
#define PI acos(-1.0)

#ifndef MFParams_class
#define MFParams_class

class MFParams{ 
public:

    // Define Fields
    Matrix<double> etheta, ephi;
    Matrix<double> etheta_avg, ephi_avg;

    // Constructor
    MFParams(Parameters& Parameters__, Coordinates&  Coordinates__, mt19937_64& Generator__)
        :Parameters_(Parameters__),Coordinates_(Coordinates__), Generator_(Generator__)
    {
        //setupm_arr();
        initialize();
    }


    double random();
    void FieldThrow(int site);
    void initialize();
    void Adjust_MCWindow();
    void Calculate_Fields_Avg();
    void Read_classical_DOFs(string filename);


    Parameters &Parameters_;
    Coordinates &Coordinates_;
    mt19937_64 &Generator_;
    int lx_,ly_,ns_;
    uniform_real_distribution<double> dis_;
    //mt19937_64 mt_rand(Parameters_.RandomSeed);


};


void MFParams::Adjust_MCWindow(){
    double ratio;
    ratio=Parameters_.AccCount[0]/(Parameters_.AccCount[0]+Parameters_.AccCount[1]);
    //cout<<"ratio= "<< ratio << "temp= "<<Parameters_.temp << endl;
    Parameters_.AccCount[0]=0; Parameters_.AccCount[1]=0;
    Parameters_.WindowSize *= abs(1.0 + 1.0*(ratio-0.5));
    //Parameters_.WindowSize =1.0;
    cout << "Ratio: " << ratio << "  window size:  "<<Parameters_.WindowSize<< endl;
    return;
} // ----------


void MFParams::FieldThrow(int site){
    int a,b;

    double Pi=Parameters_.pi;
    double MC_Window = Parameters_.WindowSize;

    a=Coordinates_.indx(site);
    b=Coordinates_.indy(site);

    ephi(a,b) += 2*Pi*(random()-0.5)*MC_Window;
    if( ephi(a,b) < 0.0) {ephi(a,b) += 2.0*Pi; }
    if( ephi(a,b) >=2.0*Pi) {ephi(a,b) -= 2.0*Pi;}


    etheta(a,b) +=  Pi*(random()-0.5)*MC_Window;
    if ( etheta(a,b) < 0.0 ) {
        etheta(a,b) = - etheta(a,b);
        ephi(a,b) = fmod( ephi(a,b)+Pi, 2.0*Pi );
    }
    if ( etheta(a,b) > Pi ) {
        etheta(a,b) -= 2.0*Pi ;
        ephi(a,b) = fmod( ephi(a,b) + Pi, 2.0*Pi );
    }

} // ----------


double MFParams::random(){

    /*
    double random_double;
    random_double=(rand()%RAND_MAX);
    random_double=random_double/RAND_MAX;

    return random_double;
    */

    return dis_(Generator_);

}

void MFParams::initialize(){


    lx_=Coordinates_.lx_;
    ly_=Coordinates_.ly_;

   // srand(Parameters_.RandomSeed);


    etheta_avg.resize(lx_,ly_);
    ephi_avg.resize(lx_,ly_);

    etheta.resize(lx_,ly_);
    ephi.resize(lx_,ly_);

    for(int j=0;j<ly_;j++){
        for(int i=0;i<lx_;i++){
            //ephi(i,j)=(0.5+0.5*pow(-1.0f,i))*Parameters_.pi + grnd()*0.2;
            //etheta(i,j)=0.5*Parameters_.pi + grnd()*0.2;

            //q=(pi,pi)
           // ephi(i,j)=0.0; //(0.5+0.5*pow(-1.0f,i))*Parameters_.pi + grnd()*0.2;
           // etheta(i,j)=0.5*(pow(-1.0,j+i)  + 1.0 )*PI ;//+ grnd()*0.2;

            //q=(0,pi)
            //ephi(i,j)=0.0; //(0.5+0.5*pow(-1.0f,i))*Parameters_.pi + grnd()*0.2;
            //etheta(i,j)=0.5*(pow(-1.0,j)  + 1.0 )*PI; //+ grnd()*0.2;

            //q=(0,0)
           // ephi(i,j)=0.0; //(0.5+0.5*pow(-1.0f,i))*Parameters_.pi + grnd()*0.2;
           // etheta(i,j)=0.0; //+ grnd()*0.2;

            //RANDOM
            ephi(i,j)=2.0*random()*PI;
            etheta(i,j)=random()*PI;
        }
    }


} // ----------

void MFParams::Calculate_Fields_Avg(){

    for(int j=0;j<ly_;j++){
        for(int i=0;i<lx_;i++){

            ephi_avg(i,j)= ephi_avg(i,j) + ephi(i,j);
            etheta_avg(i,j)=etheta_avg(i,j) + etheta(i,j) ;
        }
    }


} // ----------


void MFParams::Read_classical_DOFs(string filename){

    string tmp_str;
    double tmp_double;
    ifstream fl_in(filename.c_str());
    fl_in >> tmp_str;

    for(int i=0;i<lx_;i++){
        for(int j=0;j<ly_;j++){
            fl_in>>tmp_double>>tmp_double>>etheta(i,j)>>ephi(i,j);
        }
    }


} // ----------


#endif
