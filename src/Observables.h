#include "ParametersEngine.h"
#include "Coordinates.h"
#include "MFParams.h"
#include "Hamiltonian.h"
#include "tensor_type.h"

#ifndef OBSERVABLES_H
#define OBSERVABLES_H

#define PI acos(-1.0)

class Observables{
public:

    Observables(Parameters& Parameters__, Coordinates& Coordinates__,
                MFParams& MFParams__, Hamiltonian& Hamiltonian__)
        : Parameters_(Parameters__), Coordinates_(Coordinates__), MFParams_(MFParams__),
          Hamiltonian_(Hamiltonian__),
          lx_(Parameters_.lx), ly_(Parameters_.ly), ns_(Parameters_.ns)
    {
        Initialize();
    }

    void Initialize();

    void OccDensity();
    void Calculate_Akw();
    void Calculate_Nw();
    void Get_Non_Interacting_dispersion();
    double Lorentzian(double x, double brd);
    void TotalOccDensity();
    void DensityOfStates();
    void SiSjFULL();

    void SiSjQ_Average();
    void SiSj_Average();
    void Total_Energy_Average(double Curr_QuantE, double CurrE);

    void OccDensity(int tlabel);
    void DOSprint(int tlabel);
    complex<double> SiSjQ(int i,int j);
    double SiSj(int i,int j);
    double Omega(int i);

    complex<double> SiSjQ_Mean(int i, int j);
    complex<double> SiSjQ_square_Mean(int i, int j);

    double SiSj_Mean(int i, int j);
    double SiSj_square_Mean(int i, int j);

    double BandWidth;
    double nia_t,nib_t,nic_t,n_t;
    Matrix<complex<double>> SiSjQ_, SiSjQ_Mean_, SiSjQ_square_Mean_;
    Matrix<double> SiSj_Mean_, SiSj_square_Mean_;
    double Nematic_order_mean_, Nematic_order_square_mean_;
    Parameters& Parameters_;
    Coordinates& Coordinates_;
    MFParams& MFParams_;
    Hamiltonian& Hamiltonian_;
    int lx_,ly_,ns_;
    double dosincr_,tpi_;
    vector<double> nia_,nib_,nic_;
    Matrix<double> SiSj_,dos;
    vector<double> sx_,sy_,sz_;
    double AVG_Total_Energy, AVG_Total_Energy_sqr;
};
/*
 * ***********
 *  Functions in Class Observables ------
 *  ***********
*/


void Observables::Calculate_Akw(){


    //---------Read from input file-----------------------//
    string fileout="Akw.txt";
    double omega_min, omega_max, d_omega;
    double eta = 0.001;
    omega_min=-1.6;omega_max=2.6;d_omega=0.03;
    //---------------------------------------------------//


    int omega_index_max = int( (omega_max - omega_min)/(d_omega) );

    ofstream file_Akw_out(fileout.c_str());

    int c1,c2;

    Mat_3_Complex_doub A_up_00, A_up_11, A_up_22;
    Mat_3_Complex_doub A_dn_00, A_dn_11, A_dn_22;
    A_up_00.resize(Parameters_.ns);
    A_up_11.resize(Parameters_.ns);
    A_up_22.resize(Parameters_.ns);
    A_dn_00.resize(Parameters_.ns);
    A_dn_11.resize(Parameters_.ns);
    A_dn_22.resize(Parameters_.ns);
    for (int i=0;i<Parameters_.ns;i++){
        A_up_00[i].resize(Parameters_.ns);
        A_up_11[i].resize(Parameters_.ns);
        A_up_22[i].resize(Parameters_.ns);
        A_dn_00[i].resize(Parameters_.ns);
        A_dn_11[i].resize(Parameters_.ns);
        A_dn_22[i].resize(Parameters_.ns);
        for(int j=0;j<Parameters_.ns;j++){
            A_up_00[i][j].resize(omega_index_max);
            A_up_11[i][j].resize(omega_index_max);
            A_up_22[i][j].resize(omega_index_max);
            A_dn_00[i][j].resize(omega_index_max);
            A_dn_11[i][j].resize(omega_index_max);
            A_dn_22[i][j].resize(omega_index_max);
        }
    }


    complex<double> Nup_check(0,0);
    complex<double> Ndn_check(0,0);

    for (int j=0;j<Parameters_.ns;j++){
        for (int l=0;l<Parameters_.ns;l++){
            for(int omega_ind=0;omega_ind<omega_index_max;omega_ind++){
                A_up_00[j][l][omega_ind]=zero_complex;
                A_up_11[j][l][omega_ind]=zero_complex;
                A_up_22[j][l][omega_ind]=zero_complex;
                A_dn_00[j][l][omega_ind]=zero_complex;
                A_dn_11[j][l][omega_ind]=zero_complex;
                A_dn_22[j][l][omega_ind]=zero_complex;

                for(int n=0;n<Hamiltonian_.Ham_.n_row();n++){

                    //c= l + or1*ns_ + ns_*orbs_*spin;

                    //Hamiltonian_.Ham_(c2,n) is nth eigenvector and c2th component [checked];
                    c1 = l + ns_*Parameters_.orbs; c2 = j+ ns_*Parameters_.orbs;
                    A_dn_00[j][l][omega_ind] +=  conj(Hamiltonian_.Ham_(c1,n))*Hamiltonian_.Ham_(c2,n)*
                            Lorentzian( omega_min + (omega_ind*d_omega) - Hamiltonian_.eigs_[n], eta);

                    c1 = l+ ns_+ ns_*Parameters_.orbs; c2 = j+ ns_+ ns_*Parameters_.orbs;
                    A_dn_11[j][l][omega_ind] +=  conj(Hamiltonian_.Ham_(c1,n))*Hamiltonian_.Ham_(c2,n)*
                            Lorentzian( omega_min + (omega_ind*d_omega) - Hamiltonian_.eigs_[n], eta);

                    c1 = l + 2*ns_+ ns_*Parameters_.orbs; c2 = j + 2*ns_+ ns_*Parameters_.orbs;
                    A_dn_22[j][l][omega_ind] +=  conj(Hamiltonian_.Ham_(c1,n))*Hamiltonian_.Ham_(c2,n)*
                            Lorentzian( omega_min + (omega_ind*d_omega) - Hamiltonian_.eigs_[n], eta);



                    c1 = l; c2 = j;
                    A_up_00[j][l][omega_ind] +=  conj(Hamiltonian_.Ham_(c1,n))*Hamiltonian_.Ham_(c2,n)*
                            Lorentzian( omega_min + (omega_ind*d_omega) - Hamiltonian_.eigs_[n], eta);

                    c1 = l+ ns_; c2 = j+ ns_;
                    A_up_11[j][l][omega_ind] +=  conj(Hamiltonian_.Ham_(c1,n))*Hamiltonian_.Ham_(c2,n)*
                            Lorentzian( omega_min + (omega_ind*d_omega) - Hamiltonian_.eigs_[n], eta);

                    c1 = l + 2*ns_; c2 = j + 2*ns_;
                    A_up_22[j][l][omega_ind] +=  conj(Hamiltonian_.Ham_(c1,n))*Hamiltonian_.Ham_(c2,n)*
                            Lorentzian( omega_min + (omega_ind*d_omega) - Hamiltonian_.eigs_[n], eta);






                }

                if(j==l){
                    Nup_check += (A_up_00[j][l][omega_ind] + A_up_11[j][l][omega_ind] + A_up_22[j][l][omega_ind])*d_omega;
                    Ndn_check += (A_dn_00[j][l][omega_ind] + A_dn_11[j][l][omega_ind] + A_dn_22[j][l][omega_ind])*d_omega;
                }
            }
        }
    }

    cout << "Nup_check = "<<Nup_check<<endl;
    cout << "Ndn_check = "<<Ndn_check<<endl;

    complex<double> temp_up_00,temp_up_11,temp_up_22;
    complex<double> temp_dn_00,temp_dn_11,temp_dn_22;
    double kx,ky;
    int kx_i,ky_i;

    Mat_1_intpair k_path;
    k_path.clear();
    pair_int temp_pair;

    //--------\Gamma to X-----------------
    ky_i=0;
    for(kx_i=0;kx_i<=(Parameters_.lx/2);kx_i++){
        temp_pair.first = kx_i;
        temp_pair.second = ky_i;
        k_path.push_back(temp_pair);
    }
    //----------------------------------

    //--------X to M-----------------
    kx_i=(Parameters_.lx/2);
    for(ky_i=1;ky_i<=(Parameters_.lx/2);ky_i++){
        temp_pair.first = kx_i;
        temp_pair.second = ky_i;
        k_path.push_back(temp_pair);
    }
    //----------------------------------

    //--------M to \Gamma[with one extra point,
    //                  because in gnuplor use "set pm3d corners2color c1"
    //                  ]-----------------
    kx_i=(Parameters_.lx/2) - 1;
    ky_i=(Parameters_.lx/2) - 1;
    for(kx_i=(Parameters_.lx/2) - 1;kx_i>=-1;kx_i--){
        temp_pair.first = kx_i;
        temp_pair.second = kx_i;
        k_path.push_back(temp_pair);
    }
    //----------------------------------


    double k22_offset=0;
    for(int k_point=0;k_point<k_path.size();k_point++){

        kx_i=k_path[k_point].first;
        ky_i=k_path[k_point].second;
        kx=(2.0*PI*kx_i)/(1.0*Parameters_.lx);
        ky=(2.0*PI*ky_i)/(1.0*Parameters_.ly);

        for(int omega_ind=0;omega_ind<omega_index_max;omega_ind++){
            temp_up_00=zero_complex;temp_up_11=zero_complex;temp_up_22=zero_complex;
            temp_dn_00=zero_complex;temp_dn_11=zero_complex;temp_dn_22=zero_complex;

            for(int j=0;j<ns_;j++){
                for(int l=0;l<ns_;l++){
                    temp_up_00 += one_complex*
                            exp(iota_complex*(kx*(Coordinates_.indx(j) - Coordinates_.indx(l)) +
                                              ky*(Coordinates_.indy(j) - Coordinates_.indy(l))))*
                            A_up_00[j][l][omega_ind];
                    temp_up_11 += one_complex*
                            exp(iota_complex*(kx*(Coordinates_.indx(j) - Coordinates_.indx(l)) +
                                              ky*(Coordinates_.indy(j) - Coordinates_.indy(l))))*
                            A_up_11[j][l][omega_ind];
                    temp_up_22 += one_complex*
                            exp(iota_complex*((kx+k22_offset)*(Coordinates_.indx(j) - Coordinates_.indx(l)) +
                                              (ky+k22_offset)*(Coordinates_.indy(j) - Coordinates_.indy(l))))*
                            A_up_22[j][l][omega_ind];
                    temp_dn_00 += one_complex*
                            exp(iota_complex*(kx*(Coordinates_.indx(j) - Coordinates_.indx(l)) +
                                              ky*(Coordinates_.indy(j) - Coordinates_.indy(l))))*
                            A_dn_00[j][l][omega_ind];
                    temp_dn_11 += one_complex*
                            exp(iota_complex*(kx*(Coordinates_.indx(j) - Coordinates_.indx(l)) +
                                              ky*(Coordinates_.indy(j) - Coordinates_.indy(l))))*
                            A_dn_11[j][l][omega_ind];
                    temp_dn_22 += one_complex*
                            exp(iota_complex*((kx+k22_offset)*(Coordinates_.indx(j) - Coordinates_.indx(l)) +
                                              (ky+k22_offset)*(Coordinates_.indy(j) - Coordinates_.indy(l))))*
                            A_dn_22[j][l][omega_ind];

                }
            }
            //Use 1:6:7----for gnuplot
            file_Akw_out<< k_point<<"   "<<kx_i<<"   "<<ky_i<<"   "<<(ky_i*Parameters_.lx) + kx_i<<"    "<<
                           omega_min + (d_omega*omega_ind)<<"   "<<omega_ind<<"    "<<temp_up_00.real()<<"    "<<temp_up_11.real()
                        <<"    "<<
                          temp_up_22.real()<<"    "<<temp_dn_00.real()<<"    "<<temp_dn_11.real()<<"    "<<temp_dn_22.real()<<"    "<<temp_up_00.imag()<<"    "<<temp_up_11.imag()
                       <<"    "<<
                         temp_up_22.imag()<<"    "<<temp_dn_00.imag()<<"    "<<temp_dn_11.imag()<<"    "<<temp_dn_22.imag()<<"    "<<endl;

        }
        file_Akw_out<<endl;
    }



}


void Observables::Calculate_Nw(){

    //---------Read from input file-----------------------//
    string fileout="Nw.txt";
    double omega_min, omega_max, d_omega;
    double eta = 0.005;
    omega_min=1.0;omega_max=3.0;d_omega=0.00025;
    //---------------------------------------------------//

    int omega_index_max = int( (omega_max - omega_min)/(d_omega) );

    ofstream file_Nw_out(fileout.c_str());

    complex<double> temp_val11, temp_val22, temp_val00 ;
    int c1;


    for(int omega_ind=0;omega_ind<omega_index_max;omega_ind++){

      temp_val11=zero_complex;
      temp_val22=zero_complex;
      temp_val00=zero_complex;
      //l + 2*ns_ + ns_*orbs_*spin

        for (int j=0;j<Parameters_.ns;j++){
            for(int n=0;n<Hamiltonian_.Ham_.n_row();n++){

            c1 = j + ns_*Parameters_.orbs;
            temp_val00 +=  conj(Hamiltonian_.Ham_(c1,n))*Hamiltonian_.Ham_(c1,n)*
                    Lorentzian( omega_min + (omega_ind*d_omega) - Hamiltonian_.eigs_[n], eta);

            c1 = j+ ns_+ ns_*Parameters_.orbs;
            temp_val11 +=  conj(Hamiltonian_.Ham_(c1,n))*Hamiltonian_.Ham_(c1,n)*
                    Lorentzian( omega_min + (omega_ind*d_omega) - Hamiltonian_.eigs_[n], eta);

            c1 = j + 2*ns_+ ns_*Parameters_.orbs;
            temp_val22 +=  conj(Hamiltonian_.Ham_(c1,n))*Hamiltonian_.Ham_(c1,n)*
                    Lorentzian( omega_min + (omega_ind*d_omega) - Hamiltonian_.eigs_[n], eta);



            c1 = j;
            temp_val00 +=  conj(Hamiltonian_.Ham_(c1,n))*Hamiltonian_.Ham_(c1,n)*
                    Lorentzian( omega_min + (omega_ind*d_omega) - Hamiltonian_.eigs_[n], eta);

            c1 = j+ ns_;
            temp_val11 +=  conj(Hamiltonian_.Ham_(c1,n))*Hamiltonian_.Ham_(c1,n)*
                    Lorentzian( omega_min + (omega_ind*d_omega) - Hamiltonian_.eigs_[n], eta);

            c1 = j + 2*ns_;
            temp_val22 +=  conj(Hamiltonian_.Ham_(c1,n))*Hamiltonian_.Ham_(c1,n)*
                    Lorentzian( omega_min + (omega_ind*d_omega) - Hamiltonian_.eigs_[n], eta);


               }
        }

        file_Nw_out<<omega_min + (omega_ind*d_omega)<<"     "<<temp_val00.real()<<"     "
                  <<temp_val11.real()<<"     "<<temp_val22.real()<<endl;

    }


}

void Observables::Get_Non_Interacting_dispersion(){

    complex<double> one(1,0);
    complex<double> iota(0,1);
    complex<double> t1_,t2_,t3_,t4_,t5_,t6_,t7_,t8_;
    complex<double> Delta_xy_;
    Matrix<complex<double>> H;
    H.resize(3,3);

    char option = 'V';

    /*
    t1_ = 0.02*one;  t2_ = 0.06*one;
    t3_ = 0.03*one;  t4_ = -0.01*one;
    t5_ = 0.2*one;   t6_ = 0.3*one;
    */

    t1_ = -0.02*one;  t2_ = -0.06*one;
    t3_ = -0.03*one;  t4_ = 0.01*one;
    t5_ = -0.2*one;   t6_ = -0.3*one;
    t7_ = 0.2*one;  t8_ = 0.1*one;






    Delta_xy_=0.4*one;


    string fileout="k_vs_E_NI.txt";

    int counter_k=0;
    double dk_ =0.01;
    ofstream file_out(fileout.c_str());



    double kx_=0.0;
    double ky_=0.0;

    while(kx_<PI){

        H(0,0)=one*(-2.0*t2_*cos(kx_)  -2.0*t1_*cos(ky_)  -4.0*t3_*cos(kx_)*cos(ky_));
        H(1,1)=one*(-2.0*t1_*cos(kx_) -2.0*t2_*cos(ky_) - 4.0*t3_*cos(kx_)*cos(ky_));
        H(2,2)=one*(-2.0*t5_*(cos(kx_) + cos(ky_)) - 4.0*t6_*cos(kx_)*cos(ky_) + Delta_xy_);
        H(0,1)=one*(-4.0*t4_*sin(kx_)*sin(ky_));

        H(0,2)=iota*(-2.0*t7_*sin(kx_) + 4.0*t8_*sin(kx_)*cos(ky_) );
        H(1,2)=iota*(-2.0*t7_*sin(ky_) + 4.0*t8_*sin(ky_)*cos(kx_) );

        /*
        H(0,2)=(2.0*t7_*cos(kx_) + 4.0*iota*t8_*sin(kx_)*cos(ky_) );
       H(1,2)=(2.0*t7_*cos(ky_) + 4.0*iota*t8_*sin(ky_)*cos(kx_) );
       */






        char jobz=option;
        char uplo='U';
        int n=H.n_row();
        int lda=H.n_col();
        vector<complex<double>> work(3);
        vector<double> rwork(3*n);
        int info,lwork= -1;

        vector<double> E_;
        E_.resize(3);
        fill(E_.begin(), E_.end(),0);
        // query:
        zheev_(&jobz,&uplo,&n,&(H(0,0)),&lda,&(E_[0]),&(work[0]),&lwork,&(rwork[0]),&info);
        lwork = int(real(work[0]))+1;
        work.resize(lwork+1);
        // real work:
        zheev_(&jobz,&uplo,&n,&(H(0,0)),&lda,&(E_[0]),&(work[0]),&lwork,&(rwork[0]),&info);
        if (info!=0) {
            std::cerr<<"info="<<info<<"\n";
            perror("diag: zheev: failed with info!=0.\n");
        }

        file_out<<counter_k<<"    "<<kx_<<"   "<<ky_<<"    "<<E_[0]<<"     "<<E_[1]<<"     "<<E_[2]<<endl;
        kx_= kx_ + dk_;
        counter_k++;
    }


    kx_=PI;
    ky_=0;
    while(ky_<PI){

        H(0,0)=one*(-2.0*t2_*cos(kx_)  -2.0*t1_*cos(ky_)  -4.0*t3_*cos(kx_)*cos(ky_));
        H(1,1)=one*(-2.0*t1_*cos(kx_) -2.0*t2_*cos(ky_) - 4.0*t3_*cos(kx_)*cos(ky_));
        H(2,2)=one*(-2.0*t5_*(cos(kx_) + cos(ky_)) - 4.0*t6_*cos(kx_)*cos(ky_) + Delta_xy_);
        H(0,1)=one*(-4.0*t4_*sin(kx_)*sin(ky_));

        H(0,2)=iota*(-2.0*t7_*sin(kx_) + 4.0*t8_*sin(kx_)*cos(ky_) );
        H(1,2)=iota*(-2.0*t7_*sin(ky_) + 4.0*t8_*sin(ky_)*cos(kx_) );

        char jobz=option;
        char uplo='U';
        int n=H.n_row();
        int lda=H.n_col();
        vector<complex<double>> work(3);
        vector<double> rwork(3*n);
        int info,lwork= -1;

        vector<double> E_;
        E_.resize(3);
        fill(E_.begin(), E_.end(),0);
        // query:
        zheev_(&jobz,&uplo,&n,&(H(0,0)),&lda,&(E_[0]),&(work[0]),&lwork,&(rwork[0]),&info);
        lwork = int(real(work[0]))+1;
        work.resize(lwork+1);
        // real work:
        zheev_(&jobz,&uplo,&n,&(H(0,0)),&lda,&(E_[0]),&(work[0]),&lwork,&(rwork[0]),&info);
        if (info!=0) {
            std::cerr<<"info="<<info<<"\n";
            perror("diag: zheev: failed with info!=0.\n");
        }

        file_out<<counter_k<<"    "<<kx_<<"   "<<ky_<<"    "<<E_[0]<<"     "<<E_[1]<<"     "<<E_[2]<<endl;
        ky_= ky_ + dk_;
        counter_k++;
    }


    kx_=PI;
    ky_=PI;
    while(ky_>=0){

        kx_=ky_;

        H(0,0)=one*(-2.0*t2_*cos(kx_)  -2.0*t1_*cos(ky_)  -4.0*t3_*cos(kx_)*cos(ky_));
        H(1,1)=one*(-2.0*t1_*cos(kx_) -2.0*t2_*cos(ky_) - 4.0*t3_*cos(kx_)*cos(ky_));
        H(2,2)=one*(-2.0*t5_*(cos(kx_) + cos(ky_)) - 4.0*t6_*cos(kx_)*cos(ky_) + Delta_xy_);
        H(0,1)=one*(-4.0*t4_*sin(kx_)*sin(ky_));

        H(0,2)=iota*(-2.0*t7_*sin(kx_) + 4.0*t8_*sin(kx_)*cos(ky_) );
        H(1,2)=iota*(-2.0*t7_*sin(ky_) + 4.0*t8_*sin(ky_)*cos(kx_) );

        char jobz=option;
        char uplo='U';
        int n=H.n_row();
        int lda=H.n_col();
        vector<complex<double>> work(3);
        vector<double> rwork(3*n);
        int info,lwork= -1;

        vector<double> E_;
        E_.resize(3);
        fill(E_.begin(), E_.end(),0);
        // query:
        zheev_(&jobz,&uplo,&n,&(H(0,0)),&lda,&(E_[0]),&(work[0]),&lwork,&(rwork[0]),&info);
        lwork = int(real(work[0]))+1;
        work.resize(lwork+1);
        // real work:
        zheev_(&jobz,&uplo,&n,&(H(0,0)),&lda,&(E_[0]),&(work[0]),&lwork,&(rwork[0]),&info);
        if (info!=0) {
            std::cerr<<"info="<<info<<"\n";
            perror("diag: zheev: failed with info!=0.\n");
        }

        file_out<<counter_k<<"    "<<kx_<<"   "<<ky_<<"    "<<E_[0]<<"     "<<E_[1]<<"     "<<E_[2]<<endl;
        ky_= ky_ - dk_;
        counter_k++;
    }


}


double Observables::Lorentzian(double x, double brd){
    double temp;

    temp = (1.0/PI)*( (brd/2.0)/ ( (x*x) + ((brd*brd)/4.0) ) );

    return temp;

}

void Observables::DensityOfStates(){
    //-----------Calculate Bandwidth------//
    BandWidth=2.0;
    //-----------------------------------//

} // ----------


void Observables::OccDensity(){

} // ----------


void Observables::TotalOccDensity(){

} // ----------



complex<double> Observables::SiSjQ(int i,int j){return SiSjQ_(i,j);}

double Observables::SiSj(int i,int j){return SiSj_(i,j);}


complex<double> Observables::SiSjQ_Mean(int i,int j){return SiSjQ_Mean_(i,j);}

complex<double> Observables::SiSjQ_square_Mean(int i,int j){return SiSjQ_square_Mean_(i,j);}


double Observables::SiSj_Mean(int i,int j){return SiSj_Mean_(i,j);}

double Observables::SiSj_square_Mean(int i,int j){return SiSj_square_Mean_(i,j);}


void Observables::SiSjFULL(){

    double Cos_ij,Sin_ij,ei,ai,phase;
    int site_,site_p,ax,ay;


    for(int i=0;i<lx_;i++){
        for(int j=0;j<ly_;j++){
            site_ = Coordinates_.Nc(i,j);
            ei=MFParams_.etheta(i,j);
            ai=MFParams_.ephi(i,j);
            sx_[site_] = cos(ai) * sin(ei);
            sy_[site_] = sin(ai) * sin(ei);
            sz_[site_] = cos(ei);
        }
    }

    for(int xr=0;xr<lx_;xr++){
        for(int yr=0;yr<ly_;yr++){
            SiSj_(xr,yr)=double(0.0);
            for(int i=0;i<lx_;i++){
                for(int j=0;j<ly_;j++){
                    site_ = Coordinates_.Nc(i,j);
                    ax = (i+xr)%lx_;
                    ay = (j+yr)%ly_;
                    site_p = Coordinates_.Nc(ax,ay);
                    SiSj_(xr,yr) += sx_[site_]*sx_[site_p];
                    SiSj_(xr,yr) += sy_[site_]*sy_[site_p];
                    SiSj_(xr,yr) += sz_[site_]*sz_[site_p];
                }
            }
            SiSj_(xr,yr)*= double(1.0/(lx_*ly_));
            //cout << xr << " "<< yr<< " "<<  SiSj_(xr,yr) << endl;
        }
    }

    for(int qx=0; qx<lx_; qx++) {
        for(int qy=0; qy<ly_; qy++) {
            SiSjQ_(qx,qy)=complex<double>(0.0,0.0);
            for(int xr=0;xr<lx_;xr++){
                for(int yr=0;yr<ly_;yr++){
                    phase=2.0*Parameters_.pi*(double(qx*xr)/double(lx_)+double(qy*yr)/double(ly_));
                    Cos_ij = cos(phase);
                    Sin_ij = sin(phase);
                    SiSjQ_(qx,qy) += SiSj_(xr,yr)*complex<double>(Cos_ij,Sin_ij);
                }
            }
            SiSjQ_(qx,qy)*= double(1.0/(lx_*ly_));
            //cout << qx << " "<< qy<< " "<<  SiSjQ_(qx,qy) << endl;
        }
    }

    //     cout << 0 << " "<< 1 << " "<<  SiSj_(0,1) << endl;
    //     cout << 1 << " "<< 0 << " "<<  SiSj_(1,0) << endl;
    //     cout << 0 << " "<< 4 << " "<<  SiSjQ_(0,4) << endl;
    //     cout << 4 << " "<< 0 << " "<<  SiSjQ_(4,0) << endl;
    //     cout << 2 << " "<< 6 << " "<<  SiSjQ_(2,6) << endl;
    //     cout << 6 << " "<< 2 << " "<<  SiSjQ_(6,2) << endl;

} // ----------


void Observables::SiSjQ_Average(){

    for(int qx=0; qx<lx_; qx++) {
        for(int qy=0; qy<ly_; qy++) {
            SiSjQ_Mean_(qx,qy) += SiSjQ_(qx,qy);
            SiSjQ_square_Mean_(qx,qy) += ( SiSjQ_(qx,qy)*SiSjQ_(qx,qy) )   ;
            //cout << qx << " "<< qy<< " "<<  SiSjQ_(qx,qy) << endl;
        }
    }

    //     cout << 0 << " "<< 1 << " "<<  SiSj_(0,1) << endl;
    //     cout << 1 << " "<< 0 << " "<<  SiSj_(1,0) << endl;
    //     cout << 0 << " "<< 4 << " "<<  SiSjQ_(0,4) << endl;
    //     cout << 4 << " "<< 0 << " "<<  SiSjQ_(4,0) << endl;
    //     cout << 2 << " "<< 6 << " "<<  SiSjQ_(2,6) << endl;
    //     cout << 6 << " "<< 2 << " "<<  SiSjQ_(6,2) << endl;

} // ----------



void Observables::SiSj_Average(){

    for(int x=0; x<lx_; x++) {
        for(int y=0; y<ly_; y++) {
            SiSj_Mean_(x,y) += SiSj_(x,y);
            SiSj_square_Mean_(x,y) += ( SiSj_(x,y)*SiSj_(x,y) )   ;
            //cout << qx << " "<< qy<< " "<<  SiSjQ_(qx,qy) << endl;
        }
    }

    Nematic_order_mean_ += fabs(SiSj_(1,0) - SiSj_(0,1))*0.5;
    Nematic_order_square_mean_ += (SiSj_(1,0) - SiSj_(0,1) )*(SiSj_(1,0) - SiSj_(0,1) )*0.25;

} // ----------


void Observables::Total_Energy_Average(double Curr_QuantE, double CurrE){

    AVG_Total_Energy += Curr_QuantE + CurrE;
    AVG_Total_Energy_sqr += (Curr_QuantE + CurrE)*(Curr_QuantE + CurrE);
}

void Observables::OccDensity(int tlabel){

} // ----------


void Observables::DOSprint(int tlabel){
    double omega;
    //create name
    std::string name="Output/OrbDOS_" + to_string(tlabel) + ".dat";
    ofstream myfile;
    myfile.open(name);
    for(int ll=0;ll<=800;ll++) {
        omega=Omega(ll);
        myfile << omega-Parameters_.mus << "\t"
               << setw(12) << dos(0,ll)/(ns_*Parameters_.MCNorm) << "\t"
               << setw(12) << dos(1,ll)/(ns_*Parameters_.MCNorm) << "\t"
               << setw(12) << dos(2,ll)/(ns_*Parameters_.MCNorm) << "\t"
               << setw(12) << dos(3,ll)/(ns_*Parameters_.MCNorm) << "\t"
               << endl;
    }
    myfile.close();
} // ----------



void Observables::Initialize(){

    complex<double> zero(0.0,0.0);
    int space=2*Parameters_.orbs*ns_;
    sx_.resize(space);
    sy_.resize(space);
    sz_.resize(space);

    dos.resize(4,801);   dos.fill(0.0);

    Nematic_order_mean_ =0.0;
    Nematic_order_square_mean_ =0.0;

    SiSj_.resize(lx_,ly_);
    SiSj_Mean_.resize(lx_,ly_);
    SiSj_square_Mean_.resize(lx_,ly_);

    SiSjQ_Mean_.resize(lx_,ly_);
    SiSjQ_square_Mean_.resize(lx_,ly_);
    SiSjQ_.resize(lx_,ly_);

    for(int ix=0;ix<lx_;ix++){
        for(int iy=0;iy<ly_;iy++){
            SiSjQ_Mean_(ix,iy)=zero;
            SiSjQ_square_Mean_(ix,iy)=zero;
        }

    }
    nia_.resize(ns_); std::fill(nia_.begin(),nia_.end(),0.0);
    nib_.resize(ns_); std::fill(nib_.begin(),nib_.end(),0.0);
    nic_.resize(ns_); std::fill(nic_.begin(),nic_.end(),0.0);

    nia_t=0.0; nib_t=0.0; nic_t=0.0;
    dosincr_=0.05;
    tpi_=4.0f*atan(1.0f);
} // ----------


double Observables::Omega(int i){
    return -20.0+double(i)*dosincr_;
} // ----------









#endif // OBSERVABLES_H
