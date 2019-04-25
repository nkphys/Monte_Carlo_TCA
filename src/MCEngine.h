#include "ParametersEngine.h"
#include "Coordinates.h"
#include "MFParams.h"
#include "Hamiltonian.h"
#include "Observables.h"
#include "tensor_type.h"

#ifndef MCENGINE_H
#define MCENGINE_H


class MCEngine{
public:
    MCEngine(Parameters& Parameters__, Coordinates& Coordinates__,
             MFParams& MFParams__, Hamiltonian& Hamiltonian__,
             Observables& Observables__)
        : Parameters_(Parameters__),Coordinates_(Coordinates__),
          MFParams_(MFParams__), Hamiltonian_(Hamiltonian__),
          Observables_(Observables__),
          lx_(Parameters_.lx), ly_(Parameters_.ly), ns_(Parameters_.ns),
          orbs_(Parameters_.orbs), ED_(Parameters_.ED_)
    {

    }

    void RUN_MC();
    double Prob (double muu, double mu_new);
    double ProbCluster (double muu, double mu_new);
    Parameters &Parameters_;
    Coordinates &Coordinates_;
    MFParams &MFParams_;
    Hamiltonian &Hamiltonian_;
    Observables &Observables_;
    const int lx_, ly_, ns_, orbs_;
    bool ED_;

};

/*
 * ***********
 *  Functions in Class MCEngine ------
 *  ***********
*/

void MCEngine::RUN_MC(){

    complex<double> zero(0.0,0.0);
    bool Metropolis_Algo = Parameters_.Metropolis_Algorithm;
    bool Heat_Bath_Algo = Parameters_.Heat_Bath_Algorithm;

    int MC_sweeps_used_for_Avg=Parameters_.Last_n_sweeps_for_measurement;
    int Gap_bw_sweeps = Parameters_.Measurement_after_each_m_sweeps;

    double PrevE,CurrE,P_new,P12,muu_prev;
    double muu_prevCluster;
    double Curr_QuantE;
    double Prev_QuantE;
    double Curr_QuantECluster;
    double Prev_QuantECluster;
    int x,y,act;
    double saved_Params[2];

    string File_Out_progress;
    string File_Out_theta_phi;

    double temp_=Parameters_.temp_max;

    double initial_mu_guess;
    int n_states_occupied_zeroT;



    double Curr_Cluster_CLE;

    //starting with a random guess

    while(temp_>=Parameters_.temp_min){

        cout << "Temperature = " << temp_<<" is being done"<<endl;
        Parameters_.temp=temp_;
        Parameters_.beta=double(11604.0/temp_);

        for(int ix=0;ix<lx_;ix++){
            for(int iy=0;iy<ly_;iy++){
                Observables_.SiSjQ_Mean_(ix,iy)=zero;
                Observables_.SiSjQ_square_Mean_(ix,iy)=zero;
                Observables_.SiSj_square_Mean_(ix,iy)=0.0;
                Observables_.SiSj_Mean_(ix,iy)=0.0;

            }
        }
        Observables_.AVG_Total_Energy=0.0;
        Observables_.AVG_Total_Energy_sqr=0.0;
        Observables_.Nematic_order_square_mean_=0.0;
        Observables_.Nematic_order_mean_=0.0;

        MFParams_.etheta_avg.fill(0.0);
        MFParams_.ephi_avg.fill(0.0);



        char temp_char[50];
        sprintf(temp_char,"%.1f",temp_);

        File_Out_progress = "output_Temp" + string(temp_char) + ".txt";
        ofstream file_out_progress(File_Out_progress.c_str());

        File_Out_theta_phi = "ThetaPhi_Temp" + string(temp_char) + ".txt";
        ofstream File_Out_Theta_Phi(File_Out_theta_phi.c_str());



        file_out_progress<< "Total "<<Parameters_.IterMax<<" sweeps are performed."<<endl;
        file_out_progress<<"First "<<Parameters_.IterMax - (Gap_bw_sweeps*(MC_sweeps_used_for_Avg - 1) + MC_sweeps_used_for_Avg)<<
                           " sweeps are used for thermalization and every "<<Gap_bw_sweeps+1<<" in last "<<
                           Gap_bw_sweeps*(MC_sweeps_used_for_Avg - 1) + MC_sweeps_used_for_Avg<<
                           " sweeps are used for measurement."<<endl;
        act=1;




        Parameters_.WindowSize = 0.1; //2f + 0.003f*beta0 ;
        Parameters_.Eav=0.0;
        Parameters_.MCNorm=0;
        Parameters_.Dflag='N'; // flag to calculate only Eigenvalue
        //std::string name="Output/Conf_" + to_string(ltemp) + ".dat";
        //Parameters_.beta = double(11604.0/ (Parameters_.temp +20.0) );
        //cout << "TEMP  " << Parameters_.temp << endl;


        file_out_progress<<"I_MC"<<setw(15)<<"S(0,1)"<<setw(15)<<"S(1,0)"
                        <<setw(15)<<"S(0,Pi)"<<setw(15)<<"S(Pi,0)"<<setw(17)<<"S(0,0)"<<setw(17)<<"S(Pi,Pi)"<<setw(17)<<
                          "S(Pi/2,Pi/2)"<<setw(17)<<"< N_total >"
                       <<setw(15)<<"E_CL"<<setw(15)<<"E_QM"<<setw(15)<<"mu"<< endl;




        PrevE = Hamiltonian_.GetCLEnergy();
        Hamiltonian_.InteractionsCreate();
        Hamiltonian_.Diagonalize(Parameters_.Dflag);
        n_states_occupied_zeroT=Parameters_.ns*Parameters_.Fill*Parameters_.orbs*2.0;
        initial_mu_guess=0.5*(Hamiltonian_.eigs_[n_states_occupied_zeroT-1] + Hamiltonian_.eigs_[n_states_occupied_zeroT]);
        //initial_mu_guess=0.25;
        Parameters_.mus=Hamiltonian_.chemicalpotential(initial_mu_guess,Parameters_.Fill);
        Prev_QuantE = Hamiltonian_.E_QM();
        muu_prev=Parameters_.mus;
        Hamiltonian_.copy_eigs(1);
        cout<<"Initial Classical Energy[Full System] = "<<PrevE<<endl;
        cout<<"Initial Quantum Energy[Full System] = "<<Prev_QuantE<<endl;
        cout<<"Initial Total Energy[Full System] = "<<PrevE+Prev_QuantE<<endl;
        cout<<"Initial mu="<<muu_prev<<endl;
        //for(int i=0;i<10;i++){
        //  cout<<i<<"   "<<Hamiltonian_.eigs_[i]<<endl;
        // }



        int Confs_used=0;
        int measure_start=0;
        muu_prevCluster=muu_prev;

        if(ED_){
            Prev_QuantECluster=Prev_QuantE;
            Hamiltonian_.eigsCluster_saved_=Hamiltonian_.eigs_saved_;
        }

        for(int count=0;count<Parameters_.IterMax;count++){
            //if (count == 1){
            // Parameters_.beta = double(11604.0/ Parameters_.temp);
            // PrevE = Hamiltonian_.GetCLEnergy();
            //  Hamiltonian_.InteractionsCreate();
            //  Hamiltonian_.Diagonalize(Parameters_.Dflag);
            //  Hamiltonian_.copy_eigs(1);
            //muu = Hamiltonian_.chemicalpotential(Parameters_.mus,Parameters_.Fill);
            //Parameters_.mus = Parameters_.mus*0.4f + muu*0.6f;
            // }

            for(int i=0;i<ns_;i++) {  // For each site


                //***Before change*************//

                if(ED_==false){
                    //TCA is used
                    PrevE = Hamiltonian_.GetCLEnergy();

                    if(Parameters_.J_HUND !=0.0){
                    Hamiltonian_.InteractionsClusterCreate(i);
                    Hamiltonian_.DiagonalizeCluster(Parameters_.Dflag);
                    //n_states_occupied_zeroT=Parameters_.Fill*Hamiltonian_.eigsCluster_.size();
                    //initial_mu_guess=0.5*(Hamiltonian_.eigsCluster_[n_states_occupied_zeroT-1] + HamiltonianCluster_.eigs_[n_states_occupied_zeroT])
                    muu_prevCluster=Hamiltonian_.chemicalpotentialCluster(muu_prevCluster,Parameters_.Fill);
                    Prev_QuantECluster = Hamiltonian_.E_QMCluster();
                    Hamiltonian_.copy_eigs_Cluster(1);
                    }
                    else{
                        assert(Parameters_.J_HUND==0.0);
                        Parameters_.mus_Cluster=0.0;
                        Curr_QuantECluster=0.0;
                    }
                }
                else{
                assert(ED_);
                }


                //*******************************//



                x=Coordinates_.indx(i);
                y=Coordinates_.indy(i);
                saved_Params[0]=MFParams_.etheta(x,y);
                saved_Params[1]=MFParams_.ephi(x,y);

                MFParams_.FieldThrow(i);
                CurrE = Hamiltonian_.GetCLEnergy();

                if(Parameters_.J_HUND !=0.0){
                Hamiltonian_.InteractionsClusterCreate(i);
                Hamiltonian_.DiagonalizeCluster(Parameters_.Dflag);
                Parameters_.mus_Cluster=Hamiltonian_.chemicalpotentialCluster(muu_prevCluster,Parameters_.Fill);    
                Curr_QuantECluster = Hamiltonian_.E_QMCluster();
                }
                else{
                    assert(Parameters_.J_HUND==0.0);
                    Parameters_.mus_Cluster=0.0;
                    Curr_QuantECluster=0.0;
                }

                //Ratio of Quantum partition functions
                /*P = [ Tr(exp(-beta(Hquant_new)))/Tr(exp(-beta(Hquant_old)))]*
                      [exp(-beta*E_classical(New)) / exp(-beta*E_classical(old))]
                     * [sin(Theta_i(New)) / sin(Theta_i(Old)) ]*/
                /*exp(P12) = P
                  P12 = log (P)
                  */

                //same mu-refrence is used, otherwise engine does not work properly
                P_new = ProbCluster(muu_prev, muu_prev);
                P12 = P_new - Parameters_.beta*((CurrE)-(PrevE));
                //P12 = - Parameters_.beta*((CurrE)-(PrevE));
                //cout<<P12<<endl;
                P12 += log ((sin(MFParams_.etheta(x,y))/sin(saved_Params[0])));



                //---OR---
                //P12 = exp(-Parameters_.beta*((CurrE+Curr_QuantE)-(PrevE+Prev_QuantE)));
                //P12*= (sin(MFParams_.etheta(x,y))/sin(saved_Params[0]));

                //Heat bath algorithm [See page-129 of Prof. Elbio's Book]
                //Heat bath algorithm works for small changes i.e. when P12~1.0
                //  if (Heat_Bath_Algo){
                //     P12 =P12/(1.0+P12);
                //  }


                //Metropolis Algotithm
                // if (Metropolis_Algo){
                //    P12=min(1.0,P12);
                // }




                /*
       * VON NEUMANN's REJECTING METHOD:
       * Random number < P12 -----> ACCEPT
       * Random number > P12 -----> REJECT
       */


                //ACCEPTED
                if(P12 > 0){
                    Parameters_.AccCount[0]++;
                    act=1;
                    if(ED_){
                    PrevE=CurrE;
                    Prev_QuantECluster = Curr_QuantECluster;
                    Hamiltonian_.copy_eigs_Cluster(1);
                    muu_prevCluster=Parameters_.mus_Cluster;
                    }

                }
                else if ( exp(P12) > ( 1.0 - MFParams_.random() ) ) {
                    Parameters_.AccCount[0]++;
                    act=1;
                    if(ED_){
                    PrevE=CurrE;
                    Prev_QuantECluster = Curr_QuantECluster;
                    Hamiltonian_.copy_eigs_Cluster(1);
                    muu_prevCluster=Parameters_.mus_Cluster;
                    }

                }

                //REJECTED
                else{
                    Parameters_.AccCount[1]++;
                    act=0;
                    MFParams_.etheta(x,y) = saved_Params[0];
                    MFParams_.ephi(x,y)   = saved_Params[1];
                }


                // if ((act == 1) && (count<1000)) {

                //muu = Hamiltonian_.chemicalpotential(Parameters_.mus,Parameters_.Fill);
                //Parameters_.mus = Parameters_.mus*0.999 + muu*0.001;
                //Parameters_.mus = muu;
                //}

            }// site loop

            //      if (act == 1) {
            //       muu = Hamiltonian_.chemicalpotential(Parameters_.mus,Parameters_.Fill);
            //       Parameters_.mus = Parameters_.mus*0.99f + muu*0.01f;
            //      }

            if ( (count%10==0) ) {
                MFParams_.Adjust_MCWindow();
            }

            if(count < (Parameters_.IterMax - (Gap_bw_sweeps*(MC_sweeps_used_for_Avg - 1) + MC_sweeps_used_for_Avg)) ){
                if ( (count%10==0) ) {
                    Observables_.SiSjFULL();
                    file_out_progress << int(1.0*count) <<setw(20)<<Observables_.SiSj(0,1) <<setw(16)<< Observables_.SiSj(1,0)
                                      <<setw(16)<< Observables_.SiSjQ(0,int(lx_/2)).real() <<setw(16)<< Observables_.SiSjQ(int(lx_/2),0).real()
                                     <<setw(16)<<Observables_.SiSjQ(0,0).real() <<setw(16)<< Observables_.SiSjQ(int(lx_/2),int(lx_/2)).real() <<setw(16)<<
                                       Observables_.SiSjQ(int(lx_/4),int(lx_/4)).real() <<setw(16)<< Hamiltonian_.ClusterDensity() <<setw(16)<< CurrE
                                    <<setw(16)<< Curr_QuantECluster<<setw(15)<<Parameters_.mus_Cluster<< endl;
                }
            }
            //Average and Std. deviation is calculated is done
            else{




                if(measure_start==0){
                    measure_start++;
                    file_out_progress<<"----------Measurement is started----------"<<endl;
                    file_out_progress<<"I_MC      Avg{S(pi,0)}    Avg{S(0,pi)}    std.dev{S(pi,0)}   std.dev{S(0,pi)} Avg{S(1,0)}  Avg{S(0,1)}  std.dev{S(1,0)}  std.dev{S(0,1)}  Avg{Nematic_order}    std.dev{Nematic_order}     Avg{E_classical}  std.dev{E_classical}"<<endl;
                }
                int temp_count=count -
                        (Parameters_.IterMax - (Gap_bw_sweeps*(MC_sweeps_used_for_Avg - 1) + MC_sweeps_used_for_Avg));
                int zero_or_not = temp_count % (Gap_bw_sweeps + 1);
                if( zero_or_not==0 ){

                    if((Parameters_.Saving_Microscopic_States==true) &&
                       (Confs_used<Parameters_.No_Of_Microscopic_States)
                            ){

                        char Confs_char[50];
                        sprintf(Confs_char,"%d",Confs_used);
                        string File_Out_theta_phi_microState = "ThetaPhi_Temp" + string(temp_char) +
                                                        "MicroState" + string(Confs_char) +".txt";
                        ofstream File_Out_Theta_Phi_MicroState(File_Out_theta_phi_microState.c_str());


                        File_Out_Theta_Phi_MicroState<<"#x"<<setw(15)<<"y"<<setw(15)<<"Theta(x,y)"<<setw(15)<<"Phi(x,y)"<<endl;
                        for(int ix=0;ix<lx_;ix++){
                            for(int iy =0;iy<ly_;iy++){
                                File_Out_Theta_Phi_MicroState<<ix<<setw(15)<<iy<<setw(15)<<MFParams_.etheta(ix,iy)<<setw(15)<<MFParams_.ephi(ix,iy)<<endl;
                            }}

                    }


                    Confs_used=Confs_used+1;
                    Observables_.SiSjFULL();
                    Observables_.SiSjQ_Average();
                    Observables_.SiSj_Average();
                    //Just Classical Energy
                    Observables_.Total_Energy_Average(0.0, CurrE);

                    MFParams_.Calculate_Fields_Avg();

                    //double MC_steps_Avg_insitu = (1.0 + 1.0*(count - (Parameters_.IterMax - MC_steps_used_for_Avg)));

                    file_out_progress << int(1.0*count) <<setw(20)<<
                                         Observables_.SiSjQ_Mean(int(lx_/2),0).real()/(Confs_used*1.0)
                                      <<setw(16)<<Observables_.SiSjQ_Mean(0,int(lx_/2)).real()/(Confs_used*1.0)
                                     <<setw(16)<<


                                       sqrt(
                                           (( Observables_.SiSjQ_square_Mean(int(lx_/2),0)/(Confs_used*1.0) ) -
                                            ((Observables_.SiSjQ_Mean(int(lx_/2),0)*Observables_.SiSjQ_Mean(int(lx_/2),0) )/(Confs_used*Confs_used*1.0) ) ).real()
                                           )

                                    <<setw(16)<<
                                      sqrt(
                                          (( Observables_.SiSjQ_square_Mean(0,int(lx_/2))/(Confs_used*1.0) ) -
                                           ((Observables_.SiSjQ_Mean(0,int(lx_/2))*Observables_.SiSjQ_Mean(0,int(lx_/2)) )/(Confs_used*Confs_used*1.0) ) ).real()
                                          )
                                   <<setw(16)<<


                                     Observables_.SiSj_Mean(1,0)/(Confs_used*1.0)
                                  <<setw(16)<<Observables_.SiSj_Mean(0,1)/(Confs_used*1.0)
                                 <<setw(16)<<
                                   sqrt(
                                       (( Observables_.SiSj_square_Mean(1,0)/(Confs_used*1.0) ) -
                                        ((Observables_.SiSj_Mean(1,0)*Observables_.SiSj_Mean(1,0) )/(Confs_used*Confs_used*1.0) ) )
                                       )

                                <<setw(16)<<
                                  sqrt(
                                      (( Observables_.SiSj_square_Mean(0,1)/(Confs_used*1.0) ) -
                                       ((Observables_.SiSj_Mean(0,1)*Observables_.SiSj_Mean(0,1) )/(Confs_used*Confs_used*1.0) ) )
                                      )
                               <<setw(16)<<
                                 Observables_.Nematic_order_mean_/(Confs_used*1.0)

                              <<setw(16)<<
                                sqrt(
                                    (( Observables_.Nematic_order_square_mean_/(Confs_used*1.0) ) -
                                     ((Observables_.Nematic_order_mean_*Observables_.Nematic_order_mean_)/(Confs_used*Confs_used*1.0) ) )
                                    )


                             <<setw(32)<<
                               Observables_.AVG_Total_Energy/(Confs_used*1.0)
                            <<setw(16)<<
                              sqrt(  (Observables_.AVG_Total_Energy_sqr/(Confs_used*1.0)) -
                                     ((Observables_.AVG_Total_Energy*Observables_.AVG_Total_Energy)/(Confs_used*Confs_used*1.0))  )
                           <<endl;

                }

            }


        }// Iter Loop
        file_out_progress << "Total "<<Confs_used<< " configurations were used were measurement"<<endl;

        temp_ = temp_ - Parameters_.d_Temp;

        File_Out_Theta_Phi<<"#x"<<setw(15)<<"y"<<setw(15)<<"Theta_avg(x,y)"<<setw(15)<<"Phi_avg(x,y)"<<endl;
        for(int ix=0;ix<lx_;ix++){
            for(int iy =0;iy<ly_;iy++){
                File_Out_Theta_Phi<<ix<<setw(15)<<iy<<setw(15)<<MFParams_.etheta_avg(ix,iy)/(Confs_used*1.0)<<setw(15)<<MFParams_.ephi_avg(ix,iy)/(Confs_used*1.0)<<endl;
            }}


        MFParams_.Read_classical_DOFs(File_Out_theta_phi);
    }//Temperature loop

} // ---------



double MCEngine::Prob(double muu, double mu_new){

    double P=0.0;
    double X,Y,X2;

    for(int i=0;i<2*orbs_*ns_;i++){
        X = Parameters_.beta*( (mu_new) - Hamiltonian_.eigs_[i]);
        Y = Parameters_.beta*( (muu) - Hamiltonian_.eigs_saved_[i]);
        //P += log(1 + exp(X)) - log(1 + exp(Y));

        if(X>5){
            P +=X;
        }
        else if(fabs(X)<0.001){
            P += log(2.0 + X);
        }
        else if(X<-5){
            P +=exp(X);
        }
        else{
            P +=log(1.0 + exp(X));
        }



        if(Y>5){
            P -=Y;
        }
        else if(fabs(Y)<0.001){
            P -= log(2.0 + Y);
        }
        else if(Y<-5){
            P -=exp(Y);
        }
        else{
            P -=log(1.0 + exp(Y));
        }




    }

    return P;

} // ---------



double MCEngine::ProbCluster(double muu, double mu_new){

    double P=0.0;
    double X,Y,X2;
    int ns = (Parameters_.lx_cluster)*(Parameters_.ly_cluster);

    for(int i=0;i<2*orbs_*ns;i++){
        X = Parameters_.beta*( (mu_new) - Hamiltonian_.eigsCluster_[i]);
        Y = Parameters_.beta*( (muu) - Hamiltonian_.eigsCluster_saved_[i]);
        //P += log(1 + exp(X)) - log(1 + exp(Y));

        if(X>5){
            P +=X;
        }
        else if(fabs(X)<0.001){
            P += log(2.0 + X);
        }
        else if(X<-5){
            P +=exp(X);
        }
        else{
            P +=log(1.0 + exp(X));
        }



        if(Y>5){
            P -=Y;
        }
        else if(fabs(Y)<0.001){
            P -= log(2.0 + Y);
        }
        else if(Y<-5){
            P -=exp(Y);
        }
        else{
            P -=log(1.0 + exp(Y));
        }




    }

    return P;

} // ---------





#endif // MCENGINE_H
