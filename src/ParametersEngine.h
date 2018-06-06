#ifndef Parameters_class
#define Parameters_class

class Parameters{

public:
    int lx, ly, ns, orbs, IterMax, MCNorm, RandomSeed;
    int TBC_mx, TBC_my;
    int TBC_cellsX, TBC_cellsY;
    int lx_cluster, ly_cluster;
    double mus, mus_Cluster, Fill,pi,J_NN,J_NNN,J_HUND,k_const,lamda_12,lamda_66;

    bool Cooling_;

    bool Metropolis_Algorithm;
    bool Heat_Bath_Algorithm;

    double temp_max, beta_min;
    double temp_min, beta_max;
    double d_Temp;

    int Last_n_sweeps_for_measurement;
    int Measurement_after_each_m_sweeps;

    double temp,beta,Eav,maxmoment;
    double WindowSize, AccCount[2];
    char Dflag;

    void Initialize(string inputfile_);
    double matchstring(string file,string match);

};


void Parameters::Initialize(string inputfile_){

    maxmoment=1.0;
    double cooling_double;
    double metropolis_double;

    cout << "____________________________________" << endl;
    cout << "Reading the inputfile: " << inputfile_ << endl;
    cout << "____________________________________" << endl;
    lx = int(matchstring(inputfile_,"Xsite"));
    ly = int(matchstring(inputfile_,"Ysite"));
    TBC_mx = int(matchstring(inputfile_,"TwistedBoundaryCond_mx"));
    TBC_my = int(matchstring(inputfile_,"TwistedBoundaryCond_my"));
    TBC_cellsX = int(matchstring(inputfile_,"TBC_cellsX"));
    TBC_cellsY = int(matchstring(inputfile_,"TBC_cellsY"));
    lx_cluster = int(matchstring(inputfile_,"Cluster_lx"));
    ly_cluster = int(matchstring(inputfile_,"Cluster_ly"));


    ns = lx*ly;
    cout << "TotalNumberOfSites = "<< ns << endl;
    orbs = int(matchstring(inputfile_,"Orbitals"));
    Fill = matchstring(inputfile_,"Fill");
    cout << "TotalNumberOfParticles = "<< ns*Fill*orbs*2.0 << endl;

    IterMax = int(matchstring(inputfile_,"MaxMCsweeps"));
    MCNorm = 0.0; ; //matchstring(inputfile,"MCNorm")
    RandomSeed = matchstring(inputfile_,"RandomSeed");
    Dflag = 'N';
    J_NN=double(matchstring(inputfile_,"J_NN"));
    J_NNN=double(matchstring(inputfile_,"J_NNN"));
    J_HUND=double(matchstring(inputfile_,"J_HUND"));


    metropolis_double=double(matchstring(inputfile_,"Metropolis_Algo"));
    if(metropolis_double==1.0){
        Metropolis_Algorithm=true;
        Heat_Bath_Algorithm=false;
    }
    else if(metropolis_double==0.0){
        Metropolis_Algorithm=false;
        Heat_Bath_Algorithm=true;
    }
    else{
        cout<<"ERROR: Metropolis_Algo can be only 1 (true) or 0 (false)"<<endl;
        assert(metropolis_double==0.0);
    }

    cooling_double=double(matchstring(inputfile_,"Cooling"));
    if(cooling_double==1.0){
        Cooling_=true;

        temp_min = double(matchstring(inputfile_,"Temperature_min"));
        temp_max = double(matchstring(inputfile_,"Temperature_max"));
        d_Temp = double(matchstring(inputfile_,"dTemperature"));
        beta_max=double(11604.0/ temp_min);
        beta_min=double(11604.0/ temp_max);

    }
    else if(cooling_double==0.0){
        Cooling_=false;

        temp = double(matchstring(inputfile_,"Temperature"));   // temperature in kelvin
        beta=double(11604.0/ temp);    //Beta which is (T*k_b)^-1

        temp_min = temp;
        temp_max = temp;
        d_Temp=10.0;//arbitrary positive number

    }
    else{
        cout<<"ERROR: Cooling can be only 1 (true) or 0 (false)"<<endl;
        assert(cooling_double==0.0);
    }





    Last_n_sweeps_for_measurement=int(matchstring(inputfile_,"Last_n_sweeps_for_measurement"));
    Measurement_after_each_m_sweeps=int(matchstring(inputfile_,"Measurement_after_each_m_sweeps"));


    pi=4.00*atan(double(1.0));
    Eav=0.0;
    AccCount[0]=0;
    AccCount[1]=0;

    WindowSize=double(0.01);
    mus=0.25;
    cout << "____________________________________" << endl;
}


double Parameters::matchstring(string file,string match) {
    string test;
    string line;
    ifstream readFile(file);
    double amount;
    bool pass=false;
    while (std::getline(readFile, line)) {
        std::istringstream iss(line);
        if (std::getline(iss, test, '=') && pass==false) {
            // ---------------------------------
            if (iss >> amount && test==match) {
                // cout << amount << endl;
                pass=true;
            }
            else {
                pass=false;
            }
            // ---------------------------------
            if(pass) break;
        }
    }
    if (pass==false) {
        string errorout=match;
        errorout+="= argument is missing in the input file!";
        throw std::invalid_argument(errorout);
    }
    cout << match << " = " << amount << endl;
    return amount;
}

#endif



