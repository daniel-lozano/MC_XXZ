#include <iostream>
#include <cmath>
#include <fstream>
#include <time.h>
#include <string.h>
using namespace std;


void sweep(int* Spins, int** neighbours,double T, int N_spins  );
double get_Energy(int* Spins, int** neighbours, int N , double Temp );
double get_Cv_correction(int* Spins, int** neighbours, int N , double Temp );
double get_magnetization(int* Spins, int N);

double get_magnetization_A(int* Spins, int N);
double get_magnetization_B(int* Spins, int N);

double random_double();
double X(int* Spins,int** neighbours, int i,int j,double Temp,int sign);
double F(int* Spins,int** neighbours, int i,int j,double Temp, int sign);
double H(int* Spins,int** neighbours, int i,int j,double Temp,int sign);



#define J 1. // Positive for AFM interactions
#define Jperp 0.1
#define KB 1.
#define Z 4

int main(int argc, char** argv){
    
    /// Constants of the model ///
    
    
    /// Lattice constants ///
    int SQU=atoi(argv[1]);
    int L,H;
    L=SQU;
    H=SQU;
    int N_spins=L*H;
    
    double Tc=2.0/log(1+sqrt(2))*J; //
    
    int n_eqSweeps=2000;
    int n_measSweeps=7000;
  
    cout << "The size of the lattice is "<< SQU << endl;
    cout << "The number of equillibrium sweeps= "<< n_eqSweeps<<endl;
    cout << "The number of measureing sweeps= "<< n_measSweeps<<endl;

    
    /// Changin seed to a time dependend seed ///
    
    cout << "N_spins=" << N_spins << endl;
    
    ////------------------------------------------------------------------
    ///----------------- Creating array of neighbours ---------------- ///
    ////------------------------------------------------------------------
    
    int** neighbours =new int*[N_spins];
    for(int i=0;i<N_spins;i++){
        neighbours[i]=new int[Z];
    }
    
    /// Fealing the neighbour matrix ///
    for(int i=0;i<L;i++){
        for(int j=0;j<H;j++){
            neighbours[i*L+j][0]=i*L+(j+1)%H; //Right neighbour
            
            neighbours[i*L+j][1]=((i-1+L)%L)*L+j; //Upwards neighbour
            if(i==0){neighbours[i*L+j][1]=(L-1)*L+j; }
            
            neighbours[i*L+j][2]=i*L+(j-1+H)%H; //Left neighbour
            if(j==0){neighbours[i*L+j][2]=i*L+H-1; }
            
            neighbours[i*L+j][3]=((i+1)%L)*L+j; //Downwards neighbour
        }
    }
    
 
    
    int N_temps=50;
    double T_max=4;
    double T_min=1;
    double dt=(log(T_max/T_min))/N_temps;
    double* Temp=new double[N_temps];
    
//    Temp[0]=T_max;
//    for(int i=1; i<N_temps;i++){
//        Temp[i]=Temp[i-1]*exp(-dt);// Done in log scale for ploting
//
//    }
    
    for(int i=0; i<N_temps;i++){
        Temp[i]=T_max-i*(T_max-T_min)/N_temps;
    }
    /// Files to be use ///
    ofstream File_E_M,File_Magnettization_powers;
    
    File_E_M.open("Energy_magnetization_L"+ to_string(L)+".txt");
    File_Magnettization_powers.open("Magnetization.txt",ios::app);
  
 
    /// Starting with a fix spin configuration ///
    
    /// Random Spin configurations ////
    int* Spins=new int[N_spins];
    
    srand(time(NULL));

    for(int i=0;i<N_spins;i++){
        Spins[i]=2*(rand()%2)-1;
       
//        cout<< Spins[i] << " ";
//        if((i+1)%L==0){
//            cout<<endl;
//        }

    }
//    cout<<endl;
//    int pos1=L+int(L/2);
//    int pos2=neighbours[pos1][0];
//    cout<<"X["<< pos1 <<"," <<pos2<<"]="<< X(Spins,neighbours, pos1,pos2,1,+1) <<endl;
//    //double X(int* Spins,int** neighbours, int i,int j,double Temp,int sign);
//    string STOP;
//    cout<<"STOP!"<<endl;
//    cin>>STOP;
    

    double mean_E;
    double mean_E2;
    double mean_Cv;
    double mean_M;
    double mean_M2;
    double mean_M4;
    
    double mean_MA;
    double mean_MB;
    double E,Cv,prob,r;

    
    /// Loop over the temperatures ///
    for(int n_t=0; n_t<N_temps;n_t++){
//        break;
        
        double T=Temp[n_t];
        
        /// PRINTING CURRENT TEMPERATURE MOD 10 ///
        if(n_t%5==0){cout << "Temperature = "<< T<< endl;}
        
        ///VARIABLES FOR THE AVERAGES///
        mean_E=0;
        mean_E2=0;
        mean_M=0;
        
        mean_M2=0;
        mean_M4=0;
        
        mean_MA=0;
        mean_MB=0;

        /// Loop over equillibrium and measure sweeps ///
        for(int n_eq=0; n_eq<n_eqSweeps;n_eq++){
            sweep(Spins, neighbours,T,N_spins );
        }
        
        /// Measuring after equilibrating the system ///

            for(int n_eq=0; n_eq<n_measSweeps;n_eq++){
                
                sweep(Spins, neighbours,T,N_spins );
                
                /// Measuring after equilibrating the system ///
   
                E=get_Energy(Spins, neighbours,N_spins,T);
                
                mean_E  +=E/(n_measSweeps);
                mean_E2 +=pow(E,2)/(n_measSweeps);
                mean_Cv += get_Cv_correction(Spins, neighbours,N_spins,T)/(n_measSweeps);
                
                mean_M  +=get_magnetization(Spins, N_spins)/(n_measSweeps);
                mean_M2 +=pow(mean_M*n_measSweeps/N_spins,2)/(n_measSweeps);
                mean_M4 +=pow(mean_M*n_measSweeps/N_spins,4)/(n_measSweeps);
                
                //Bipartite lattice magnetization
                mean_MA +=get_magnetization_A(Spins, N_spins)/(n_measSweeps);
                mean_MB +=get_magnetization_B(Spins, N_spins)/(n_measSweeps);

        }
     // This function needs a correction due to the temperature dependence of the effectiva Hamiltonian
     Cv=(mean_E2-pow(mean_E,2))/(N_spins*pow(T,2))+ mean_Cv/(N_spins*pow(T,2));

     File_E_M << T << " " << mean_E/N_spins << " "<< Cv << " "<< mean_M/N_spins <<" "<< mean_MA/N_spins<< " "<< mean_MB/N_spins  <<  endl;
        
     File_Magnettization_powers<<  T << " " << pow(mean_M2,2) << " "<< mean_M4<<" "<< SQU <<endl;
        
    }
    File_E_M.close();
    File_Magnettization_powers.close();
    
//    for(int i=0;i<N_spins;i++){
//
//        cout<< Spins[i] << " ";
//
//
//        if((i+1)%SQU==0){
//            cout<<endl;
//        }
//
//    }
//    cout<<endl;
    
    return 0;
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////FUNCTIONS OF THE SYSTEM ////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

double X(int* Spins,int** neighbours, int i,int j,double Temp,int sign){
    double sum=-4*J*(sign*Spins[i])*Spins[j];
    
    double Neig_i=0;
    double Neig_j=0;
   
    for(int k=0;k<Z;k++){
        Neig_i+=Spins[neighbours[i][k]];
        Neig_j+=Spins[neighbours[j][k]];
        
        if(sign==-1){// If the sign is used then the ith sign must be added separatelly
            Neig_j-=2*Spins[neighbours[j][k]];
            
        }
    }
   
    sum+=2*J*(sign*Spins[i]*Neig_i+Spins[j]*Neig_j);
    sum/=(Temp*KB);// Taking away the - signs
    
    return sum;
    
    
}

double F(int* Spins,int** neighbours, int i,int j,double Temp,int sign){
    double sum=-4*J*(sign*Spins[i])*Spins[j];
    
    double Neig_i=0;
    double Neig_j=0;
//    cout<<"Sum before="<< sum<<endl;
//    cout<< "Spin[i="<<i<<"]="<<Spins[i]<<endl;
//    cout<< "Spin[j="<<j<<"]="<<Spins[j]<<endl;
    for(int k=0;k<Z;k++){
        Neig_i+=Spins[neighbours[i][k]];
        Neig_j+=Spins[neighbours[j][k]];
        
        if(sign==-1){// If the sign is used then the ith sign must be added separatelly
            Neig_j-=2*Spins[neighbours[j][k]];
            
        }
//        sum+=(Spins[i]*Spins[neighbours[j][k]]+ Spins[j]*Spins[neighbours[i][k]]);
    }
//    cout<<"Neigh_i "<< i<< "="<< Neig_i<<endl;
//    cout<<"Neigh_j "<< j<< "="<< Neig_j<<endl;
    sum+=2*J*(sign*Spins[i]*Neig_i+Spins[j]*Neig_j);
    sum/=(Temp*KB);// Taking away the - sign
    
    if(sum==0){
        //Using the Taylor expansion
        return 0.5;
    }
    return (exp(sum)-1-sum)/pow(sum,2);
    
    
}

double H(int* Spins,int** neighbours, int i,int j,double Temp,int sign){
    double sum=-4*J*(sign*Spins[i])*Spins[j];
    
    double Neig_i=0;
    double Neig_j=0;
    //    cout<<"Sum before="<< sum<<endl;
    //    cout<< "Spin[i="<<i<<"]="<<Spins[i]<<endl;
    //    cout<< "Spin[j="<<j<<"]="<<Spins[j]<<endl;
    for(int k=0;k<Z;k++){
        Neig_i+=Spins[neighbours[i][k]];
        Neig_j+=Spins[neighbours[j][k]];
        
        if(sign==-1){// If the sign is used then the i-th sign must be added separatelly
            Neig_j-=2*Spins[neighbours[j][k]];
            
        }
    }

    sum+=2*J*(sign*Spins[i]*Neig_i+Spins[j]*Neig_j);
    sum/=(Temp*KB);// taking away the - sign
    
    if(sum==0){
        //Using the Taylor expansion
        return 1.;
    }
    if((exp(sum)-1)/sum<0){
        cout<<"Negative Value! Something is wrong in H!!!"<<endl;
    }
    return (exp(sum)-1)/sum;
    
    
}


double get_Energy(int* Spins, int** neighbours,  int N, double Temp ){
    /// Energy variable to be computed ///
    double E=0;
    
    /// Loop over sites ///
    for(int i=0; i<N; i++){
        /// Loop over neighbours ///
        for(int j=0;j<Z; j++ ){
            
            E+=J*Spins[i]*Spins[neighbours[i][j]];
            
            if(Jperp !=0 && Spins[i]!=Spins[neighbours[i][j]]){
                E-=4*(pow(Jperp,2)/(Temp*KB))*(H(Spins, neighbours,  i, neighbours[i][j],Temp,+1));//*(1-Spins[i]*Spins[neighbours[i][j]])
            }// Get the energy only if the term contributes, i.e., if the interaction is AFM
        }
    }
    return E/2.;// It is still symmetric!!!

}

double get_Cv_correction(int* Spins, int** neighbours,  int N, double Temp ){
    
    /// Energy variable to be computed ///
    double Cv=0;
    
    /// Loop over sites ///
    for(int i=0; i<N; i++){
        /// Loop over neighbours ///
        for(int j=0;j<Z; j++ ){
            
            
            if(Jperp !=0 && Spins[i]!=Spins[neighbours[i][j]]){
                Cv+=2*pow(Jperp,2)*(1-Spins[i]*Spins[neighbours[i][j]])*exp(X(Spins, neighbours,  i, neighbours[i][j],Temp,+1));
                
            }// Get the energy only if the term contributes
        }
    }
    return Cv/2.;// divide by two due to overcounting
    
}

void sweep(int* Spins, int** neighbours, double T, int N_spins  ){
    int site;
    double DE;
    double prob,r;
    double H_pi_j;
    double H_ni_j;
    
    for(int N=0; N<N_spins; N++){
        /// site to be flipped ///
        site=(rand()%N_spins);

        /// Defining Energy of the flip ///
        DE=0;
        
        for(int j=0; j<Z;j++){
            DE -= 2*J*Spins[site]*Spins[neighbours[site][j]] ;
            
            if(Jperp!=0){
                H_pi_j=H(Spins, neighbours,  site, neighbours[site][j],T,+1);
                H_ni_j=H(Spins, neighbours,  site, neighbours[site][j],T,-1);
                
                DE -= 2*(pow(Jperp,2)/(T*KB))*(H_ni_j-H_pi_j);
                DE -= 2*(pow(Jperp,2)/(T*KB))*(Spins[site]*Spins[neighbours[site][j]])*(H_ni_j+H_pi_j);
                
            }
        }
        /// Defining probability
        prob=exp(-DE/(KB*T));
        r=random_double();
        
        /// Flipp the site if the conditions are fulfilled
        if( (r<prob) || (DE<0) ){
            Spins[site]*=-1;
        }
    }
}

double get_magnetization(int* Spins, int N){
    
    double M=0;
    for(int i=0; i<N; i++){
        M+=Spins[i];
    }
    return M;
    
}
double get_magnetization_A(int* Spins, int N){
    
    double M=0;
    for(int i=0; i<N; i++){
        
        if(int(sqrt(N))%2==1 && i%2==0){M+=Spins[i];}
        
        else{
            if(int(i/sqrt(N))%2==0 && i%2==0){
                
                M+=Spins[i];
            }
            if(int(i/sqrt(N))%2==1 && i%2==1){
                M+=Spins[i];
            }
        }
        
        
    }
    return M;
    
}
double get_magnetization_B(int* Spins, int N){
    
    double M=0;
    for(int i=0; i<N; i++){
        if(int(sqrt(N))%2==1 && i%2==1){M+=Spins[i];}
        
        else{
            if(int(i/sqrt(N))%2==0 && i%2==1){
                
                M+=Spins[i];
            }
            if(int(i/sqrt(N))%2==1 && i%2==0){
                M+=Spins[i];
            }
        }
    }
    return M;
    
}

double random_double(){
    
    const long max_rand = 1000000L;
    double x1 =0, x2 = 1., x;
    
    x = x1 + ( x2 - x1) * (random() % max_rand) / max_rand;
    
    
    return x;
    
}



