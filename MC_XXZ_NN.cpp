#include <iostream>
#include <cmath>
#include <fstream>
#include <time.h>
#include <string.h>
#include <iomanip>
#include <sstream>
using namespace std;


void sweep(int* Spins, int** neighbours,int** neighbours2,int** neighbours3, double T, int N_spins  );
double get_Energy(int* Spins, int** neighbours , int** neighbours2, int** neighbours3,  int N, double Temp );
double get_magnetization(int* Spins, int N);
double get_Cv_correction(int* Spins, int** neighbours , int** neighbours2, int** neighbours3,  int N, double Temp );
double get_magnetization_A(int* Spins, int N);
double get_magnetization_B(int* Spins, int N);

double random_double();

#define J 1.00 // Positive for AFM interactions
#define Jperp 0.20
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
    int** neighbours2 =new int*[N_spins];
    int** neighbours3 =new int*[N_spins];
    for(int i=0;i<N_spins;i++){
        neighbours[i]=new int[Z];
        neighbours2[i]=new int[Z];
        neighbours3[i]=new int[Z];
    }
    
    /// Fealing the neighbour matrix ///
    for(int i=0;i<L;i++){
        
        for(int j=0;j<H;j++){
            
            neighbours[i*L+j][0]=i*L+(j+1)%H; //Right neighbour
            neighbours[i*L+j][1]=((i-1+L)%L)*L+j; //Upwards neighbour
            neighbours[i*L+j][2]=i*L+(j-1+H)%H; //Left neighbour
            neighbours[i*L+j][3]=((i+1)%L)*L+j; //Downwards neighbour
            
            neighbours2[i*L+j][0]=((i-1+L)%L)*L+(j+1+H)%H; //Right neighbour
            neighbours2[i*L+j][1]=((i-1+L)%L)*L+(j-1+H)%H; //Upwards neighbour
            neighbours2[i*L+j][2]=((i+1+L)%L)*L+(j-1+H)%H; //Left neighbour
            neighbours2[i*L+j][3]=((i+1+L)%L)*L+(j+1+H)%H;; //Downwards neighbour
            
            neighbours3[i*L+j][0]=i*L+(j+2)%H; //Right neighbour
            neighbours3[i*L+j][1]=((i-2+L)%L)*L+j; //Upwards neighbour
            neighbours3[i*L+j][2]=i*L+(j-2+H)%H; //Left neighbour
            neighbours3[i*L+j][3]=((i+2)%L)*L+j; //Downwards neighbour
        }
    }
//    int positions=0;
//    cout<< "Initial position="<< positions<<endl;
//    for(int i=0; i<Z; i++){
//        cout<< "1st NN "<<neighbours[positions][i]<< endl;
//        cout<< "2nd NN "<<neighbours2[positions][i]<< endl;
//        cout<< "3rd NN "<<neighbours3[positions][i]<< endl;
//    }
//
//
    
    int N_temps=50;
    double T_max=4;
    double T_min=1.5;
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
    
    //Setting the name of the file to contain the information of Jperp//
    std::stringstream stream;
    stream << std::fixed << std::setprecision(2) << Jperp ;
    std::string s = stream.str();

    cout<<"J="<< J<<endl;
    cout<<"Jperp="<< Jperp<<endl;
    
    File_E_M.open("Energy_magnetization_J123_L"+ to_string(L)+"_Jperp"+s+".txt");
    File_Magnettization_powers.open("Magnetization.txt",ios::app);
  
 
    /// Starting with a fix spin configuration ///
    
    /// Random Spin configurations ////
    int* Spins=new int[N_spins];
    
    srand(time(NULL));

    for(int i=0;i<N_spins;i++){
        Spins[i]=2*(rand()%2)-1;
       
        cout<< Spins[i] << " ";
        if((i+1)%L==0){
            cout<<endl;
        }

    }

    double mean_E;
    double mean_E2;
    double mean_Cv;
    double mean_M;
    double mean_M2;
    double mean_M4;
    
    double mean_MA;
    double mean_MB;
    double E,Cv,prob,r;

    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////



    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    /// Loop over the temperatures ///
    for(int n_t=0; n_t<N_temps;n_t++){
//        break;
        
        double T=Temp[n_t];
        
        /// PRINTING CURRENT TEMPERATURE MOD 10 ///
        if(n_t%10==0){cout << "Temperature = "<< T<< endl;}
        
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
//            sweep(Spins, neighbours,T,N_spins );
            sweep(Spins, neighbours,neighbours2,neighbours3, T, N_spins  );

        }
        
        /// Measuring after equilibrating the system ///

            for(int n_eq=0; n_eq<n_measSweeps;n_eq++){
                
                sweep(Spins, neighbours,neighbours2,neighbours3, T, N_spins  );

                /// Measuring after equilibrating the system ///
   
                E=get_Energy(Spins, neighbours,neighbours2,neighbours3,N_spins,T);
                
                mean_E  +=E/(n_measSweeps);
                mean_E2 +=pow(E,2)/(n_measSweeps);
                mean_Cv += get_Cv_correction(Spins, neighbours,neighbours2,neighbours3,N_spins,T)/(n_measSweeps);
                
                mean_M  +=get_magnetization(Spins, N_spins)/(n_measSweeps);
                mean_M2 +=pow(mean_M*n_measSweeps/N_spins,2)/(n_measSweeps);
                mean_M4 +=pow(mean_M*n_measSweeps/N_spins,4)/(n_measSweeps);
                
                //Bipartite lattice magnetization
                mean_MA +=get_magnetization_A(Spins, N_spins)/(n_measSweeps);
                mean_MB +=get_magnetization_B(Spins, N_spins)/(n_measSweeps);

        }
     // This function needs a correction due to the temperature dependence of the effectiva Hamiltonian
     Cv=(mean_E2-pow(mean_E,2))/(N_spins*pow(T,2)) - mean_Cv/(N_spins*pow(T,2));

//     File_E_M << T << " " << mean_E/N_spins << " "<< Cv << " "<< mean_M/N_spins <<" "<< mean_MA/N_spins<< " "<< mean_MB/N_spins  <<  endl;
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
    cout<<"Jperp="<<s<<endl;
    return 0;
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////FUNCTIONS OF THE SYSTEM ////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


double get_Energy(int* Spins, int** neighbours , int** neighbours2, int** neighbours3,  int N, double Temp ){
    /// Energy variable to be computed ///
    double E=0;
    double J1,J2,J3;
    double beta=1./(KB*Temp);
    J1=J+pow(Jperp,2)*beta-(2./3.)*(2*Z-6)*J*pow(Jperp*beta,2);
    J2=4*J*pow(Jperp*beta,2)/3;
    J3=4*J*pow(Jperp*beta,2)/3;
    
    J1+=pow(Jperp,2)*beta-(4./3.)*(2*Z-6)*J*pow(Jperp*beta,2);
    J2+=8*J*pow(Jperp*beta,2)/3;
    J3+=8*J*pow(Jperp*beta,2)/3;
    
    /// Loop over sites ///
    for(int i=0; i<N; i++){
        /// Loop over neighbours ///
        for(int j=0;j<Z; j++ ){
            
            E+=J1*Spins[i]*Spins[neighbours[i][j]];
            E+=J2*Spins[i]*Spins[neighbours2[i][j]];
            E+=J3*Spins[i]*Spins[neighbours3[i][j]];
        
        }
    }
    
    return E/2.-pow(Jperp,2)*beta*N*Z;// It is still symmetric!!!

}

double get_Cv_correction(int* Spins, int** neighbours , int** neighbours2, int** neighbours3,  int N, double Temp ){
    
    /// Energy variable to be computed ///
    double Cv=0;
    double J1,J2,J3;
    double beta=1./(KB*Temp);
    
    J1=2*(pow(Jperp,2)-(4./3.)*J*(2*Z-6)*beta*pow(Jperp,2) );
    J2=2*8*J*beta*pow(Jperp,2)/3;
    J3=2*8*J*beta*pow(Jperp,2)/3;
    
    J1+= -(4./3.)*(2*Z-6)*J*beta*pow(Jperp,2);
    J2+=8*J*beta*pow(Jperp,2)/3;
    J3+=8*J*beta*pow(Jperp,2)/3;
    
    /// Loop over sites ///
    for(int i=0; i<N; i++){
        /// Loop over neighbours ///
        for(int j=0;j<Z; j++ ){
            
            Cv+=J1*Spins[i]*Spins[neighbours[i][j]];
            Cv+=J2*Spins[i]*Spins[neighbours2[i][j]];
            Cv+=J3*Spins[i]*Spins[neighbours3[i][j]];
            
        }
    }
    
    return Cv/2.- pow(Jperp,2)*N*Z;// It is still symmetric!!!
    
}

void sweep(int* Spins, int** neighbours,int** neighbours2,int** neighbours3, double T, int N_spins  ){
    int site;
    double DE;
    double prob,r;
    int posj,posj2,posj3;

    
    double J1,J2,J3;
    double beta=1./(KB*T);
    
    J1=J+pow(Jperp,2)*beta-(2./3.)*J*pow(Jperp*beta,2);
    J2=4*J*pow(Jperp*beta,2)/3;
    J3=4*J*pow(Jperp*beta,2)/3;
    //Correction
    J1+=pow(Jperp,2)*beta-(4./3.)*J*pow(Jperp*beta,2);
    J2+=8*J*pow(Jperp*beta,2)/3;
    J3+=8*J*pow(Jperp*beta,2)/3;
    
    
    for(int N=0; N<N_spins; N++){
        /// site to be flipped ///
        site=(rand()%N_spins);

        DE=0;
        
        for(int j=0; j<Z;j++){
            posj=neighbours[site][j];//Position of the jth neighbour
            posj2=neighbours2[site][j];//Position of the jth neighbour
            posj3=neighbours3[site][j];//Position of the jth neighbour
            
            DE -= 2*J1*Spins[site]*Spins[posj];
            DE -= 2*J2*Spins[site]*Spins[posj];
            DE -= 2*J3*Spins[site]*Spins[posj];
            //        cout<<neighbours[site1][j]<<endl;
            //        cout<<"Only J DE= "<<DE<<endl;
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

