#include <iostream>
#include <cmath>
#include <fstream>
#include <time.h>
#include <string.h>
using namespace std;


void sweep(int* Spins, int** neighbours,double T, int N_spins  );
double get_Energy(int* Spins, int** neighbours, int N  );
double get_magnetization(int* Spins, int N);
double random_double();
double F(int* Spins,int** neighbours, int i,int j,double Temp);



#define J 1.
#define Jperp 0.001
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
    
    int n_eqSweeps=1000;
    int n_measSweeps=5000;
  
    cout << "The size of the lattice is "<< SQU << endl;
    
    
    /// Changin seed to a time dependend seed ///
    srand(time(NULL));
    
    cout << "N_spins=" << N_spins << endl;
    cout << "N_bonds=" << 2*L*L << endl;
    
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
            
            neighbours[i*L+j][1]=((i-1)%L)*L+j; //Upwards neighbour
            if(i==0){neighbours[i*L+j][1]=(L-1)*L+j; }
            
            neighbours[i*L+j][2]=i*L+(j-1)%H; //Left neighbour
            if(j==0){neighbours[i*L+j][2]=i*L+H-1; }
            
            neighbours[i*L+j][3]=((i+1)%L)*L+j; //Downwards neighbour
        }
    }
    
 
    
    int N_temps=50;
    double T_max=5;
    double T_min=1;
    double* Temp=new double[N_temps];
    
    for(int i=0; i<N_temps;i++){
        Temp[i]=T_max-i*(T_max-T_min)/N_temps;
    }
    /// Files to be use ///
    ofstream File_E_M;
    
    File_E_M.open("Energy_magnetization_L"+ to_string(L)+".txt");
  
 
    /// Starting with a fix spin configuration ///
//    cout<< "Fealing spin configurations "<< t << endl;
    
    /// Random Spin configurations ////
    int* Spins=new int[N_spins];
    
    
    for(int i=0;i<N_spins;i++){
        Spins[i]=2*(rand()%2)-1;
       
        cout<< Spins[i] << " ";
        if((i+1)%10==0){
            cout<<endl;
        }

    }
    cout<<endl;
    int pos1=5;
    int pos2=neighbours[pos1][0];
    cout<<"F["<< pos1 <<"," <<pos2<<"]" <<endl;
    cout<<F(Spins,neighbours, pos1,pos2,1)<<endl;
    
    double mean_E=0;
    double mean_E2=0;
    double mean_M=0;
    double E,Cv,prob,r;

    
    /// Loop over the temperatures ///
    for(int n_t=0; n_t<N_temps;n_t++){
        break;
        
        double T=Temp[n_t];
        
        /// PRINTING CURRENT TEMPERATURE MOD 10 ///
        if(n_t%10==0){cout << "Temperature = "<< T<< endl;}
        
        ///VARIABLES FOR THE AVERAGES///
        mean_E=0;
        mean_E2=0;
        mean_M=0;
        
        /// Loop over equillibrium and measure sweeps ///
        for(int n_eq=0; n_eq<n_eqSweeps+n_measSweeps;n_eq++){

            sweep(Spins, neighbours,T,N_spins );
            
            /// Measuring after equilibrating the system ///
            
            if(n_eq>n_eqSweeps){

                E=get_Energy(Spins, neighbours,N_spins);

                mean_E+=E/(n_measSweeps);
                mean_M+=get_magnetization(Spins, N_spins)/(N_spins*(n_measSweeps));
                mean_E2+=pow(E,2)/(n_measSweeps);
                
                /// Saving configurations ///
               

            }

        }

     Cv=( mean_E2-pow(mean_E,2))/(N_spins*pow(T,2));

     File_E_M << T << " " << mean_E/N_spins << " "<< mean_M << " "<< Cv << endl;

    }
    File_E_M.close();

    return 0;
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////FUNCTIONS OF THE SYSTEM ////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

double F(int* Spins,int** neighbours, int i,int j,double Temp){
    double sum=-2*Spins[i]*Spins[j];
    double Neig_i=0;
    double Neig_j=0;
    cout<<"Sum before="<< sum<<endl;
    cout<< "Spin[i="<<i<<"]="<<Spins[i]<<endl;
    cout<< "Spin[j="<<j<<"]="<<Spins[j]<<endl;
    for(int k=0;k<Z;k++){
        Neig_i+=Spins[neighbours[i][k]];
        Neig_j+=Spins[neighbours[j][k]];
//        sum+=(Spins[i]*Spins[neighbours[j][k]]+ Spins[j]*Spins[neighbours[i][k]]);
    }
    cout<<"Neigh_i "<< i<< "="<< Neig_i<<endl;
    cout<<"Neigh_j "<< j<< "="<< Neig_j<<endl;
    sum+=(Spins[i]*Neig_j+Spins[j]*Neig_i);
    sum*=-2*J/Temp;
    
    if(sum==0){
        cout<<"sum=0!!!"<<endl;
        return 0.5;
    }
    return (exp(sum)-1-sum)/pow(sum,2);
    
    
}

double get_Energy(int* Spins, int** neighbours,  int N ){
    /// Energy variable to be computed ///
    double E=0;
    
    /// Loop over sites ///
    for(int i=0; i<N; i++){
        /// Loop over neighbours ///
        for(int j=0;j<Z; j++ ){
            E-=J*Spins[i]*Spins[neighbours[i][j]];//-2*(pow(Jperp,2)/T)*(1-Spins[i]*Spins[neighbours[i][j]])*F(Spins, neighbours,  i, j,Temp)
        }
    }
    return E/2.;// It is still symmetric!!!

}
void sweep(int* Spins, int** neighbours, double T, int N_spins  ){
    int site;
    double DE;
    double prob,r;
    
    
    for(int N=0; N<N_spins; N++){
        /// site to be flipped ///
        site=(rand()%N_spins);

        /// Defining Energy of the flip ///
        DE=0;
        
        for(int j=0; j<Z;j++){
            DE+=2*J*Spins[site]*Spins[neighbours[site][j]];
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

double random_double(){
    
    const long max_rand = 1000000L;
    double x1 =0, x2 = 1., x;
    
    x = x1 + ( x2 - x1) * (random() % max_rand) / max_rand;
    
    
    return x;
    
}



