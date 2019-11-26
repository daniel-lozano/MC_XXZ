#include<iostream>
#include<iomanip>
#include<assert.h>
#include<ctime>
#include <time.h>
#include<random>
#include<string>
#include<sstream>
#include<functional>
#include<vector>
#include<algorithm>
#include<cmath> // complex math
#include<complex> //complex numbers
#include <fstream> // Files
#include <string>

using namespace std;
std::random_device rd;  //Will be used to obtain a seed for the random number engine
std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
std::uniform_real_distribution<> dis(0., 1.0);
//std::uniform_real_distribution<> dis_2(0.8,1.2);
std::normal_distribution<> dis_2(1.0 ,0.2);
//std::normal_distribution<> dis_3(0.1,0.2);


void sweep(int* Spins,int* Bonds,double* Bonds_vals, int** neighbours, int** neighbours2, double T, int N_spins  );
double get_Energy(int* Spins, int* Bonds,double* Bonds_vals, int** neighbours, int** neighbours2,  int N );
double get_magnetization(int* Spins, int* Bonds, int N);
void print_matrix(double* Matrix,int N,int L );
double plaquette_density(int* Spins,int** neighbours,int* Spins_bonds , int N_spins, int L, int H,double T);
double jackknife(double* array,int size);

#define J 1// Positive for antiferro
#define J2 1
#define KB 1.
#define Z 6

#define H_field 0.0 // when it's positive it should promote having spins pointing up, i.e. -H_field*\sum spins


int main(int argc, char** argv){
    
    /// Constants of the model ///
   
    
    /// Lattice constants ///
    int SQU=atoi(argv[1]);
    int L,H;
    L=SQU;
    H=SQU;
    int N_spins=L*H;
    double Tc=2.0/log(1+sqrt(2))*J; //
    
    int n_eqSweeps=5000;
    int n_measSweeps=10000;
    int trials=1;
    cout << "The size of the lattice is "<< SQU << endl;
    cout << "The external field is H_field="<< H_field << endl;
    
    
    /// Changin seed to a time dependend seed ///
    srand(time(NULL));
    
    cout << "N_spins=" << N_spins << endl;
    cout << "N_bonds=" << 2*L*L << endl;

    /// Creating matrix of bonds ///

    int M=2*L;

    
    /// Creating array of neighbours ///
    int** neighbours =new int*[N_spins];
    int** neighbours2 =new int*[N_spins];
    for(int i=0;i<N_spins;i++){
        neighbours[i]=new int[Z];
        neighbours2[i]=new int[Z];
    }
    /// Fealing the neighbour matrix ///
    for(int i=0;i<H;i++){
        for(int j=0;j<L;j++){
            
            ///First NN
            neighbours[i*L+j][0]=i*L+(j+1)%L; //Right neighbour
                
            neighbours[i*L+j][1]=((i-1+H)%H)*L+j; //Upwards neighbour
            
            neighbours[i*L+j][2]=((i-1+H)%H)*L+(j-1+L)%L; //Upwards diagonal
            
            neighbours[i*L+j][3]=i*L+(j-1+L)%L; //Left neighbour
                
            neighbours[i*L+j][4]=((i+1)%H)*L+j; //Downwards neighbour
            
            neighbours[i*L+j][5]=((i+1)%H)*L+(j+1)%L; //Downwards diagonal
            
            ///Second NN
           
            neighbours2[i*L+j][0]=((i+1)%H)*L+(j+1)%L; //Right up
            
            neighbours2[i*L+j][1]=((i+2)%H)*L+(j-1+L)%L; //Up
            
            neighbours2[i*L+j][2]=((i+1)%H)*L+(j-2+L)%L; //Left up
            
            neighbours2[i*L+j][3]=((i-1+H)%H)*L+(j-1+L)%L; //Left down
        
            neighbours2[i*L+j][4]=((i-2+H)%H)*L+(j+1)%L; //Downwards
            
            neighbours2[i*L+j][5]=((i-1+H)%H)*L+(j+2)%L; //Right down
//            cout<<"Position i="<<i*L+j<<endl;
//            cout<<((i+1)%H)*L+(j+1)%L<<endl;
//            cout<<((i+2)%H)*L+(j-1+L)%L<<endl;
//            cout<<((i+1)%H)*L+(j-2+L)%L<<endl;
//            cout<<((i-1+H)%H)*L+(j-1+L)%L<<endl;
//            cout<<((i-2+H)%H)*L+(j+1)%L<<endl;
//            cout<<((i-1+H)%H)*L+(j+2)%L<<endl;
        }
    }
    
 
    
    int N_temps=50;
    double T_max=5.0;//5
    double T_min=0.1;//2
    double* Temp=new double[N_temps];
    
    for(int i=0; i<N_temps;i++){
        Temp[i]=T_max-i*(T_max-T_min)/N_temps;
    }
    /// Files to be use ///
    ofstream File_E_M, File_configurations,File_labels, File_temps;
    
    File_E_M.open("Energy_magnetization_L"+ to_string(L)+".txt");
    File_configurations.open("spin_configurations_mattis_L"+ to_string(L)+".txt");

    
    cout << "Number of trials " << trials << endl;
    
    
    for(int t=0; t<trials; t++){
        /// Starting with a fix spin configuration ///
        cout<< "Fealing spin configurations "<< t << endl;
        
        /// Random Spin configurations ////
        int* Spins=new int[N_spins];
        int* Spins_bonds=new int[2*N_spins];

        /// Random Spin bonds configurations ////
        int* Random_variable=new int[N_spins];
        /// Random values of the bonds ///
        double* Random_values=new double[N_spins];
        
        for(int i=0;i<N_spins;i++){
            Spins[i]=2*(rand()%2)-1;
            Random_variable[i]=1.;
            Random_values[i]=1.;
        }
        
        for(int i=0;i<2*N_spins;i++){
            Spins_bonds[i]=0;
            
        }
                    
        
      
//        normalize_bonds(Random_values, neighbours, N_spins, L);
//        cout<< "Random_values after"<<endl;
//        print_matrix(Random_values,N_spins,L);
//        cout<<endl;
        ofstream File_bonds;

        double mean_E=0;
        double mean_E2=0;
        double mean_M=0;
        double mean_M2=0;
        double mean_M3=0;
        double mean_M4=0;
        double M=0;
        double mean_droplet_density=0;

        double E,Cv,prob,r, X1, X3;
        double error_E=0;
        double* Array_Energy=new double[n_measSweeps];
        
       
        /// Loop over the temperatures ///
        for(int n_t=0; n_t<N_temps;n_t++){
            //        break;
            
            double T=Temp[n_t];
            
            /// PRINTING CURRENT TEMPERATURE MOD 10 ///
            if(n_t%10==0){
                
                
                cout << "Temperature = "<< T<< endl;

            }
            
            
            
            
            ///VARIABLES FOR THE AVERAGES///
            mean_E=0;
            mean_E2=0;
            mean_M=0;
            mean_M2=0;
            mean_M3=0;
            mean_M4=0;
            mean_droplet_density=0;
            
            /// Loop over equillibrium and measure sweeps ///
            for(int n_eq=0; n_eq<n_eqSweeps;n_eq++){
                
                sweep(Spins,Random_variable,Random_values, neighbours,neighbours2,T,N_spins );
                
            }
                /// Measuring after equilibrating the system ///
            for(int n_eq=0; n_eq<n_measSweeps;n_eq++){
                
                sweep(Spins,Random_variable,Random_values, neighbours,neighbours2,T,N_spins );
                
                E=get_Energy(Spins, Random_variable, Random_values, neighbours,neighbours2,N_spins);
                
                Array_Energy[n_eq]=E/N_spins;
                mean_E+=E/(n_measSweeps);
                
                M=get_magnetization(Spins, Random_variable,N_spins)/N_spins;//
                mean_M+=M/((n_measSweeps));
                mean_M2+=pow(M,2)/((n_measSweeps));
                mean_M3+=pow(M,3)/((n_measSweeps));
                mean_M4+=pow(M,4)/((n_measSweeps));
                mean_E2+=pow(E,2)/(n_measSweeps);
//                mean_droplet_density+=plaquette_density(Spins,neighbours,Spins_bonds , N_spins,L,  H,T)/(n_measSweeps);
                
                /// Saving configurations ///
                    
            }
            
            Cv=( mean_E2-pow(mean_E,2))/(N_spins*pow(T,2));
            X1=( mean_M2-pow(mean_M,2))/(pow(T,1));
            X3=( mean_M4-4*(mean_M*mean_M3)- 3*pow(mean_M2,2)+12*pow(mean_M,2)*mean_M2-6*pow(mean_M,4))/(pow(T,3));
            error_E=jackknife(Array_Energy,n_measSweeps);
            
            File_E_M << T << " " << mean_E/N_spins << " "<< mean_M << " "<< Cv<< " "<< mean_droplet_density<< " "<< mean_M2 << " "<< mean_M4 << " "<< X1 << " "<< X3<< " "<< error_E<< endl;
           
        }
        if(1){
            int config=1;//2*(rand()%2)-1;
            /// Printing fewer configurations ///
            for(int i=0; i< N_spins;i++){
                File_configurations<< Spins[i]*config << " ";
            }
            File_configurations<< endl;
           
        }
        File_E_M.close();
        File_configurations.close();
        
        
    }
    
    

    return 0;
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////FUNCTIONS OF THE SYSTEM ////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


double get_Energy(int* Spins, int* Bonds,double* Bonds_vals, int** neighbours, int** neighbours2,  int N ){
    /// Energy variable to be computed ///
    double E=0;
    double M=0;
    
    /// Loop over sites ///
    for(int i=0; i<N; i++){
        M+=Spins[i];
        
        /// Loop over neighbours ///
        for(int j=0;j<Z; j++ ){
            if(J!=0){
                E+= J*Spins[i]*Spins[neighbours[i][j]];
            }
            if(J2!=0){
                E+= J2*Spins[i]*Spins[neighbours2[i][j]];
            }
        }
    }

    return E/2.- H_field*M;
    
}
void sweep(int* Spins,int* Bonds,double* Bonds_vals, int** neighbours,int** neighbours2, double T, int N_spins  ){
    int site;
    double DE;
    double prob,r;
    double Jij=1;
    double Jij2=1;
    
    
    for(int N=0; N<N_spins; N++){
        /// site to be flipped ///
        site=(rand()%N_spins);
        
        DE=0;
        DE+= 2*H_field*Spins[site];//Energy change by the field
        
        for(int j=0; j<Z;j++){
           
            
            DE+= -2*J*Spins[site]*Spins[neighbours[site][j]];
            
            DE+= -2*J2*Spins[site]*Spins[neighbours2[site][j]];
            
        }
        /// Defining probability
        prob=exp(-DE/(KB*T));
        r=dis(gen);
        
        /// Flipp the site if the conditions are fulfilled
        if( (r<prob) || (DE<0) ){
            Spins[site]*=-1;
        }
    }
}

double get_magnetization(int* Spins, int* Bonds, int N){
    
    double M=0;
    for(int i=0; i<N; i++){
        M+=Spins[i];
    }
    return M;
    
}


void print_matrix(double* Matrix,int N,int L ){
    for(int i=0;i<N;i++){
        cout<< Matrix[i]<<" ";
        if((i+1)%L==0){
            cout<<endl;
        }
    }
    
}

double plaquette_density(int* Spins,int** neighbours,int* Spins_bonds , int N_spins, int L, int H,double T){
    
  
    int pos=0;
    double rand1;
    double rand2;
    for(int i=0;i<H;i++){
        
         for(int j=0;j<L;j++){
             
             pos=i*L+j;
             //Initialice positions in bonds
             Spins_bonds[2*i*L+j]=0;
             Spins_bonds[2*i*L+L+j]=0;
             rand1=dis(gen);
             rand2=dis(gen);
             //Generating bonds for the droplets
             
             if(Spins[pos]*Spins[neighbours[pos][0]]==1 && Spins[pos]==-1 && rand1<(1-exp(-2*J/T)) ){//Vertical bonds
                 Spins_bonds[2*i*L+j]=1;
                 
             }
             if(Spins[pos]*Spins[neighbours[pos][3]]==1 && Spins[pos]==-1 && rand2<(1-exp(-2*J/T)) ){//Vertical bonds
                Spins_bonds[2*i*L+j]=1;
                             
            }
             
         }
    }
    
    double val=0;
    //Counting dimmers bonds
    for(int i=0;i<2*H;i++){
        for(int j=0;j<L;j++){
            
            //checking for horizontal bonds
            if( i%2==0 && Spins_bonds[i*L+j]==1 &&
               Spins_bonds[i*L+(j+1)%L]==0 && // right
               Spins_bonds[i*L+(j-1+L)%L]==0 && // left
               
               Spins_bonds[((i-1+2*H)%(2*H))*L+(j)%L]==0 && // up
               Spins_bonds[((i+1)%(2*H))*L+(j)%L]==0 && // down
               
               Spins_bonds[((i-1+2*H)%(2*H))*L+(j+1)%L]==0 && // up right
               Spins_bonds[((i+1)%(2*H))*L+(j+1)%L]==0 ){ // down right
                val+=1;
            }
            if( i%2==1 && Spins_bonds[i*L+j]==1 &&
               
               Spins_bonds[((i-2+2*H)%(2*H))*L+j]==0 && // up
               Spins_bonds[((i+2)%(2*H))*L+j]==0 && // down
               
               Spins_bonds[((i-1+2*H)%2*H)*L+j]==0 && // up right
               Spins_bonds[((i+1)%2*H)*L+j]==0 && // down right
               
               Spins_bonds[((i-1+2*H)%2*H)*L+(j-1)%L]==0 && // up left
               Spins_bonds[((i+1)%2*H)*L+(j-1)%L]==0 ){ // down left
                val+=1;
            }
               
               
            }
        }
    
//    cout<<"Val="<<val<<endl;
//
//    for(int i=0;i<2*N_spins;i++){
//
//        cout<< Spins_bonds[i]<<" ";
//        if(i%L==L-1){
//           cout<<endl;
//       }
//    }
    return val/(N_spins);
}

double jackknife(double* array,int size){
    double av_mean=0;
    double mean_j=0;
    double error=0;
    int sub_bin_size=size/10;
    int counter=0;
    
    // Compute the mean value of the array
    for(int i=0; i<size; i++){
        av_mean+=array[i]/size;
    }
    
    
    for(int i=0;i<10;i++){
        counter=0;
        mean_j=0;
        
        for(int j=0;j<size/10;j++){
            
            if(j+i*sub_bin_size<size){
                counter++;
                mean_j+=array[j+i*sub_bin_size];
            }
        }
        error+=(sub_bin_size)*pow(av_mean-(mean_j/counter),2)/size;
    }
    
//    for(int i=0; i<size; i++){
//        mean_j=0;
//        //Compute the mean value of the array with out the i-th component
//        for(int j=0; j<size; j++){
//
//            if(j!=i){
//                mean_j+=array[j]/(size-1);
//            }
//        }
//        //Compute the error of the mean as (N-1)\sum_i (mean_i-mean)**2 /N
//        error+=(size-1)*pow(av_mean-mean_j,2)/size;
//    }
    
    return sqrt(error);
    
}


