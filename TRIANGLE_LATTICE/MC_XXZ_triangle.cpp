#include <iostream>
#include <cmath>
#include <fstream>
#include <time.h>
#include <string.h>
#include <iomanip>
#include <sstream>
using namespace std;


void sweep(int* Spins, int** neighbours,double T, int N_spins  );
double get_Energy(int* Spins, int** neighbours, int N , double Temp );
double get_Cv_correction(int* Spins, int** neighbours, int N , double Temp );
double get_magnetization(int* Spins, int N);

double get_magnetization_A(int* Spins, int N);
double get_magnetization_B(int* Spins, int N);

double jackknife(double* array,int size);



double random_double();
double X(int* Spins,int** neighbours, int i,int j,double Temp,int sign);
double F(int* Spins,int** neighbours, int i,int j,double Temp, int sign);
double H_func(int* Spins,int** neighbours, int i,int j,double Temp,int sign);
double H_DE(int* Spins,int** neighbours, int i,int j,double Temp,int sign, int pos);

#define J 1.00 // Positive for AFM interactions
double Jperp =0.05;
#define KB 1.
#define Z 6

int main(int argc, char** argv){
    
    /// Constants of the model ///
    
    /// Lattice constants ///
    int SQU=atoi(argv[1]);
    Jperp=atof(argv[2]);
    int L,H;
    L=SQU;
    H=SQU;//
    int N_spins=L*H;
    
//    double Tc=2.0/log(1+sqrt(2))*J; //
    
    int n_eqSweeps=5000;//1000000;
    int n_measSweeps=10000;//000000;
  
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
     for(int i=0;i<H;i++){
           for(int j=0;j<L;j++){
               neighbours[i*L+j][0]=i*L+(j+1)%L; //Right neighbour
                   
               neighbours[i*L+j][1]=((i-1+H)%H)*L+j; //Upwards neighbour
               
               neighbours[i*L+j][2]=((i-1+H)%H)*L+(j-1+L)%L; //Upwards diagonal
               
               neighbours[i*L+j][3]=i*L+(j-1+L)%L; //Left neighbour
                   
               neighbours[i*L+j][4]=((i+1)%H)*L+j; //Downwards neighbour
               
               neighbours[i*L+j][5]=((i+1)%H)*L+(j+1)%L; //Downwards diagonal
               
           }
       }
    
 
//    cout<<"Problem(?)"<<endl;

    int N_temps=50;
    double T_max=2.5;
    double T_min=0.8;
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
    
    File_E_M.open("Energy_magnetization_L"+ to_string(L)+"_Jperp"+s+".txt");
    File_Magnettization_powers.open("Magnetization_Jperp"+s+".txt",ios::app);
  
 
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

    double mean_E;
    double mean_E2;
    double mean_Cv;
    double mean_M;
    double mean_M2;
    double mean_M4;
    
    double mean_MA;
    double mean_MB;
    double E,Cv,prob,r;
    double error_M2=0;
    double error_M4=0;
    double* Array_M2=new double[n_measSweeps];
    double* Array_M4=new double[n_measSweeps];


    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//    double DE_check=0;
//    int site1=0;
//    int posj,posk;
//    double T=0.5;
//    double E_check=0;
//    double H_ij_pi,H_ij_ni;
//    double H_jk_pi,H_jk_ni;
//    
//   
//    
//    for(int j=0; j<Z;j++){
//        posj=neighbours[site1][j];//Position of the jth neighbour
//        
//        DE_check -= 2*J*Spins[site1]*Spins[posj];
////        cout<<neighbours[site1][j]<<endl;
////        cout<<"Only J DE= "<<DE<<endl;
//        if(Jperp!=0){
//            
////            H_ij_pi=H_func(Spins, neighbours,  site1, posj,T,+1);
////            H_ij_ni=H_func(Spins, neighbours,  site1, posj,T,-1);
////
//            H_ij_pi=H_DE(Spins, neighbours,  site1, posj,T,+1,site1);
//            H_ij_ni=H_DE(Spins, neighbours,  site1, posj,T,-1,site1);
//            
//            
//            DE_check -= 2.*(pow(Jperp,2)/(T*KB))*(H_ij_ni-H_ij_pi);
//            
//            DE_check -= 2.*(pow(Jperp,2)/(T*KB))*(Spins[site1]*Spins[posj])*(H_ij_ni+H_ij_pi);
//            
//            for(int k=0; k<Z; k++){
//                posk=neighbours[posj][k];
//                
//                if(posk!=site1 && Spins[posj]!=Spins[posk] ){
//                    
//                    H_jk_pi=H_DE(Spins,neighbours, posj,posk,T,+1, site1);
//                    H_jk_ni=H_DE(Spins,neighbours, posj,posk,T,-1, site1);
//                    DE_check -=2.*(pow(Jperp,2)/(T*KB))*(1-Spins[posj]*Spins[posk])*( H_jk_ni-H_jk_pi );
//                }
//            }
////            cout<<"Full, DE= "<<DE<<endl;
//        }
//    }
//    cout<<"with function DE="<<DE_check<<endl;
//    E_check=-get_Energy(Spins, neighbours,N_spins,T);
//    Spins[site1]*=-1;
//    E_check+=get_Energy(Spins, neighbours,N_spins,T);
//    cout<< "with energy function DE="<<E_check<<endl;
//    cout<<"Overal difference is= "<< E_check-DE_check<<endl;
//    string STOP;
//    cout<<"STOP"<<endl;
//    cin>>STOP;
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    /// Loop over the temperatures ///
    for(int n_t=0; n_t<N_temps;n_t++){
//        break;
        
        double T=Temp[n_t];
        double M=0;
        /// PRINTING CURRENT TEMPERATURE MOD 10 ///
        if(n_t%10==0){cout << "Temperature = "<< T<< endl;}
        
        ///VARIABLES FOR THE AVERAGES///
        mean_E=0;
        mean_E2=0;
        mean_M=0;
        
        mean_M2=0;
        mean_M4=0;
        
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
                M=get_magnetization(Spins, N_spins);
                mean_M  +=M/(n_measSweeps);
                mean_M2 +=pow(M/N_spins,2)/((n_measSweeps));
                mean_M4 +=pow(M/N_spins,4)/((n_measSweeps));
                Array_M2[n_eq]=pow(M/N_spins,2);
                Array_M4[n_eq]=pow(M/N_spins,4);
                
                
        }
     // This function needs a correction due to the temperature dependence of the effectiva Hamiltonian
     Cv=(mean_E2-pow(mean_E,2))/(N_spins*pow(T,2))+ mean_Cv/(N_spins*pow(T,2));
     error_M2=jackknife(Array_M2,n_measSweeps);
     error_M4=jackknife(Array_M4,n_measSweeps);
        

     File_E_M << T << " " << mean_E/N_spins << " "<< Cv << " "<< mean_M/N_spins <<  endl;
        
     File_Magnettization_powers<<  T << " " << pow(mean_M2,2) << " "<< mean_M4<< " "<< error_M2<< " " << error_M4 <<" "<< SQU <<endl;
        
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

double X(int* Spins,int** neighbours, int i,int j,double Temp,int sign){
    double sum=-4.*J*(sign*Spins[i])*Spins[j];
    
    double Neig_i=0;
    double Neig_j=0;
   
    for(int k=0;k<Z;k++){
        Neig_i+=Spins[neighbours[i][k]];
        Neig_j+=Spins[neighbours[j][k]];
    }
    
    if(sign==-1){// If the sign is used then the ith sign must be added separatelly
        Neig_j-=2*Spins[i];
    }
   
    sum+=2*J*(sign*Spins[i]*Neig_i+Spins[j]*Neig_j);
    sum/=(Temp*KB);// Taking away the - signs
    
    return sum;
    
    
}

double F(int* Spins,int** neighbours, int i,int j,double Temp,int sign){
    double sum=-4.*J*(sign*Spins[i])*Spins[j];
    
    double Neig_i=0;
    double Neig_j=0;
//    cout<<"Sum before="<< sum<<endl;
//    cout<< "Spin[i="<<i<<"]="<<Spins[i]<<endl;
//    cout<< "Spin[j="<<j<<"]="<<Spins[j]<<endl;
    for(int k=0;k<Z;k++){
        Neig_i+=Spins[neighbours[i][k]];
        Neig_j+=Spins[neighbours[j][k]];
        
    }
    if(sign==-1){// If the sign is used then the ith sign must be added separatelly
        Neig_j-=2*Spins[i];
        
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

double H_func(int* Spins,int** neighbours, int i,int j,double Temp,int sign){
    
    double sum=-4.*J*(sign*Spins[i])*Spins[j];
    double Neig_i=0;
    double Neig_j=0;
    //    cout<<"Sum before="<< sum<<endl;
    //    cout<< "Spin[i="<<i<<"]="<<Spins[i]<<endl;
    //    cout<< "Spin[j="<<j<<"]="<<Spins[j]<<endl;
    for(int k=0;k<Z;k++){
        Neig_i+=Spins[neighbours[i][k]];
        Neig_j+=Spins[neighbours[j][k]];
    }
    if(sign==-1){// If the sign is used then the i-th sign must be added separatelly
        Neig_j-=2*Spins[i];
        
    }

    sum+=2.*J*(sign*Spins[i]*Neig_i+Spins[j]*Neig_j);
    sum/=(Temp*KB);// taking away the - sign
    
    if(sum==0){
        //Using the Taylor expansion
        return 1.0;
    }
    if((exp(sum)-1.)/sum<0){
        cout<<"Negative Value! Something is wrong in H!!!"<<endl;
    }
    return (exp(sum)-1)/sum;
    
    
}

double H_DE(int* Spins,int** neighbours, int i,int j,double Temp,int sign, int pos){
    
    double sum=-4.*J*Spins[i]*Spins[j];
    
    if( (pos==i || pos ==j) && sign==-1 ){
        sum*=sign;
    }
    
    double Neig_i=0;
    double Neig_j=0;
    
    int pik,pjk;
    
    //    cout<<"Sum before="<< sum<<endl;
    //    cout<< "Spin[i="<<i<<"]="<<Spins[i]<<endl;
    //    cout<< "Spin[j="<<j<<"]="<<Spins[j]<<endl;
    for(int k=0;k<Z;k++){
        pik=neighbours[i][k];
        pjk=neighbours[j][k];
        
        Neig_i+=Spins[pik];
        Neig_j+=Spins[pjk];
        
        //Additional step which flips the spin in the site pos
        if(pos==pik && sign==-1 ){// Substract if the neighbour was suppose to be flipped
            Neig_i-=2*Spins[pik];
        }
        if(pos==pjk && sign==-1 ){
            Neig_j-=2*Spins[pjk];
        }
    }
    if(pos==i && sign==-1 ){
        sum+=2.*J*(sign*Spins[i]*Neig_i);
    }
    else{
        sum+=2.*J*(Spins[i]*Neig_i);
        
    }
    if(pos==j && sign==-1){
        sum+=2.*J*(sign*Spins[j]*Neig_j);
    }
    else{
        sum+=2.*J*(Spins[j]*Neig_j);
    }
    
    sum/=(Temp*KB);// taking away the - sign
    
    if(sum==0){
        //Using the Taylor expansion
        return 1.0;
    }
    if((exp(sum)-1.)/sum<0){
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
            
            if(Jperp !=0 ){
                
                E+= -2.*(pow(Jperp,2)/(Temp*KB))*(1.-Spins[i]*Spins[neighbours[i][j]])*(H_func(Spins, neighbours,  i, neighbours[i][j],Temp,+1));//
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
    int posj,posk;
    double H_ij_pi,H_ij_ni;
    double H_jk_pi,H_jk_ni;
    
    
    for(int N=0; N<N_spins; N++){
        /// site to be flipped ///
        site=(rand()%N_spins);

        /// Defining Energy of the flip ///
//        DE=0;
//
//        for(int j=0; j<Z;j++){
//            DE -= 2*J*Spins[site]*Spins[neighbours[site][j]] ;
//
//            if(Jperp!=0){
//                H_pi_j=H_func(Spins, neighbours,  site, neighbours[site][j],T,+1);
//                H_ni_j=H_func(Spins, neighbours,  site, neighbours[site][j],T,-1);
//
//                DE -= 2*(pow(Jperp,2)/(T*KB))*(H_ni_j-H_pi_j);
//                DE -= 2*(pow(Jperp,2)/(T*KB))*(Spins[site]*Spins[neighbours[site][j]])*(H_ni_j+H_pi_j);
//
//            }
        DE=0;
        
        for(int j=0; j<Z;j++){
            posj=neighbours[site][j];//Position of the jth neighbour
            
            DE -= 2*J*Spins[site]*Spins[posj];
            //        cout<<neighbours[site1][j]<<endl;
            //        cout<<"Only J DE= "<<DE<<endl;
            if(Jperp!=0){
                
                //            H_ij_pi=H_func(Spins, neighbours,  site1, posj,T,+1);
                //            H_ij_ni=H_func(Spins, neighbours,  site1, posj,T,-1);
                //
                H_ij_pi=H_DE(Spins, neighbours,  site, posj,T,+1,site);
                H_ij_ni=H_DE(Spins, neighbours,  site, posj,T,-1,site);
                
                
                DE -= 2.*(pow(Jperp,2)/(T*KB))*(H_ij_ni-H_ij_pi);
                
                DE -= 2.*(pow(Jperp,2)/(T*KB))*(Spins[site]*Spins[posj])*(H_ij_ni+H_ij_pi);
                
                for(int k=0; k<Z; k++){
                    posk=neighbours[posj][k];
                    
                    if(posk!=site && Spins[posj]!=Spins[posk] ){
                        
                        H_jk_pi=H_DE(Spins,neighbours, posj,posk,T,+1, site);
                        H_jk_ni=H_DE(Spins,neighbours, posj,posk,T,-1, site);
                        DE -=2.*(pow(Jperp,2)/(T*KB))*(1-Spins[posj]*Spins[posk])*( H_jk_ni-H_jk_pi );
                    }
                }
                //            cout<<"Full, DE= "<<DE<<endl;
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


