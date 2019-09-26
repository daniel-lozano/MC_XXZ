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

using namespace std;

std::random_device rd;  //Will be used to obtain a seed for the random number engine
std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
std::uniform_real_distribution<> dis(0.0, 1.0);

#define Lx 3
#define Ly 3
#define Lz 3
#define z_nn 6
#define L_sub_lattice 4
#define J 1 //negative for Ferromagnetic
#define Jperp 0.01
#define J2 0.2
#define J3a 0.2
#define J3b 0.0

// Set to 1 to have periodic boundary conditions (PBC)
#define PBC 1

struct Lattice{
    // Values for the spins in the lattice
    int value[Lx][Ly][Lz][L_sub_lattice];
    
    // Values for the magnetization of each tetrahedron
    int Magnetization[Lx][Ly][Lz];
    
    // Marking of the tetrahedra
    int Marking[Lx][Ly][Lz];
};

int delta(int x, int y);
double get_energy(Lattice Spin,double T);
void sweep_lattice(Lattice & Spin, double Temp); // Modifies the structure values
int get_fractions(Lattice Spin, int ans);

// Neighbour functions //
double get_1_NN(Lattice Spin,int nx, int ny, int nz, int mu);
double get_2_NN(Lattice Spin,int nx, int ny, int nz, int mu);
double get_3a_NN(Lattice Spin,int nx, int ny, int nz, int mu);
double get_3b_NN(Lattice Spin,int nx, int ny, int nz, int mu);


//Loop algorithm functions//

void init_list(int list[][4], int L);
void print_list(int list[][4], int start ,int L);
int find_index(int list[][4],int L,int x, int y, int z, int mu);
void sweep_closed_loops(Lattice & Spin, double Temp);

// Neighbour functions //
int get_1_NN_loop(Lattice Spin,// Lattice spins
                  int nx_init,// initial position x
                  int ny_init,// initial position y
                  int nz_init,// initial position z
                  int mu_init,// initial position mu
                  int nx_prev,// previous position x
                  int ny_prev,// previous position y
                  int nz_prev,// previous position z
                  int mu_prev,// previous position mu
                  int neighbours_list[][4], // list of neighbours
                  int other// Number of times a tetrahedra has being cosecutively visited
);

int main(){
    
    int lattice_spins[Lx][Ly][Lz][L_sub_lattice];
    int N_thermalice=1000;
    int N_measure=10000;
    
    // Temperature array
    int N_temps=50;
    double T_max=5;
    double T_min=0.5;//0.001;
    double dt=(log(T_max/T_min))/N_temps;
    double* T_array=new double[N_temps];
    double* E_array=new double[N_temps];
    
    // Array where the energy and temperature will be store
    ofstream File_E_T;
    File_E_T.open("Energy_"+ to_string(Lx)+ "_"+ to_string(Ly)+ "_"+ to_string(Lz)+"_H.txt",ios::app);//"_J2_"+"_J3a_"
    
    // Inicializing temperature and Energy arrays
    T_array[0]=T_max;
    E_array[0]=0;
    
    for(int i=1; i<N_temps;i++){
        T_array[i]=T_array[i-1]*exp(-dt);// Done in log scale for ploting
//        cout<< T_array[i] << endl;
        E_array[i]=0;
    }
    cout << "Lx=" << Lx<< endl;
    cout << "Ly=" << Ly<< endl;
    cout << "Lz=" << Lz<< endl;
    cout<< "T_min=" << T_array[0]<<endl;
    cout << "T_max=" << T_array[N_temps-1]<<endl;

    Lattice Spins;
    srand(time(NULL));
    
    cout << "Doing the calculus with a structure" << endl;
    
    int mag_cell=0;
    for(int i=0; i<Lx; i++){
        for(int j=0; j<Ly; j++){
            for(int k=0; k<Lz; k++){
                mag_cell=0;
//                cout<< i<<j<<k << " ";
                for(int mu=0; mu<L_sub_lattice; mu++){
                    Spins.value[i][j][k][mu]=2*(rand()%2)-1;
//                    cout<< Spins.value[i][j][k][mu]<< " ";
                    mag_cell+= Spins.value[i][j][k][mu];
                }
//                cout<< endl;
                Spins.Marking[i][j][k]=0;
                Spins.Magnetization[i][j][k]=mag_cell;
            }
        }
    }
    cout<< "Energy=" << get_energy(Spins,1)/(Lx*Ly*Lz*L_sub_lattice)<<endl; // Energy per site

    cout << "Sweeps\n" << endl;
    
    
    double T,E,E_2;
    double Energy=0;
    int trigger=1;
    
//    for(int i=0;i<10;i++){
//        sweep_closed_loops(Spins, T_array[0]);
//        cout<< "Energy=" << get_energy(Spins)/(Lx*Ly*Lz*L_sub_lattice)<<endl;
//    }
//

//    for(int i=0; i<Lx; i++){
//        for(int j=0; j<Ly; j++){
//            for(int k=0; k<Lz; k++){
//                mag_cell=0;
//                cout<< i<<j<<k << " ";
//                for(int mu=0; mu<L_sub_lattice; mu++){
////                    Spins.value[i][j][k][mu]=2*(rand()%2)-1;
//                    cout<< Spins.value[i][j][k][mu]<< " ";
////                    mag_cell+= Spins.value[i][j][k][mu];
//                }
//                cout<< endl;
////                Spins.Marking[i][j][k]=0;
////                Spins.Magnetization[i][j][k]=mag_cell;
//            }
//        }
//    }
    
    
    int N_sites=Lx*Ly*Lz*L_sub_lattice;
    
//    double s2s2=0;
//    double s2s3=0;
//    double* S2S3=new double[N_temps];
//    for(int i=0; i<N_temps;i++){
//        S2S3[i]=0;// Done in log scale for ploting
//        //        cout<< T_array[i] << endl;
//    }

    
    for(int i=0; i<N_temps; i++){
//        break;
        T=T_array[i];
        E=0;
        E_2=0;
        
//        s2s2=0;
//        s2s3=0;

        get_fractions(Spins, 1);
        cout<< "i=" << i << " out of " << N_temps<< endl;
        for(int j =0; j< N_thermalice; j++){
            sweep_lattice(Spins,T);
            sweep_closed_loops(Spins, T);
        }

        for(int j=0; j< N_measure; j++){
            sweep_lattice(Spins,T);
            sweep_closed_loops(Spins, T);
            Energy=get_energy(Spins,T);
            E+=Energy/N_measure;
            E_2+=pow(Energy,2)/N_measure;
//            s2s2+=(1.*Spins.value[0][0][0][2]*Spins.value[0][1][0][2])/N_measure;//[1][1][1]
//            s2s3+=(1.*Spins.value[0][0][0][2]*Spins.value[0][1][0][3])/N_measure;//[0][1][0]
//
            
            
        }
        E_array[i]=E/(N_sites);
        
        cout<<T << " " << E_array[i]  << endl;
//        cout<<T << " <S2S2>= " << (s2s2)/N_sites  << endl;
//        cout<<T << " <S2S3>= " << (s2s3)/N_sites  << endl;
//        cout<<T << " <S2S2-S2S3>= " << (s2s2-s2s3)/N_sites  << endl;
//
//        File_E_T<< T << " " << E_array[i] << " " << (E_2-pow(E,2))/(N_sites*pow(T,2))<<" "<< (s2s2-s2s3)/N_sites << endl;

        File_E_T<< T << " " << E_array[i] << " " << (E_2-pow(E,2))/(N_sites*pow(T,2))<< endl;

    }
    File_E_T.close();

    

    return 0;
}

//// ------------------------------------------------------------------------------------------------
//// ------------------------------------------------------------------------------------------------
//// --------------------------------------- FUNCTIONS ----------------------------------------------
//// ------------------------------------------------------------------------------------------------
//// ------------------------------------------------------------------------------------------------


int delta(int x, int y){
    int response=0;
    if(x==y){response=1;}
    return response;
}

void sweep_lattice(Lattice & Spin, double Temp){
    
    int N=Lx*Ly*Lz*L_sub_lattice;
    int rx,ry,rz,r_lat;
    int i,j,k,flip,random;
    double DE,prob,alpha;
    double Kb=8.6173E-5;
    
    double J1=J+2*pow(Jperp,2)/Temp-(2.)*(2*z_nn-6)*J*pow(Jperp/Temp,2);
    double Jprime=4*J*pow(Jperp/Temp,2);
    
    
    for(int count=0; count<N; count++){
        
        rx=rand()%Lx;
        ry=rand()%Ly;
        rz=rand()%Lz;
        r_lat=rand()%L_sub_lattice;
        DE=0;
        flip=rand()%2;
        DE -= 2*J1*Spin.value[rx][ry][rz][r_lat]*get_1_NN(Spin,rx, ry, rz, r_lat);// Positive due to antiferromagnetic bond
        
        if(Jprime!=0){
            DE -= 2*Jprime*Spin.value[rx][ry][rz][r_lat]*get_2_NN(Spin,rx, ry, rz, r_lat);
        }
        if(Jprime!=0){
            DE -= 2*Jprime*Spin.value[rx][ry][rz][r_lat]*get_3a_NN(Spin,rx, ry, rz, r_lat);
        }
        if(J3b!=0){
            DE -= 2*J3b*Spin.value[rx][ry][rz][r_lat]*get_3b_NN(Spin,rx, ry, rz, r_lat);
        }
        prob=exp(-DE/(Temp));
        alpha=dis(gen);
        
        if(alpha<prob || DE<0){
            Spin.value[rx][ry][rz][r_lat]*=-1;// New state is better...
            Spin.Magnetization[rx][ry][rz]=Spin.value[rx][ry][rz][0]+
                                           Spin.value[rx][ry][rz][1]+
                                           Spin.value[rx][ry][rz][2]+
            Spin.value[rx][ry][rz][3];
        }
        
    }
    
}

double get_energy(Lattice Spin,double T){
    double sum=0;
    
    double J1=J+2*pow(Jperp,2)/T-(2.)*(2*z_nn-6)*J*pow(Jperp/T,2);
    double Jprime=4*J*pow(Jperp/T,2);
    
    for(int i=0; i<Lx; i++){
        for(int j=0; j<Ly; j++){
            for(int k=0; k<Lz; k++){
                for(int mu=0; mu<L_sub_lattice; mu++ ){
                    sum += J1*Spin.value[i][j][k][mu]*get_1_NN(Spin,i, j, k,mu);
                    if(Jperp!=0){sum+=Jprime*Spin.value[i][j][k][mu]*get_2_NN(Spin,i, j, k,mu);}
                    if(Jperp!=0){sum+=Jprime*Spin.value[i][j][k][mu]*get_3a_NN(Spin,i, j, k,mu);}
                    if(J3b!=0){sum+=J3b*Spin.value[i][j][k][mu]*get_3b_NN(Spin,i, j, k,mu);}

                }
            }
        }
    }
    return sum*0.5;// The extra factor of 1/2 is added due to the over counting
    
}

int get_fractions(Lattice Spin, int ans){
    int mag0=0;
    int mag2=0;
    int mag4=0;
    int sum;
    for(int i=0; i<Lx; i++){
        for(int j=0; j<Ly; j++){
            for(int k=0; k<Lz; k++){
                
                sum=Spin.value[i][j][k][0]+Spin.value[i][j][k][1]+Spin.value[i][j][k][2]+Spin.value[i][j][k][3];
                
                if(abs(sum)==0){
                    mag0+=1;
                }
                else if(abs(sum)==2){
                    mag2+=1;
                }
                else if(abs(sum)==4){
                    mag4+=1;
                }
                else{ cout<<"Sum="<<sum<< endl;}
                
            }
        }
    }
    if(ans==1){
        cout<< "num 0= "<< mag0<< ", num 2= "<< mag2 << ", num 4= "<< mag4<<endl;
    }
    return mag0;
}
//------------------------------------------------------------------------------
//------------- First, second and Third nearest Neighbors functions-------------
//------------------------------------------------------------------------------
double get_1_NN(Lattice Spin,int nx, int ny, int nz, int mu){
    double sum=0;
//    cout << "Starting Spin position" << endl;
//    cout<< (nx)<< " " << (ny) << " " << (nz) << " "<< mu<< endl;
//    cout << "Neighbours" << endl;

    int dx=0;
    int dy=0;
    int dz=0;
//    cout << "1 Neighbours"<< endl;
    for(int nu=0; nu<L_sub_lattice; nu++){
        
        // Define the advance in the lattice given a lattice site
        dx=0;
        dy=0;
        dz=0;
        
        if(mu!=nu){
            //Conditions for the movement in dx
            
            if(mu==0 && nu==1){dx=-1;}
            if(mu==0 && nu==2){dy=-1;}
            if(mu==0 && nu==3){dz=-1;}
            
            if(mu==1){dx=1;}
            if(nu==1){dx=-1;}
            
            if(mu==2){dy=1;}
            if(nu==2){dy=-1;}
            
            if(mu==3){dz=1;}
            if(nu==3){dz=-1;}
            
     
//             Sums the neighbours in the same pyrochlore and the spins in the lower (upper) pyrochlore
            
            // open boundary conditions
            if(PBC==1){
                sum+=Spin.value[(nx+dx+Lx)%Lx][(ny+dy+Ly)%Ly][(nz+dz+Lz)%Lz][nu]+Spin.value[nx][ny][nz][nu];
//                cout<< (nx+dx+Lx)%Lx << " " << (ny+dy+Ly)%Ly << " " << (nz+dz+Lz)%Lz << " "<< nu<< endl;
            }
            
            else{
                if((nx+dx)>=0 && (ny+dy)>=0 && (nz+dz)>=0 && (nx+dx)<Lx && (ny+dy)<Ly && (nz+dz)<Lz){
//                    cout<<"nx+dx="<<nx+dx<<endl;
//                    cout<<"ny+dy="<<ny+dy<<endl;
//                    cout<<"nz+dz="<<nz+dz<<endl;
//                    cout<< (nx+dx+Lx)<< " " << (ny+dy+Ly) << " " << (nz+dz+Lz) << " "<< nu<< endl;
//
//                    cout<< (nx+dx+Lx)%Lx << " " << (ny+dy+Ly)%Ly << " " << (nz+dz+Lz)%Lz << " "<< nu<< " S="<< Spin.value[(nx+dx+Lx)%Lx][(ny+dy+Ly)%Ly][(nz+dz+Lz)%Lz][nu]<< endl;
//                    cout << endl;
//
                    sum+=Spin.value[(nx+dx+Lx)%Lx][(ny+dy+Ly)%Ly][(nz+dz+Lz)%Lz][nu];
                }
                sum+=Spin.value[nx][ny][nz][nu];
//                cout<< (nx)%Lx << " " << (ny)%Ly << " " << (nz)%Lz << " "<< nu<<" S="<< Spin.value[nx][ny][nz][nu]<< endl;
            }
            
        }
    }
//    cout<< "Sum="<<sum << endl;
    return sum;
}


double get_2_NN(Lattice Spin,int nx, int ny, int nz, int mu){
    double sum=0;
    int dx1,dx2,dx3,dx4;
    int dy1,dy2,dy3,dy4;
    int dz1,dz2,dz3,dz4;
    
//    cout << "2 Neighbours"<< endl;
//    cout << "Starting Spin position" << endl;
//    cout<< (nx)%Lx << " " << (ny)%Ly << " " << (nz)%Lz << " "<< mu<< endl;
//    cout << endl;
    
    for(int nu=0; nu<L_sub_lattice; nu++){
        
        if(nu!=mu){
            
            // Note: d\vec{x}_{1}=-d\vec{x}_{3} and d\vec{x}_{2}=-d\vec{x}_{4}, use for optimization later.
            dx1=delta(mu,1) - delta(mu,2)*delta(nu,3)- delta(mu,3)*delta(nu,2)-delta(mu,3)*delta(nu,0)-delta(mu,0)*delta(nu,3);
            dy1=delta(mu,2) - delta(mu,1)*delta(nu,3)- delta(mu,3)*delta(nu,1)-delta(mu,1)*delta(nu,0)-delta(mu,0)*delta(nu,1);
            dz1=delta(mu,3) - delta(mu,1)*delta(nu,2)- delta(mu,2)*delta(nu,1)-delta(mu,2)*delta(nu,0)-delta(mu,0)*delta(nu,2);
            
            dx2=delta(mu,1)-delta(mu,2)*delta(nu,0)-delta(mu,0)*delta(nu,2);
            dy2=delta(mu,2)-delta(mu,3)*delta(nu,0)-delta(mu,0)*delta(nu,3);
            dz2=delta(mu,3)-delta(mu,1)*delta(nu,0)-delta(mu,0)*delta(nu,1);
            
            dx3=-delta(nu,1)+delta(mu,2)*delta(nu,3)+delta(mu,3)*delta(nu,2)+delta(mu,3)*delta(nu,0)+delta(mu,0)*delta(nu,3);
            dy3=-delta(nu,2)+delta(mu,1)*delta(nu,3)+delta(mu,3)*delta(nu,1)+delta(mu,1)*delta(nu,0)+delta(mu,0)*delta(nu,1);
            dz3=-delta(nu,3)+delta(mu,1)*delta(nu,2)+delta(mu,2)*delta(nu,1)+delta(mu,2)*delta(nu,0)+delta(mu,0)*delta(nu,2);
            
            dx4=-delta(nu,1)+delta(mu,2)*delta(nu,0)+delta(mu,0)*delta(nu,2);
            dy4=-delta(nu,2)+delta(mu,3)*delta(nu,0)+delta(mu,0)*delta(nu,3);
            dz4=-delta(nu,3)+delta(mu,1)*delta(nu,0)+delta(mu,0)*delta(nu,1);
            
            
            if(PBC==1){
                sum+=Spin.value[(nx+dx1+Lx)%Lx][(ny+dy1+Ly)%Ly][(nz+dz1+Lz)%Lz][nu];
                sum+=Spin.value[(nx+dx2+Lx)%Lx][(ny+dy2+Ly)%Ly][(nz+dz2+Lz)%Lz][nu];
                sum+=Spin.value[(nx+dx3+Lx)%Lx][(ny+dy3+Ly)%Ly][(nz+dz3+Lz)%Lz][nu];
                sum+=Spin.value[(nx+dx4+Lx)%Lx][(ny+dy4+Ly)%Ly][(nz+dz4+Lz)%Lz][nu];
                
            }
            
            else{ // Only within the system interactions
                if((nx+dx1)>=0 && (nx+dx1)< Lx && (ny+dy1)>=0 && (ny+ dy1)< Ly && (nz+dz1)>=0 && (nz+dz1)< Lz){
                    sum+=Spin.value[(nx+dx1+Lx)%Lx][(ny+dy1+Ly)%Ly][(nz+dz1+Lz)%Lz][nu];
//                    cout << (nx+dx1+Lx)%Lx << " " << (ny+dy1+Ly)%Ly<< " " << (nz+dz1+Lz)%Lz<< " " << nu<< endl;
                }
                if((nx+dx2)>=0 && (nx+dx2)< Lx && (ny+dy2)>=0 && (ny+ dy2)< Ly && (nz+dz2)>=0 && (nz+dz2)< Lz){
                    sum+=Spin.value[(nx+dx2+Lx)%Lx][(ny+dy2+Ly)%Ly][(nz+dz2+Lz)%Lz][nu];
//                    cout << (nx+dx2+Lx)%Lx << " " << (ny+dy2+Ly)%Ly<< " " << (nz+dz2+Lz)%Lz<< " " << nu<< endl;
                }
                if((nx+dx3)>=0 && (nx+dx3)< Lx && (ny+dy3)>=0 && (ny+ dy3)< Ly && (nz+dz3)>=0 && (nz+dz3)< Lz){
                    sum+=Spin.value[(nx+dx3+Lx)%Lx][(ny+dy3+Ly)%Ly][(nz+dz3+Lz)%Lz][nu];
//                    cout << (nx+dx3+Lx)%Lx << " " << (ny+dy3+Ly)%Ly<< " " << (nz+dz3+Lz)%Lz<< " " << nu<< endl;
                }
                if((nx+dx4)>=0 && (nx+dx4)< Lx && (ny+dy4)>=0 && (ny+ dy4)< Ly && (nz+dz4)>=0 && (nz+dz4)< Lz){
                    sum+=Spin.value[(nx+dx4+Lx)%Lx][(ny+dy4+Ly)%Ly][(nz+dz4+Lz)%Lz][nu];
//                    cout << (nx+dx4+Lx)%Lx << " " << (ny+dy4+Ly)%Ly<< " " << (nz+dz4+Lz)%Lz<< " " << nu<< endl;
                }
            }
            
            ////////// Printing neighbouring lattice sites
//            cout << "["<< nx+dx1 << "," << ny+dy1 << ","<< nz+dz1<< "]"<< endl;
//            cout << "["<< nx+dx2 << "," << ny+dy2 << ","<< nz+dz2<< "]"<< endl;
//            cout << "["<< nx+dx3 << "," << ny+dy3 << ","<< nz+dz3 << "]"<< endl;
//            cout << "["<< nx+dx4 << "," << ny+dy4 << ","<< nz+dz4 << "]" << endl;
        }
        
    }
//    cout << "End of 2 Neighbours"<< endl;
//    cout << endl;

    return sum;
}

double get_3a_NN(Lattice Spin,int nx, int ny, int nz, int mu){
    double sum=0;
    
    int dx1,dx2,dx3,dx4,dx5,dx6;
    int dy1,dy2,dy3,dy4,dy5,dy6;
    int dz1,dz2,dz3,dz4,dz5,dz6;
    
//
//    cout << "Starting Spin position" << endl;
//    cout<< (nx)%Lx << " " << (ny)%Ly << " " << (nz)%Lz << " "<< mu<< endl;
//    cout << "3a Neighbours"<< endl;
//
    dx1=delta(mu,0)+delta(mu,1);
    dy1=delta(mu,2);
    dz1=delta(mu,3);
    
    dx2=-delta(mu,0)-delta(mu,1);
    dy2=-delta(mu,2);
    dz2=-delta(mu,3);
    
    dx3=delta(mu,1)-delta(mu,2);
    dy3=delta(mu,0)-delta(mu,1)+delta(mu,2)-delta(mu,3);
    dz3=delta(mu,3);
    
    dx4=delta(mu,1)-delta(mu,3);
    dy4=-delta(mu,0)+delta(mu,2);
    dz4=-delta(mu,1)-delta(mu,2)+delta(mu,3);
    
    dx5=-delta(mu,1)+delta(mu,2)+delta(mu,3);
    dy5=delta(mu,1)-delta(mu,2);
    dz5=delta(mu,0)-delta(mu,3);
    
    dx6=-delta(mu,1);
    dy6=-delta(mu,2)+delta(mu,3);
    dz6=-delta(mu,0)+delta(mu,1)+delta(mu,2)-delta(mu,3);
    
    
    if(PBC==1){
        sum+=Spin.value[(nx+dx1+Lx)%Lx][(ny+dy1+Ly)%Ly][(nz+dz1+Lz)%Lz][mu];
        sum+=Spin.value[(nx+dx2+Lx)%Lx][(ny+dy2+Ly)%Ly][(nz+dz2+Lz)%Lz][mu];
        sum+=Spin.value[(nx+dx3+Lx)%Lx][(ny+dy3+Ly)%Ly][(nz+dz3+Lz)%Lz][mu];
        sum+=Spin.value[(nx+dx4+Lx)%Lx][(ny+dy4+Ly)%Ly][(nz+dz4+Lz)%Lz][mu];
        sum+=Spin.value[(nx+dx5+Lx)%Lx][(ny+dy5+Ly)%Ly][(nz+dz5+Lz)%Lz][mu];
        sum+=Spin.value[(nx+dx6+Lx)%Lx][(ny+dy6+Ly)%Ly][(nz+dz6+Lz)%Lz][mu];
    }
    
    else{
        // Only within the system interactions
        if((nx+dx1)>=0 && (nx+dx1)< Lx && (ny+dy1)>=0 && (ny+ dy1)< Ly && (nz+dz1)>=0 && (nz+dz1)< Lz){
            sum+=Spin.value[(nx+dx1+Lx)%Lx][(ny+dy1+Ly)%Ly][(nz+dz1+Lz)%Lz][mu];
//            cout << (nx+dx1+Lx)%Lx << " " << (ny+dy1+Ly)%Ly<< " " << (nz+dz1+Lz)%Lz<< " " << mu<< endl;
        }
        // Only within the system interactions
        if((nx+dx2)>=0 && (nx+dx2)< Lx && (ny+dy2)>=0 && (ny+ dy2)< Ly && (nz+dz2)>=0 && (nz+dz2)< Lz){
            sum+=Spin.value[(nx+dx2+Lx)%Lx][(ny+dy2+Ly)%Ly][(nz+dz2+Lz)%Lz][mu];
//            cout << (nx+dx2+Lx)%Lx << " " << (ny+dy2+Ly)%Ly<< " " << (nz+dz2+Lz)%Lz<< " " << mu<< endl;
        }
        // Only within the system interactions
        if((nx+dx3)>=0 && (nx+dx3)< Lx && (ny+dy3)>=0 && (ny+ dy3)< Ly && (nz+dz3)>=0 && (nz+dz3)< Lz){
            sum+=Spin.value[(nx+dx3+Lx)%Lx][(ny+dy3+Ly)%Ly][(nz+dz3+Lz)%Lz][mu];
//            cout << (nx+dx3+Lx)%Lx << " " << (ny+dy3+Ly)%Ly<< " " << (nz+dz3+Lz)%Lz<< " " << mu<< endl;
        }
        // Only within the system interactions
        if((nx+dx4)>=0 && (nx+dx4)< Lx && (ny+dy4)>=0 && (ny+ dy4)< Ly && (nz+dz4)>=0 && (nz+dz4)< Lz){
            sum+=Spin.value[(nx+dx4+Lx)%Lx][(ny+dy4+Ly)%Ly][(nz+dz4+Lz)%Lz][mu];
//            cout << (nx+dx4+Lx)%Lx << " " << (ny+dy4+Ly)%Ly<< " " << (nz+dz4+Lz)%Lz<< " " << mu<< endl;
        }
        if((nx+dx5)>=0 && (nx+dx5)< Lx && (ny+dy5)>=0 && (ny+ dy5)< Ly && (nz+dz5)>=0 && (nz+dz5)< Lz){
            sum+=Spin.value[(nx+dx5+Lx)%Lx][(ny+dy5+Ly)%Ly][(nz+dz5+Lz)%Lz][mu];
//            cout << (nx+dx5+Lx)%Lx << " " << (ny+dy5+Ly)%Ly<< " " << (nz+dz5+Lz)%Lz<< " " << mu<< endl;
        }
        if((nx+dx6)>=0 && (nx+dx6)< Lx && (ny+dy6)>=0 && (ny+ dy6)< Ly && (nz+dz6)>=0 && (nz+dz6)< Lz){
            sum+=Spin.value[(nx+dx6+Lx)%Lx][(ny+dy6+Ly)%Ly][(nz+dz6+Lz)%Lz][mu];
//            cout << (nx+dx6+Lx)%Lx << " " << (ny+dy6+Ly)%Ly<< " " << (nz+dz6+Lz)%Lz<< " " << mu<< endl;
        }
    }
    
//    cout << "End of 3a Neighbours"<< endl;
    return sum;
}

double get_3b_NN(Lattice Spin,int nx, int ny, int nz, int mu){
    double sum=0;
    
    int dx1,dx2,dx3,dx4,dx5,dx6;
    int dy1,dy2,dy3,dy4,dy5,dy6;
    int dz1,dz2,dz3,dz4,dz5,dz6;
    
    //
    //    cout << "Starting Spin position" << endl;
    //    cout<< (nx)%Lx << " " << (ny)%Ly << " " << (nz)%Lz << " "<< mu<< endl;
    //    cout << "3a Neighbours"<< endl;
    //
    dx1=delta(mu,2)+delta(mu,3);
    dy1=-delta(mu,0)-delta(mu,1);
    dz1=delta(mu,0)+delta(mu,1);
    
    dx2=-delta(mu,2)-delta(mu,3);
    dy2=delta(mu,0)+delta(mu,1);
    dz2=-delta(mu,0)+delta(mu,1);
    
    dx3=delta(mu,0);
    dy3=-delta(mu,0)+delta(mu,3);
    dz3=delta(mu,1)+delta(mu,2);
    
    dx4=delta(mu,0);
    dy4=delta(mu,1)-delta(mu,3);
    dz4=-delta(mu,0)-delta(mu,2);
    
    dx5=-delta(mu,0)-delta(mu,2)+delta(mu,3);
    dy5=delta(mu,0)-delta(mu,1)-delta(mu,3);
    dz5=delta(mu,2);
    
    dx6=-delta(mu,0)+delta(mu,2)-delta(mu,3);
    dy6=delta(mu,3);
    dz6=delta(mu,0)-delta(mu,1)-delta(mu,2);
    
    
    if(PBC==1){
        sum+=Spin.value[(nx+dx1+Lx)%Lx][(ny+dy1+Ly)%Ly][(nz+dz1+Lz)%Lz][mu];
        sum+=Spin.value[(nx+dx2+Lx)%Lx][(ny+dy2+Ly)%Ly][(nz+dz2+Lz)%Lz][mu];
        sum+=Spin.value[(nx+dx3+Lx)%Lx][(ny+dy3+Ly)%Ly][(nz+dz3+Lz)%Lz][mu];
        sum+=Spin.value[(nx+dx4+Lx)%Lx][(ny+dy4+Ly)%Ly][(nz+dz4+Lz)%Lz][mu];
        sum+=Spin.value[(nx+dx5+Lx)%Lx][(ny+dy5+Ly)%Ly][(nz+dz5+Lz)%Lz][mu];
        sum+=Spin.value[(nx+dx6+Lx)%Lx][(ny+dy6+Ly)%Ly][(nz+dz6+Lz)%Lz][mu];
    }
    
    else{
        // Only within the system interactions
        if((nx+dx1)>=0 && (nx+dx1)< Lx && (ny+dy1)>=0 && (ny+ dy1)< Ly && (nz+dz1)>=0 && (nz+dz1)< Lz){
            sum+=Spin.value[(nx+dx1+Lx)%Lx][(ny+dy1+Ly)%Ly][(nz+dz1+Lz)%Lz][mu];
            //            cout << (nx+dx1+Lx)%Lx << " " << (ny+dy1+Ly)%Ly<< " " << (nz+dz1+Lz)%Lz<< " " << mu<< endl;
        }
        // Only within the system interactions
        if((nx+dx2)>=0 && (nx+dx2)< Lx && (ny+dy2)>=0 && (ny+ dy2)< Ly && (nz+dz2)>=0 && (nz+dz2)< Lz){
            sum+=Spin.value[(nx+dx2+Lx)%Lx][(ny+dy2+Ly)%Ly][(nz+dz2+Lz)%Lz][mu];
            //            cout << (nx+dx2+Lx)%Lx << " " << (ny+dy2+Ly)%Ly<< " " << (nz+dz2+Lz)%Lz<< " " << mu<< endl;
        }
        // Only within the system interactions
        if((nx+dx3)>=0 && (nx+dx3)< Lx && (ny+dy3)>=0 && (ny+ dy3)< Ly && (nz+dz3)>=0 && (nz+dz3)< Lz){
            sum+=Spin.value[(nx+dx3+Lx)%Lx][(ny+dy3+Ly)%Ly][(nz+dz3+Lz)%Lz][mu];
            //            cout << (nx+dx3+Lx)%Lx << " " << (ny+dy3+Ly)%Ly<< " " << (nz+dz3+Lz)%Lz<< " " << mu<< endl;
        }
        // Only within the system interactions
        if((nx+dx4)>=0 && (nx+dx4)< Lx && (ny+dy4)>=0 && (ny+ dy4)< Ly && (nz+dz4)>=0 && (nz+dz4)< Lz){
            sum+=Spin.value[(nx+dx4+Lx)%Lx][(ny+dy4+Ly)%Ly][(nz+dz4+Lz)%Lz][mu];
            //            cout << (nx+dx4+Lx)%Lx << " " << (ny+dy4+Ly)%Ly<< " " << (nz+dz4+Lz)%Lz<< " " << mu<< endl;
        }
        if((nx+dx5)>=0 && (nx+dx5)< Lx && (ny+dy5)>=0 && (ny+ dy5)< Ly && (nz+dz5)>=0 && (nz+dz5)< Lz){
            sum+=Spin.value[(nx+dx5+Lx)%Lx][(ny+dy5+Ly)%Ly][(nz+dz5+Lz)%Lz][mu];
            //            cout << (nx+dx5+Lx)%Lx << " " << (ny+dy5+Ly)%Ly<< " " << (nz+dz5+Lz)%Lz<< " " << mu<< endl;
        }
        if((nx+dx6)>=0 && (nx+dx6)< Lx && (ny+dy6)>=0 && (ny+ dy6)< Ly && (nz+dz6)>=0 && (nz+dz6)< Lz){
            sum+=Spin.value[(nx+dx6+Lx)%Lx][(ny+dy6+Ly)%Ly][(nz+dz6+Lz)%Lz][mu];
            //            cout << (nx+dx6+Lx)%Lx << " " << (ny+dy6+Ly)%Ly<< " " << (nz+dz6+Lz)%Lz<< " " << mu<< endl;
        }
    }
    
    //    cout << "End of 3a Neighbours"<< endl;
    return sum;
}


//------------------------------------------------------------------------------
//---------------------------------- Loop algorithm------------------------------
//------------------------------------------------------------------------------

void sweep_closed_loops(Lattice & Spin, double Temp){
    
    //-----------------------------------------------------------------------------
    //------------------    getting path of 1nn   ---------------------------------
    //-----------------------------------------------------------------------------
    
    int neighbours_list[z_nn][4];
    int lattice_list[Lx*Ly*Lz*L_sub_lattice][4];// x,y,z,mu
    
    //    Initialicing neighbours list
    init_list(neighbours_list,z_nn);
    init_list(lattice_list,Lx*Ly*Lz*L_sub_lattice);
    
//    print_list(neighbours_list, 0 ,z_nn);
    
    int current_x,current_y,current_z,current_mu; //current position of the path
    int last_x,last_y,last_z,last_mu; // Last position of the path
    int paths; // posible paths to choose in a site
    int path_chosen; // Random path chosen
    int index; // Index where the path starts
    int num_points=1; // Binary variable counting times a tetrahedro have being accessed
    int length=0;
    current_x=rand()%Lx;//0;
    current_y=rand()%Ly;//0;
    current_z=rand()%Lz;//0;
    current_mu=rand()%L_sub_lattice;//3;
//    cout<< "Initial position=" << current_x << current_y << current_z << current_mu <<endl;
    
    //Starting position should be the same
    last_x=current_x;
    last_y=current_y;
    last_z=current_z;
    last_mu=current_mu;
    
//    cout<< "\nInitiating loop algorithm!\n"<<endl;
    
    for(int i=0;i<Lx*Ly*Lz*L_sub_lattice;i++){
        
        paths=get_1_NN_loop(Spin, current_x, current_y, current_z, current_mu, last_x, last_y, last_z, last_mu, neighbours_list,num_points);// Gets number of possible paths and updates the neighbours list
        
        index=find_index(lattice_list,Lx*Ly*Lz*L_sub_lattice,current_x, current_y, current_z, current_mu);//Finds a repeated position
        
        if(index!=-1){
            //            cout<< "Path have being closed"<<endl;
            //            cout<<"Index="<<index<<endl;
            break;
        }
        
        lattice_list[i][0]=current_x;
        lattice_list[i][1]=current_y;
        lattice_list[i][2]=current_z;
        lattice_list[i][3]=current_mu;
        length++;// Adding a site to the path
        num_points=1;
        
        if(paths==0){
            //            cout<<"Loop interrupted! No closed loop found"<<endl;
            break;
        }
        path_chosen=rand()%paths;
        
        //        cout<< "Possible paths"<< endl;
        //        print_list(neighbours_list,0,z_nn);
        //
        //        cout<< "Path chosen="<< path_chosen<<endl;
        //        cout<< neighbours_list[path_chosen][0]<< " ";
        //        cout<< neighbours_list[path_chosen][1]<< " ";
        //        cout<< neighbours_list[path_chosen][2]<< " ";
        //        cout<< neighbours_list[path_chosen][3]<< endl;
        //
        // Adding the position to the path
        last_x=current_x;
        last_y=current_y;
        last_z=current_z;
        last_mu=current_mu;
        
        // Getting the new starting position
        current_x=neighbours_list[path_chosen][0];
        current_y=neighbours_list[path_chosen][1];
        current_z=neighbours_list[path_chosen][2];
        current_mu=neighbours_list[path_chosen][3];
        
        if(last_x==current_x && last_y==current_y && last_z==current_z){
            num_points=2;
        }
    }
    
    //Getting the loop positions
    if(paths!=0){//
    double DE=0;
    double beta=dis(gen);
    double prob=0;

    if(paths==0){
        index=0;
        
    }//if an open path is found the starting point is the first random site
    
    //Going through the path
    for(int i=index; i<length;i++){
        
        if(lattice_list[i][0]==-1){
            break;
        }
        DE -= 2*J*Spin.value[lattice_list[i][0]][lattice_list[i][1]][lattice_list[i][2]][lattice_list[i][3]]*get_1_NN(Spin,lattice_list[i][0], lattice_list[i][1], lattice_list[i][2], lattice_list[i][3]);// Positive due to antiferromagnetic bond
        
        if(J2!=0){//Getting the energy term for every path, Second neighbours
            DE -= 2*J2*Spin.value[lattice_list[i][0]][lattice_list[i][1]][lattice_list[i][2]][lattice_list[i][3]]*get_2_NN(Spin,lattice_list[i][0], lattice_list[i][1], lattice_list[i][2], lattice_list[i][3]);
        }
        if(J3a!=0){ //Getting the energy term for every path, Third neighbours
            DE-= 2*J3a*Spin.value[lattice_list[i][0]][lattice_list[i][1]][lattice_list[i][2]][lattice_list[i][3]]*get_3a_NN(Spin,lattice_list[i][0], lattice_list[i][1], lattice_list[i][2], lattice_list[i][3]);
        }
        if(J3b!=0){ //Getting the energy term for every path, Third neighbours
            DE-= 2*J3b*Spin.value[lattice_list[i][0]][lattice_list[i][1]][lattice_list[i][2]][lattice_list[i][3]]*get_3b_NN(Spin,lattice_list[i][0], lattice_list[i][1], lattice_list[i][2], lattice_list[i][3]);
        }
        //Flipping path to check the DE produce by it
        Spin.value[lattice_list[i][0]][lattice_list[i][1]][lattice_list[i][2]][lattice_list[i][3]]*=-1;
       
    }
    prob=exp(-DE/Temp);
//    cout<< "Prob="<<prob<< endl;
//    cout<< "beta="<< beta<<endl;
//    cout<< "T="<<Temp<<endl;
//    cout << "Lenght=" << length <<endl;
//    cout<< "DE=" << DE/(Lx*Ly*Lz*L_sub_lattice)<< endl;
//    cout<< "Index="<<index<< endl;
    //Unflip the path if the Energy increases or the MC step doesnt permit it
    if( prob<beta ){//
        for(int i=index; i<length;i++){
            
            if(lattice_list[i][0]==-1){
                break;
            }
            //Unflipping
            Spin.value[lattice_list[i][0]][lattice_list[i][1]][lattice_list[i][2]][lattice_list[i][3]]*=-1;
            
        }
        
    }
}
    
    
    
}




void init_list(int list[][4], int L){
    for(int i=0;i<L;i++){
        list[i][0]=-1;
        list[i][1]=-1;
        list[i][2]=-1;
        list[i][3]=-1;
    }
}
void print_list(int list[][4], int start ,int L){
    for(int i=start;i<L;i++){
        if(list[i][0]==-1){
            break;
        }
        cout<< list[i][0]<< " ";
        cout<< list[i][1]<<  " ";
        cout<< list[i][2]<< " ";
        cout<< list[i][3]<< endl;
    }
}

int find_index(int list[][4],int L,int x, int y, int z, int mu){
    int index=-1;
    for(int i=0;i<L;i++){
        if(list[i][0]==x && list[i][1]==y && list[i][2]==z && list[i][3]==mu ){
            index=i;
            break;
        }
    }
    return index;
}

//Working code! Finds possible paths to cross

int get_1_NN_loop(Lattice Spin, int nx_init, int ny_init, int nz_init, int mu_init,int nx_prev, int ny_prev, int nz_prev, int mu_prev, int neighbours_list[][4], int other){
    
    //    cout << "Starting Spin position" << endl;
    //    cout<< (nx_init)<< " " << (ny_init) << " " << (nz_init) << " "<< mu_init<< "= "<< Spin.value[nx_init][ny_init][nz_init][mu_init] << endl;
    //    cout<< "Previous Spin position" << endl;
    //    cout<< (nx_prev)<< " " << (ny_prev) << " " << (nz_prev) << " "<< mu_prev<< "= "<< Spin.value[nx_prev][ny_prev][nz_prev][mu_prev] << endl;
    //
    //    cout << "Neighbours " << endl;
    
    int dx=0;
    int dy=0;
    int dz=0;
    int counter=0;
    //    cout << "1 Neighbours"<< endl;
    for(int nu=0; nu<L_sub_lattice; nu++){
        
        // Define the advance in the lattice given a lattice site
        dx=0;
        dy=0;
        dz=0;
        
        if(mu_init!=nu){
            //Conditions for the movement in dx
            
            if(mu_init==0 && nu==1){dx=-1;}
            if(mu_init==0 && nu==2){dy=-1;}
            if(mu_init==0 && nu==3){dz=-1;}
            
            if(mu_init==1){dx=1;}
            if(nu==1){dx=-1;}
            
            if(mu_init==2){dy=1;}
            if(nu==2){dy=-1;}
            
            if(mu_init==3){dz=1;}
            if(nu==3){dz=-1;}
            
            
            
            // open boundary conditions
            if(PBC==1 ){
                
                if(
                   ( (nx_init+dx+Lx)%Lx!=nx_prev || (ny_init+dy+Ly)%Ly!=ny_prev  || (nz_init+dz+Lz)%Lz!=nz_prev || nu!=mu_prev) &&//it's not the previous position
                   Spin.value[(nx_init+dx+Lx)%Lx][(ny_init+dy+Ly)%Ly][(nz_init+dz+Lz)%Lz][nu]*
                   Spin.value[nx_init][ny_init][nz_init][mu_init]==-1 &&//it's an inverse spin
                   other!=1) {
                    
                    neighbours_list[counter][0]=(nx_init+dx+Lx)%Lx;
                    neighbours_list[counter][1]=(ny_init+dy+Ly)%Ly;
                    neighbours_list[counter][2]=(nz_init+dz+Lz)%Lz;
                    neighbours_list[counter][3]=nu;
                    counter++;
                }
                
                if(
                   (nx_init != nx_prev || ny_init!=ny_prev || nz_init!=nz_prev  ||nu!=mu_prev) && //it's not the previous position
                   Spin.value[nx_init][ny_init][nz_init][nu]*Spin.value[nx_init][ny_init][nz_init][mu_init]==-1 &&//it's an inverse spin
                   other!=2
                   ){
                    
                    neighbours_list[counter][0]=nx_init;
                    neighbours_list[counter][1]=ny_init;
                    neighbours_list[counter][2]=nz_init;
                    neighbours_list[counter][3]=nu;
                    counter++;
                    
                }
                
            }
            
        }
    }
    return counter;
}

