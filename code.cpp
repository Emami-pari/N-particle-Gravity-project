#include  <iostream>
#include  <cmath>                          
#include  <cstdlib>                       
#include  <unistd.h>                       
#include <fstream>
#include <stdio.h>
#include <string.h>
#include <time.h>
using namespace std;
typedef double  real;                     
const int NDIM = 3;             
real  dt_tot =30;

void correct_step(real pos[][NDIM], real vel[][NDIM], 
                  const real acc[][NDIM], const real jerk[][NDIM],
                  const real old_pos[][NDIM], const real old_vel[][NDIM], 
                  const real old_acc[][NDIM], const real old_jerk[][NDIM],
                  int n, real dt);
void evolve(const real mass[], real pos[][NDIM], real vel[][NDIM],real v[],real p[],int n, real & t, real dt_param, real dt_dia, real dt_out,
            real dt_tot, bool init_out, bool x_flag);
void evolve_step(const real mass[], real pos[][NDIM], real vel[][NDIM],
                 real acc[][NDIM], real jerk[][NDIM], int n, real & t,
                 real dt, real & epot, real & coll_time);
void get_acc_jerk_pot_coll(const real mass[], const real pos[][NDIM],
                           const real vel[][NDIM], real acc[][NDIM],
                           real jerk[][NDIM], int n, real & epot,
                           real & coll_time);
void get_snapshot(real mass[], real pos[][NDIM], real vel[][NDIM],real v[],real p[], int n);
void predict_step(real pos[][NDIM], real vel[][NDIM], 
                  const real acc[][NDIM], const real jerk[][NDIM],
                  int n, real dt);
void put_snapshot(const real mass[], const real pos[][NDIM],
                  const real vel[][NDIM], real v[],real p[],int n, real t);
bool read_options(int argc, char *argv[], real & dt_param, real & dt_dia,
                  real & dt_out, real & dt_tot, bool & i_flag, bool & x_flag);
void write_diagnostics(const real mass[], const real pos[][NDIM],
                       const real vel[][NDIM], const real acc[][NDIM],
                       const real jerk[][NDIM], int n, real t, real epot,
                       int nsteps, real & einit, bool init_flag,
                       bool x_flag);



int main(int argc, char *argv[])
{
    real  dt_param = 0.03;    
    real  dt_dia = 1;         
    real  dt_out = 1;          
    bool  init_out = false;  
    bool  x_flag = false;     

    if (! read_options(argc, argv, dt_param, dt_dia, dt_out, dt_tot, init_out,
                       x_flag))
        return 1;                

    int n;                      
    cout<<"enter the number of body:";cin >> n;
    real t;                   
    cout<<"enter time:";cin >> t;

    real * mass = new real[n];             
    real (* pos)[NDIM] = new real[n][NDIM];  
    real (* vel)[NDIM] = new real[n][NDIM];    
    real (* v) = new real[n]; 
    real (* p)=new real[n];
    get_snapshot(mass, pos, vel,v,p, n);
   
    evolve(mass, pos, vel, v,p, n, t, dt_param, dt_dia, dt_out, dt_tot, init_out,
           x_flag);

    delete[] mass;
    delete[] pos;
    delete[] vel;
    delete[] v;
    delete[] p;
}


bool read_options(int argc, char *argv[], real & dt_param, real & dt_dia,
                  real & dt_out, real & dt_tot, bool & i_flag, bool & x_flag)
{
    int c;
    while ((c = getopt(argc, argv, "hd:e:o:t:ix")) != -1)
        switch(c){
            case 'h': cerr << "usage: " << argv[0]
                           << " [-h (for help)]"
                           << " [-d step_size_control_parameter]\n"
                           << "         [-e diagnostics_interval]"
                           << " [-o output_interval]\n"
                           << "         [-t total_duration]"
                           << " [-i (start output at t = 0)]\n"
                           << "         [-x (extra debugging diagnostics)]"
                           << endl;
                      return false;      
            case 'd': dt_param = atof(optarg);
                      break;
            case 'e': dt_dia = atof(optarg);
                      break;
            case 'i': i_flag = true;
                      break;
            case 'o': dt_out = atof(optarg);
                      break;
            case 't': dt_tot = atof(optarg);
                      break;
            case 'x': x_flag = true;
                      break;
            case '?': cerr << "usage: " << argv[0]
                           << " [-h (for help)]"
                           << " [-d step_size_control_parameter]\n"
                           << "         [-e diagnostics_interval]"
                           << " [-o output_interval]\n"
                           << "         [-t total_duration]"
                           << " [-i (start output at t = 0)]\n"
                           << "         [-x (extra debugging diagnostics)]"
                           << endl;
                      return false;     
            }

    return true;                  
}


void get_snapshot(real mass[], real pos[][NDIM], real vel[][NDIM],real v[],real p[], int n)
{
int h=0;
cout<<"if you want to use random generation,press 1: ";cin >> h;
if (h==1)
  {
srand(time(0));
ofstream m;
m.open("input.txt",ios::app);
double mrms;
double d=0;
for (int i = 0; i < n ; i++)
{
            mass[i]=rand()%200; 
            double a;
            for(int k=0;k<NDIM;k++)
             {
                a=(rand()%81)/80.0;
               	pos[i][k]=20*a-10;
		p[i]+=pow(pos[i][k],2);
	     }
		if(sqrt(p[i])>10) 
{
p[i]=0;
mass[i]=0;
i=i-1;		
}
		else
		{
		d+=pow(mass[i],2);
		for(int k=0;k<NDIM;k++)
                   {
                m<<"x"<<i+1<<k+1<<":\t"<<pos[i][k]<<"\t";
                   }  
		m<<endl;     
      	        for (int k = 0; k < NDIM; k++)
                  { 
                vel[i][k]=0;
                m<<"v"<<i+1<<k+1<<":\t"<<vel[i][k]<<"\t";
                  }
		m<<endl;
		}
}
mrms=sqrt(d/n);
for (int i = 0; i < n ; i++)
{
m<<i+1<<"\t"<<mass[i]<<endl;
}
m<<"mrms:"<<"\t"<<mrms<< endl;
m.close();


}


else
{
int f;
f=0;

double sumv[3];
cout<<"if you want to exert vcm=0, press 1: ";cin >> f;
if (f==1)
  {
for (int i = 0; i < n-1 ; i++)
{
       { cout<<"enter mass"<<i+1<<":";cin >> mass[i];  }  
               
        for (int k = 0; k < NDIM; k++)
           { cout<<"enter x"<<i+1<<k+1<<":";cin >> pos[i][k];}               
        for (int k = 0; k < NDIM; k++)
           {  cout<<"enter v"<<i+1<<k+1<<":";cin >> vel[i][k];
              sumv[k]+=vel[i][k]*mass[i]; 
            }
}
cout<<"enter mass"<<n<<":";cin >> mass[n-1];
 for (int k = 0; k < NDIM; k++)
           { cout<<"enter x"<<n<<k+1<<":";cin >> pos[n-1][k];}               
for(int k=0;k<3;k++)
{
vel[n-1][k]=-sumv[k]/mass[n-1];
}
  }


else 
{
double d;
double mrms;
ofstream m;
m.open("input.txt",ios::app);
  for (int i = 0; i < n ; i++)
   {
        cout<<"enter mass"<<i+1<<":";cin >> mass[i];
	d+=pow(mass[i],2);
	m<<"mass"<<i+1<<":\t"<<mass[i]<<endl;                     
        for (int k = 0; k < NDIM; k++)
           { 
	cout<<"enter x"<<i+1<<k+1<<":";cin >> pos[i][k];
	m<<"x"<<i+1<<k+1<<":\t"<<pos[i][k]<<"\t"; 
	   }
	m<<endl;              
        for (int k = 0; k < NDIM; k++)
           {  
	cout<<"enter v"<<i+1<<k+1<<":";cin >> vel[i][k];             
	m<<"v"<<i+1<<k+1<<":\t"<<vel[i][k]<<"\t";
	   }
	m<<endl;
    }
mrms=sqrt(d/n);
m<<"mrms:"<<"\t"<<mrms<< endl;
m.close();
}
}
}
    

void put_snapshot(const real mass[], const real pos[][NDIM],
                  const real vel[][NDIM], real v[],real p[],int n, real t)
{
    cout.precision(16);                 
double volume;
char name[]="property";
char end[]=".ods";
    ofstream o;                                                         
    for (int i = 0; i < n ; i++)
{                  
	char tp[10];
	sprintf(tp, "%d", i+1);
strcat(name,tp);
strcat(name,end);
o.open(name,ios::app); 
o<<t<<"\t"; 
        for (int k = 0; k < NDIM; k++)
          {              
            o<<pos[i][k]<<"\t"; 
	    p[i]+=pow(pos[i][k],2);                                       
           }
            o<<sqrt(p[i])<<"\t";
    	    
        for (int k= 0; k< NDIM; k++)
	{
            v[i]+=pow(vel[i][k],2);  
	}
 	    o<<sqrt(v[i])<<"\t"<<endl;
	    o.close();
 	    strcpy(name,"property");
}

ofstream d; 
d.open("volume.ods",ios::app);
d<<t<<"\t";
int c=0;
double l=0;
double rms;
for(int i=0;i<n;i++)
{

	if(sqrt(p[i])>50)
	 {
	p[i]=0;
	c+=1;
	 }
l+=p[i];
}
rms=sqrt(l/n);
volume=(4*M_PI/3)*pow(rms,3);
d<<volume<<" \t"<<c<<"\t";
for(int i=0;i<n;i++)
{
if(p[i]==0) 
d<<"mass"<<i+1<<"="<<mass[i]<<" is out of tidal radius";	
}
d<<endl;
d.close();



ofstream s;
s.open("v distribution.txt",ios::app);
double a,b;
double num[100];
double max=sqrt(v[1]);
for(int i=0;i<n;i++)
if(sqrt(v[i])>max) max=sqrt(v[i]);	
	for(int j=0;j<100;j++)
	{ num[j]=0;
a=(j*max/100);
b=((j+1)*max/100);
for(int i=0;i<n;i++)
{
if(a<=sqrt(v[i]) && sqrt(v[i])<=b) num[j]=num[j]+1;
}
 if (0<=t && t<2 || dt_tot<=t && t<=dt_tot+1 ||  ((dt_tot/2-1)<=t && t<=dt_tot/2)) {s<<j<<"\t"<<a<<"_"<<b<<"\t  "<<num[j]<<"\t"<<t<<endl;
}
         }
        s.close();
	           
}
    


void write_diagnostics(const real mass[], const real pos[][NDIM],
                       const real vel[][NDIM], const real acc[][NDIM],
                       const real jerk[][NDIM], int n, real t, real epot,
                       int nsteps, real & einit, bool init_flag,
                       bool x_flag)
{
    real ekin = 0;                  
    for (int i = 0; i < n ; i++)
        for (int k = 0; k < NDIM ; k++)
            ekin += 0.5 * mass[i] * vel[i][k] * vel[i][k];

    real etot = ekin + epot;          

    if (init_flag)                      
        einit = etot;                

    cerr << "at time t = " << t << " , after " << nsteps
         << " steps :\n  E_kin = " << ekin
         << " , E_pot = " << epot
         << " , E_tot = " << etot << endl;
    cerr << "                "
         << "absolute energy error: E_tot - E_init = "
         << etot - einit << endl;
    cerr << "                "
         << "relative energy error: (E_tot - E_init) / E_init = "
         << (etot - einit) / einit << endl;

    if (x_flag){
        cerr << "  for debugging purposes, here is the internal data "
             << "representation:\n";
        for (int i = 0; i < n ; i++){
            cerr << "    internal data for particle " << i+1 << " : " << endl;
            cerr << "      ";
            cerr << mass[i];
            for (int k = 0; k < NDIM; k++)
                cerr << ' ' << pos[i][k];
            for (int k = 0; k < NDIM; k++)
                cerr << ' ' << vel[i][k];
            for (int k = 0; k < NDIM; k++)
                cerr << ' ' << acc[i][k];
            for (int k = 0; k < NDIM; k++)
                cerr << ' ' << jerk[i][k];
            cerr << endl;
        }
    }
}
    
 

void evolve(const real mass[], real pos[][NDIM], real vel[][NDIM],real v[],real p[],int n, real & t, real dt_param, real dt_dia, real dt_out,
            real dt_tot, bool init_out, bool x_flag)
{
    cerr << "Starting a Hermite integration for a " << n
         << "-body system,\n  from time t = " << t 
         << " with time step control parameter dt_param = " << dt_param
         << "  until time " << t + dt_tot 
         << " ,\n  with diagnostics output interval dt_dia = "
         << dt_dia << ",\n  and snapshot output interval dt_out = "
         << dt_out << "." << endl;

    real (* acc)[NDIM] = new real[n][NDIM];          
    real (* jerk)[NDIM] = new real[n][NDIM];        
    real epot;                       
    real coll_time;                

    get_acc_jerk_pot_coll(mass, pos, vel, acc, jerk, n, epot, coll_time);

    int nsteps = 0;              
    real einit;               

    write_diagnostics(mass, pos, vel, acc, jerk, n, t, epot, nsteps, einit,
                      true, x_flag);
    if (init_out)                                    
        put_snapshot(mass, pos, vel, v,p, n, t);

    real t_dia = t + dt_dia;          
    real t_out = t + dt_out;          
    real t_end = t + dt_tot;       

    while (true){
        while (t < t_dia && t < t_out && t < t_end){
            real dt = dt_param * coll_time;
            evolve_step(mass, pos, vel, acc, jerk, n, t, dt, epot, coll_time);
            nsteps++;
        }
        if (t >= t_dia){
            write_diagnostics(mass, pos, vel, acc, jerk, n, t, epot, nsteps,
                              einit, false, x_flag);
            t_dia += dt_dia;
        }
        if (t >= t_out){
            put_snapshot(mass, pos, vel, v,p, n, t);
            t_out += dt_out;
        }
        if (t >= t_end)
            break;
    }

    delete[] acc;
    delete[] jerk;
}


void evolve_step(const real mass[], real pos[][NDIM], real vel[][NDIM],
                 real acc[][NDIM], real jerk[][NDIM], int n, real & t,
                 real dt, real & epot, real & coll_time)
{
    real (* old_pos)[NDIM] = new real[n][NDIM];
    real (* old_vel)[NDIM] = new real[n][NDIM];
    real (* old_acc)[NDIM] = new real[n][NDIM];
    real (* old_jerk)[NDIM] = new real[n][NDIM];

    for (int i = 0; i < n ; i++)
        for (int k = 0; k < NDIM ; k++){
          old_pos[i][k] = pos[i][k];
          old_vel[i][k] = vel[i][k];
          old_acc[i][k] = acc[i][k];
          old_jerk[i][k] = jerk[i][k];
        }

    predict_step(pos, vel, acc, jerk, n, dt);
    get_acc_jerk_pot_coll(mass, pos, vel, acc, jerk, n, epot, coll_time);
    correct_step(pos, vel, acc, jerk, old_pos, old_vel, old_acc, old_jerk,
                 n, dt);
    t += dt;

    delete[] old_pos;
    delete[] old_vel;
    delete[] old_acc;
    delete[] old_jerk;
}



void predict_step(real pos[][NDIM], real vel[][NDIM], 
                  const real acc[][NDIM], const real jerk[][NDIM],
                  int n, real dt)
{
    for (int i = 0; i < n ; i++)
        for (int k = 0; k < NDIM ; k++){
            pos[i][k] += vel[i][k]*dt + acc[i][k]*dt*dt/2
                                      + jerk[i][k]*dt*dt*dt/6;
            vel[i][k] += acc[i][k]*dt + jerk[i][k]*dt*dt/2;
        }
}


void correct_step(real pos[][NDIM], real vel[][NDIM], 
                  const real acc[][NDIM], const real jerk[][NDIM],
                  const real old_pos[][NDIM], const real old_vel[][NDIM], 
                  const real old_acc[][NDIM], const real old_jerk[][NDIM],
                  int n, real dt)
{
    for (int i = 0; i < n ; i++)
        for (int k = 0; k < NDIM ; k++){
            vel[i][k] = old_vel[i][k] + (old_acc[i][k] + acc[i][k])*dt/2
                                      + (old_jerk[i][k] - jerk[i][k])*dt*dt/12;
            pos[i][k] = old_pos[i][k] + (old_vel[i][k] + vel[i][k])*dt/2
                                      + (old_acc[i][k] - acc[i][k])*dt*dt/12;
        }
}



void get_acc_jerk_pot_coll(const real mass[], const real pos[][NDIM],
                           const real vel[][NDIM], real acc[][NDIM],
                           real jerk[][NDIM], int n, real & epot,
                           real & coll_time)
{
    for (int i = 0; i < n ; i++)
        for (int k = 0; k < NDIM ; k++)
            acc[i][k] = jerk[i][k] = 0;
    epot = 0;
    const real VERY_LARGE_NUMBER = 1e300;
    real coll_time_q = VERY_LARGE_NUMBER;    
    real coll_est_q;                           
                                              
    for (int i = 0; i < n ; i++){
        for (int j = i+1; j < n ; j++){            
            real rji[NDIM];                       
            real vji[NDIM];                      
            for (int k = 0; k < NDIM ; k++){
                rji[k] = pos[j][k] - pos[i][k];
                vji[k] = vel[j][k] - vel[i][k];
            }
            real r2 = 0;                           
            real v2 = 0;                           
            real rv_r2 = 0;                        
            for (int k = 0; k < NDIM ; k++){
                r2 += rji[k] * rji[k];
                v2 += vji[k] * vji[k];
                rv_r2 += rji[k] * vji[k];
            }
            rv_r2 /= r2;
            real r = sqrt(r2);                     
            real r3 = r * r2;                      



            epot -= mass[i] * mass[j] / r;



            real da[3];                            
            real dj[3];                           
            for (int k = 0; k < NDIM ; k++){
                da[k] = rji[k] / r3;                           
                dj[k] = (vji[k] - 3 * rv_r2 * rji[k]) / r3;  
            }
            for (int k = 0; k < NDIM ; k++){
                acc[i][k] += mass[j] * da[k];                 
                acc[j][k] -= mass[i] * da[k];                 
                jerk[i][k] += mass[j] * dj[k];               
                jerk[j][k] -= mass[i] * dj[k];              
            }



            coll_est_q = (r2*r2) / (v2*v2);
            if (coll_time_q > coll_est_q)
                coll_time_q = coll_est_q;



            real da2 = 0;                                   
            for (int k = 0; k < NDIM ; k++)                
                da2 += da[k] * da[k];                     
            double mij = mass[i] + mass[j];              
            da2 *= mij * mij;                             

            coll_est_q = r2/da2;
            if (coll_time_q > coll_est_q)
                coll_time_q = coll_est_q;
        }                                     
    }                                              
    coll_time = sqrt(sqrt(coll_time_q));            
}                                             


