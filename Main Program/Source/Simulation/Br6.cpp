#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <ctime>
#include <sys/time.h>
//#include <unistd.h>
double coverage=71.95;
double A=4e-4;
#define be (0.00225)  //Bond energy
#define steps (5000001)      //Total MC steps
#define S (200)  //Size of the plane
#define bl (1.)     //Size of 1 unit of plane(Bond length)
#define kb (8.6173324e-5)
#define trials (1)
char firstname[10]="Br-";
#define T (5.)
#define N 624
#define M 397
#define MATRIX_A 0x9908b0dfUL   /* constant vector a */
#define UPPER_MASK 0x80000000UL /* most significant w-r bits */
#define LOWER_MASK 0x7fffffffUL /* least significant r bits */
static unsigned long mt[N]; /* the array for the state vector  */
static int mti=N+1; /* mti==N+1 means mt[N] is not initialized */
void init_genrand(unsigned long s)
{
    mt[0]= s & 0xffffffffUL;
    for (mti=1; mti<N; mti++) {
        mt[mti] =
	    (1812433253UL * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti);
        mt[mti] &= 0xffffffffUL;
    }
}
unsigned long genrand_int32(void)
{
    unsigned long y;
    static unsigned long mag01[2]={0x0UL, MATRIX_A};
    if (mti >= N) { /* generate N words at one time */
        int kk;
        if (mti == N+1)   /* if init_genrand() has not been called, */
            init_genrand(5489UL); /* a default initial seed is used */
        for (kk=0;kk<N-M;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+M] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        for (;kk<N-1;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+(M-N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        y = (mt[N-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
        mt[N-1] = mt[M-1] ^ (y >> 1) ^ mag01[y & 0x1UL];
        mti = 0;
    }
    y = mt[mti++];
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);
    return y;
}
int ran(int x){
  return (genrand_int32()%x);
}
using namespace std;

int main(){

  int sample[S][S];
  double Energy(int sample[][S],int,int);
  int ran(int);
  int xo,yo,xp,yp,index,step;
  long remaintime;
  double Ei,Ef,rn;
  timeval t_start, t_end;
  char imanum[3];
  char Astr[6];
  char cover[6];
  char input;
  //cout<<"Please input the prefix of the output file:"<<endl;
  //cin>>firstname;
  cout<<"Please input the starting index of the trial: ";
  cin>>index;
  //cout<<"Please input the coverage: ";
  //cin>>coverage;
  const int n=S*S*coverage/100;   //Number of atoms
  //cout<<"Please input the repulsion constant A: ";
  //cin>>A;

  for(int trial=0;trial<trials;trial++){

    step=0;
  char filename[30];
  sprintf(imanum,"%d",(index+trial));
  cout<<"The index of the current trial is "<<imanum<<endl;
  sprintf(Astr,"%.1f",be*1000);
  sprintf(cover,"%.2f",coverage);
  strcpy(filename,firstname);
  strcat(strcat(strcat(filename,cover),"-"),Astr);
  strcat(strcat(strcat(filename,"-"),imanum),".txt");

  fstream check;
  check.open(filename,ios::in);
  if(!check)
    cout<<"The data is to be written into "<<filename<<'.'<<endl;
  else{
    cout<<"The data is to be written into "<<filename<<", which is already existed."<<endl<<"Input c to continue, input o to overwrite, or input other char to exit: ";
    cin>>input;
    if(input=='c'){
      ifstream fin(filename);
      fin>>step;
      if(step==0||step==1){
      cout<<endl<<"The file you choose may be wrong. The program is exiting now. No change is made to the file."<<endl;
      return 1;
      }
    }
    else if(input=='o');
    else
      return 1;
  }

  cout<<"Processing........Please wait......"<<endl<<endl;

  for(int j=0;j<S;j++)
    for(int i=0;i<S;i++)
      sample[j][i]=0;           //Initialize the plane.

  init_genrand((unsigned)time(NULL));    //Changing seed.

  for(int k=0;k<n;){
    xo=ran(S-1);
    yo=ran(S-1);
    if(sample[yo][xo]==0){       //Randomly distribute the atoms on the plane.
      sample[yo][xo]=1;
      k++;
    }
  }

  for(int s=step; s<steps;s++){

    if(s==0)
      gettimeofday(&t_start,NULL);

  for(xo=ran(S-1),yo=ran(S-1);sample[yo][xo]==0;){     //Randomly choose an atom
     xo=ran(S-1);
     yo=ran(S-1);
  }

  for(xp=ran(S-1),yp=ran(S-1);sample[yp][xp]==1;){ //Randomly choose an space
     xp=ran(S-1);
     yp=ran(S-1);
  }

  Ei=Energy(sample,xo,yo);   //Initial Energy

  sample[yo][xo]=0;
  sample[yp][xp]=1;              //Do the swaping

  Ef=Energy(sample,xp,yp);

  if(Ef>=Ei){
    rn=ran(LOWER_MASK)/(LOWER_MASK+0.);
    if(rn>=exp(-(Ef-Ei)/kb/T)){
      sample[yo][xo]=1;
      sample[yp][xp]=0;              //Cancel the swaping
    }
  }

  cout<<fixed;
  cout.precision(2);

  if(s%(steps/10000)==0){
    gettimeofday(&t_end,NULL);
    remaintime=(steps*(trials-trial)-s)/((steps/10000)/((1000000*(t_end.tv_sec-t_start.tv_sec)+t_end.tv_usec-t_start.tv_usec)/1000000.));
    cout<<"\r"<<(s*100.)/steps<<'%'<<"  Time remaining: "<<remaintime/3600<<'h'<<(remaintime/60-60*(remaintime/3600))<<"min"<<remaintime%60<<"s    ";
    cout.flush();
    gettimeofday(&t_start,NULL);
  }

  if(s%10000==0){
    if((s/10000)!=(steps/10000)){
      ofstream mcfile1;
      mcfile1.open(filename);
      mcfile1<<s<<endl;
      for(int j=0;j<S;j++){            //Output results
         for(int i=0;i<S;i++)
          mcfile1<<sample[j][i]<<" ";
         mcfile1<<endl;
      }
    mcfile1.close();
    }
    else{
      ofstream mcfile2;
      mcfile2.open(filename);
      for(int j=0;j<S;j++){            //Output results
         for(int i=0;i<S;i++)
           mcfile2<<sample[j][i]<<" ";
         mcfile2<<endl;
      }
      mcfile2<<endl;
      mcfile2.close();
    }
  }
  }
  cout<<endl<<endl<<"Finished! The result has been output to "<<filename<<endl<<endl;
  }
  return 0;
}

double BondEnergy(int s[][S],int i,int j){

  int bn=0;

  if(i==0){
    if(j==0){
      if(s[0][1]==1)
	bn++;
      if(s[1][0]==1)
	bn++;
      if(s[0][S-1]==1)
	bn++;
      if(s[S-1][0]==1)
	bn++;
      if(s[1][S-1]==1)
	bn++;
      if(s[S-1][1]==1)
	bn++;
    }
    else if(j==S-1){
      if(s[0][0]==1)
	bn++;
      if(s[0][S-1]==1)
	bn++;
      if(s[S-1][S-1]==1)
	bn++;
      if(s[S-2][0]==1)
	bn++;
      if(s[S-2][1]==1)
	bn++;
      if(s[S-1][1]==1)
	bn++;
    }
    else{
      if(s[j+1][0]==1)
	bn++;
      if(s[j+1][S-1]==1)
	bn++;
      if(s[j][S-1]==1)
	bn++;
      if(s[j-1][0]==1)
	bn++;
      if(s[j-1][1]==1)
	bn++;
      if(s[j][1]==1)
	bn++;
    }
  }

  else if(i==S-1){
    if(j==0){
      if(s[1][S-1]==1)
	bn++;
      if(s[1][S-2]==1)
	bn++;
      if(s[0][S-2]==1)
	bn++;
      if(s[S-1][S-1]==1)
	bn++;
      if(s[S-1][0]==1)
	bn++;
      if(s[0][0]==1)
	bn++;
    }
    else if(j==S-1){
      if(s[0][S-1]==1)
	bn++;
      if(s[0][S-2]==1)
	bn++;
      if(s[S-1][S-2]==1)
	bn++;
      if(s[S-2][S-1]==1)
	bn++;
      if(s[S-2][0]==1)
	bn++;
      if(s[S-1][0]==1)
	bn++;
    }
    else{
      if(s[j+1][S-1]==1)
	bn++;
      if(s[j+1][S-2]==1)
	bn++;
      if(s[j][S-2]==1)
	bn++;
      if(s[j-1][S-1]==1)
	bn++;
      if(s[j-1][0]==1)
	bn++;
      if(s[j][0]==1)
	bn++;
    }
  }

  else{
    if(j==0){
      if(s[1][i]==1)
	bn++;
      if(s[1][i-1]==1)
	bn++;
      if(s[0][i-1]==1)
	bn++;
      if(s[S-1][i]==1)
	bn++;
      if(s[S-1][i+1]==1)
	bn++;
      if(s[0][i+1]==1)
	bn++;
    }
    else if(j==S-1){
      if(s[0][i]==1)
	bn++;
      if(s[0][i-1]==1)
	bn++;
      if(s[S-1][i-1]==1)
	bn++;
      if(s[S-2][i]==1)
	bn++;
      if(s[S-2][i+1]==1)
	bn++;
      if(s[S-1][i+1]==1)
	bn++;
    }
    else{
      if(s[j+1][i]==1)
	bn++;
      if(s[j+1][i-1]==1)
	bn++;
      if(s[j][i-1]==1)
	bn++;
      if(s[j-1][i]==1)
	bn++;
      if(s[j-1][i+1]==1)
	bn++;
      if(s[j][i+1]==1)
	bn++;
    }
  }
  return -bn*be;
}

double re(int a, int b, int c, int d){

  int x=a-c;
  int y=b-d;
  double sr=0.;

  if(x==0&&y==0)
    return sr;
  else{
    if(x>-2*y && 2*x<-y && x+S<y)
      sr=(x+S)*(x+S)+(y-S)*(y-S)+(x+S)*(y-S);
    else if(2*x>=-y && x+2*y>S && x<=y && 2*x+y<=2*S)
      sr=x*x+(y-S)*(y-S)+x*(y-S);
    else if(2*x+y>2*S && x+2*y>2*S)
      sr=(x-S)*(x-S)+(y-S)*(y-S)+(x-S)*(y-S);
    else if(x<=-2*y && x+2*y>=-2*S && x<=y && 2*x+y<-S)
      sr=(x+S)*(x+S)+y*y+(x+S)*y;
    else if(x>y && x+2*y<=2*S && 2*x+y>S && x>=-2*y)
      sr=(x-S)*(x-S)+y*y+(x-S)*y;
    else if(x+2*y<-2*S && 2*x+y<-2*S)
      sr=(x+S)*(x+S)+(y+S)*(y+S)+(x+S)*(y+S);
    else if(2*x+y>=-2*S && x>y && x+2*y<-S && 2*x<=-y)
      sr=x*x+(y+S)*(y+S)+x*(y+S);
    else if(2*x>-y && x-y>S && x<-2*y)
      sr=(x-S)*(x-S)+(y+S)*(y+S)+(x-S)*(y+S);
    else
      sr=x*x+y*y+x*y;
  }
  return A*pow(sr,-3/2)/pow(bl,3);
}

double RepulsiveEnergy(int sample[][S],int i,int j){

  double sum=0.;

  for(int k=0;k<S*S;k++)
    if(*(&sample[0][0]+k)==1)
      sum+=re(i,j,k%S,k/S);

  return sum;
}

double Energy(int sample[][S],int x,int y){
  return BondEnergy(sample,x,y)+RepulsiveEnergy(sample,x,y);
}
