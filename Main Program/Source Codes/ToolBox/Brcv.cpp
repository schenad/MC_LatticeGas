#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <ctime>
#include <sys/time.h>
double coverage=71.95;
#define A (4e-4)
#define be (0.00225)  //Bond energy
const unsigned long long int steps=100000001ULL;      //Total MC steps
#define trials (1)
#define interval (1)
#define temps (1)
#define S (100)  //Size of the plane
#define bl (1.)     //Size of 1 unit of plane (Bond length)
#define kb (8.6173324e-5)
char firstname[10]="Br-";

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
  int atcount(int sample[][S]);
  int ran(int);
  int xo,yo,xp,yp,trial;
  double Ei,Ef,rn,T,Te,Etot;
  char Astr[6];
  char Tstr[6];
  char cover[6];
  double Ekp[1000];
  cout<<"Please input the starting temperature: ";
  cin>>Te;
  cout<<"Please input the starting trial: ";
  cin>>trial;
  timeval t_start, t_end;
  unsigned long long remaintime=0ULL;
  const int n=S*S*coverage/100;
  char filename[40];
  if(Te==(int)Te)
    sprintf(Tstr,"%d",(int)(Te));
  else
    sprintf(Tstr,"%.1f",Te);

  cout<<"Processing........Please wait......"<<endl;
  gettimeofday(&t_start,NULL);

  ofstream recordE;
  ofstream recordE1;
  ofstream recordE2;
  ofstream mcfile;

  for(int tr=0; tr<trials; tr++){
    T=Te;
    for(int temp=0;temp<temps;temp++,T+=interval){

  if(T==(int)T)
    sprintf(Astr,"%d",(int)(T));
  else
    sprintf(Astr,"%.1f",T);
  sprintf(cover,"%.2f",coverage);
  for(int i=0;i<5;i++){
    if(cover[i]=='.')
      cover[i]='_';
     if(Astr[i]=='.')
      Astr[i]='_';
     if(Tstr[i]=='.')
      Tstr[i]='_';
  }
  strcpy(filename,firstname);
  strcat(strcat(strcat(filename,cover),"-"),Astr);
  strcat(filename,".txt");

  char num[50];
  sprintf(num,"%d",tr+trial);
  strcat(strcat(strcat(strcat(num,"-"),firstname),cover),"-");
  strcat(num,Tstr);
  cout<<endl<<endl<<"Results are to be output to "<<num<<endl;

  for(int j=0;j<S;j++)
    for(int i=0;i<S;i++)
      sample[j][i]=0;           //Initialize the plane.

  init_genrand((unsigned)time(NULL));    //Changing seed.

  for(int k=0;k<n;){
    xo=ran(S);
    yo=ran(S);
    if(sample[yo][xo]==0){       //Randomly distribute the atoms on the plane.
      sample[yo][xo]=1;
      k++;
    }
  }

  for(int k=0;k<S*S;k++)
    if(*(&sample[0][0]+k)==1)
      Etot=Energy(sample,k%S,k/S)/2;

  for(unsigned int s=0; s<steps;s++){

   for(xo=ran(S),yo=ran(S);sample[yo][xo]==0;){     //Randomly choose an atom
     xo=ran(S);
     yo=ran(S);
  }

  for(xp=ran(S),yp=ran(S);sample[yp][xp]==1;){ //Randomly choose an space
     xp=ran(S);
     yp=ran(S);
  }

  Ei=Energy(sample,xo,yo);   //Initial Energy

  sample[yo][xo]=0;
  sample[yp][xp]=1;              //Do the swaping

  Ef=Energy(sample,xp,yp);
  Etot+=Ef-Ei;

  if(Ef>=Ei){
    rn=ran(LOWER_MASK)/(LOWER_MASK+0.);
    if(rn>=exp(-(Ef-Ei)/kb/T)){
      sample[yo][xo]=1;
      sample[yp][xp]=0;              //Cancel the swaping
      Etot-=Ef-Ei;
    }
  }

  cout<<fixed;
  cout.precision(2);

  if(s%(steps/10000)==0){
    gettimeofday(&t_end,NULL);
    remaintime=(temps*(trials-tr)-temp-s/steps)*((1000000*(t_end.tv_sec-t_start.tv_sec)+t_end.tv_usec-t_start.tv_usec)/100.);
    cout<<"\r"<<(((double)s/steps+temp+temps*tr)*100.)/trials/temps<<'%'<<"  Time remaining: "<<remaintime/3600<<'h'<<(remaintime/60-60*(remaintime/3600))<<"min"<<remaintime%60<<"s ";
    cout<<cout.precision(4)<<fixed<<'\b'<<"dEp="<<(Ef-Ei)*1000<<"meV "<<"Ep="<<Etot<<"eV     ";
    cout.flush();
    gettimeofday(&t_start,NULL);
  }

  if(s!=0&&s%(steps/1000)==0)
    Ekp[s/(steps/1000)-1]=Etot+atcount(sample)*kb*T;

  }
  mcfile.open(filename,ios::app);
  for(int j=0;j<S;j++){            //Output results
     for(int i=0;i<S;i++)
      mcfile<<sample[j][i]<<" ";
     mcfile<<endl;
  }
  mcfile<<endl;
  mcfile.close();

  double Esum=0;
  double Esum2=0;
  char numE[40];
  char numE2[40];
  char allE[40];
  strcpy(numE,num);
  strcpy(numE2,num);
  strcpy(allE,filename);
  strcat(strcat(numE,"E"),".txt");
  strcat(strcat(numE2,"E2"),".txt");
  strcat(strcat(allE,"E"),".txt");
  recordE.precision(18);
  recordE1.precision(18);
  recordE2.precision(18);
  recordE.open(allE,ios::app);
  recordE1.open(numE,ios::app);
  recordE2.open(numE2,ios::app);
  for(int i=900;i<1000;i++){
    Esum+=Ekp[i];
    Esum2+=Ekp[i]*Ekp[i];
    recordE<<Ekp[i]<<endl;
  }
  recordE1<<Esum/100.<<endl;
  recordE2<<Esum2/100.<<endl;
  recordE.close();
  recordE1.close();
  recordE2.close();

    }
  }

  cout<<endl<<endl<<"Finished! "<<endl<<endl;

  return 0;
}

double BondEnergy(int s[][S],int i,int j){

  int bn=0;

  if(s[j+1-S*(j>=S-1)][i]==1)
    bn++;
  if(s[j+1-S*(j>=S-1)][i-1+S*(i<=0)]==1)
    bn++;
  if(s[j][i-1+S*(i<=0)]==1)
    bn++;
  if(s[j-1+S*(j<=0)][i]==1)
    bn++;
  if(s[j-1+S*(j<=0)][i+1-S*(i>=S-1)]==1)
    bn++;
  if(s[j][i+1-S*(i>=S-1)]==1)
    bn++;

  return -bn*be;
}

int atcount(int sample[][S]){

  int numsg=0;

    for(int pt=0;pt<S*S;pt++){
      if(*(&sample[0][0]+pt)==1)
	if(BondEnergy(sample,pt%S,pt/S)==0)
    		 numsg++;
    }

    return numsg;

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
