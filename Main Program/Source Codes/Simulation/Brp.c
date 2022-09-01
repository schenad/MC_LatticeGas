#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <ctime>
#define steps (10000000)      //Total MC steps
#define interval (100000)
#define bl (1.)     //Size of 1 unit of plane(Bond length)
#define S (200)  //Size of the plane
#define n (13460)   //Number of atoms
#define be (0.002)  //Bond energy
#define A (be/2.5)     //Constant of the repulsive energy.
#define kb (8.6173324e-5)
#define T (5.)
#define N 624
#define M 397
#define MATRIX_A 0x9908b0dfUL   /* constant vector a */
#define UPPER_MASK 0x80000000UL /* most significant w-r bits */
#define LOWER_MASK 0x7fffffffUL /* least significant r bits */

static unsigned long mt[N]; /* the array for the state vector  */
static int mti=N+1; /* mti==N+1 means mt[N] is not initialized */

/* initializes mt[N] with a seed */
void init_genrand(unsigned long s)
{
    mt[0]= s & 0xffffffffUL;
    for (mti=1; mti<N; mti++) {
        mt[mti] =
	    (1812433253UL * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti);
        /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
        /* In the previous versions, MSBs of the seed affect   */
        /* only MSBs of the array mt[].                        */
        /* 2002/01/09 modified by Makoto Matsumoto             */
        mt[mti] &= 0xffffffffUL;
        /* for >32 bit machines */
    }
}

/* initialize by an array with array-length */
/* init_key is the array for initializing keys */
/* key_length is its length */
/* slight change for C++, 2004/2/26 */
void init_by_array(unsigned long init_key[], int key_length)
{
    int i, j, k;
    init_genrand(19650218UL);
    i=1; j=0;
    k = (N>key_length ? N : key_length);
    for (; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1664525UL))
          + init_key[j] + j; /* non linear */
        mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++; j++;
        if (i>=N) { mt[0] = mt[N-1]; i=1; }
        if (j>=key_length) j=0;
    }
    for (k=N-1; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1566083941UL))
          - i; /* non linear */
        mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++;
        if (i>=N) { mt[0] = mt[N-1]; i=1; }
    }

    mt[0] = 0x80000000UL; /* MSB is 1; assuring non-zero initial array */
}

/* generates a random number on [0,0xffffffff]-interval */
unsigned long genrand_int32(void)
{
    unsigned long y;
    static unsigned long mag01[2]={0x0UL, MATRIX_A};
    /* mag01[x] = x * MATRIX_A  for x=0,1 */

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

    /* Tempering */
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
  double Be(int s[][S],int,int);
  int ran(int);
  int xo,yo,xp,yp;
  double Ei,Ef,rn;

  ofstream mcfile;
  mcfile.open("Brp.txt");

  for(int j=0;j<S;j++)
    for(int i=0;i<S;i++)
      sample[j][i]=0;           //Initialize the plane.

  for(int k=0;k<n;){
    xo=ran(S-1);
    yo=ran(S-1);
    if(sample[yo][xo]==0){       //Randomly distribute the atoms on the plane.
      sample[yo][xo]=1;
      k++;
    }
  }

  for(int s=0; s<steps;s++){

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

  if(Ef==Ei){
    if(ran(1)==0){
      sample[yo][xo]=1;
      sample[yp][xp]=0;              //Cancel the swaping
    }
  }
  else if(Ef>Ei){
    rn=ran(LOWER_MASK)/(LOWER_MASK+0.);
    if(rn>exp(-(Ef-Ei)/kb/T)){
      sample[yo][xo]=1;
      sample[yp][xp]=0;              //Cancel the swaping
    }
  }
    if(s%interval==0){
	  for(int j=0;j<S;j++){            //Output results
	     for(int i=0;i<S;i++)
	      mcfile<<sample[j][i]<<" ";
	     mcfile<<endl;
	     }
	  mcfile<<endl;
   }
  }
  mcfile.close();
  return 0;
}

double Be(int s[][S],int i,int j){

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

  return bn;

}

double BondEnergy(int s[][S],int i,int j){

  return -(Be(s,i,j)*be);

}

double re(int a, int b, int c, int d){

  int x=a-c;
  int y=b-d;
  double s=0.;

  if(x==0&&y==0)
    return s;

  else{
    if(x>-2*y && 2*x<-y && x+S<y)
      s=(x+S)*(x+S)+(y-S)*(y-S)+(x+S)*(y-S);
    else if(2*x>=-y && x+2*y>S && x<=y && 2*x+y<=2*S)
      s=x*x+(y-S)*(y-S)+x*(y-S);
    else if(2*x+y>2*S && x+2*y>2*S)
      s=(x-S)*(x-S)+(y-S)*(y-S)+(x-S)*(y-S);
    else if(x<=-2*y && x+2*y>=-2*S && x<=y && 2*x+y<-S)
      s=(x+S)*(x+S)+y*y+(x+S)*y;
    else if(x>y && x+2*y<=2*S && 2*x+y>S && x>=-2*y)
      s=(x-S)*(x-S)+y*y+(x-S)*y;
    else if(x+2*y<-2*S && 2*x+y<-2*S)
      s=(x+S)*(x+S)+(y+S)*(y+S)+(x+S)*(y+S);
    else if(2*x+y>=-2*S && x>y && x+2*y<-S && 2*x<=-y)
      s=x*x+(y+S)*(y+S)+x*(y+S);
    else if(2*x>-y && x-y>S && x<-2*y)
      s=(x-S)*(x-S)+(y+S)*(y+S)+(x-S)*(y+S);
    else
      s=x*x+y*y+x*y;
  }

  return A*pow(s,-3/2)/pow(bl,3);

}

double RepulsiveEnergy(int s[][S],int i,int j){

  double sum=0.;

  for(int k=0;k<S*S;k++)
    if(*(&s[0][0]+k)==1)
      sum+=re(i,j,k%S,k/S);

  return sum;
}

double Energy(int sample[][S],int x,int y){
  return BondEnergy(sample,x,y)+RepulsiveEnergy(sample,x,y);
}
