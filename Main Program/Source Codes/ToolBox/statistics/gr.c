#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#define S (200)  //Size of the plane
#define trials (10)
#define bl (0.44) // Bond length
#define range (30.)
#define bars (150)
#define interval (range/bars)

using namespace std;

const int n=S*S*trials;    //Number of atoms approximately(can be more than the exact number, but cannot be less)
int size[n];
int n_i=0;
int ip[trials][n];
int np;
int trial=0;
double A;

double F(double A, double d){
  return A*d*(1-d*(4*S*bl-4*d+M_PI*d)/M_PI/S/bl/S/bl);
}

int main(){

  int samples[trials][S][S];
  int sample[S][S];
  int getsize(int s[][S],int);
  int check(int);
  double dis(int, int);
  int isboundary(int s[][S],int);
  int isl=0;
  int isla=0;
  int isno;
  double distance=S*bl;
  double r[trials][n];
  int counts[trials][bars];
  double g[trials][bars];
  double G[bars];

  ifstream fin("Br.txt");
  for(int k=0;k<trials;k++)
    for(int i=0;i<S;i++)
       for(int j=0;j<S;j++)
         fin>>samples[k][i][j];
  fin.close();

  for(int t=0;t<trials;t++,trial++){

   np=0;
   isla=0;
   A=0.;
   isno=n_i;
  for(int i=0;i<bars;i++)
    counts[trial][i]=0;

   for(int i=0;i<S;i++)
       for(int j=0;j<S;j++) 
         sample[i][j]=samples[trial][i][j];

   for(int pt=0;pt<S*S;pt++){
     if(*(&sample[0][0]+pt)==1)
       if(check(pt)==1){
 	size[n_i]++;
        getsize(sample,pt);
        n_i++;
      }
   }
   
  int islp=0;
  int dislp;
  int pi[trials][n_i][1000];
  for(;isl<n_i;isl++,isla++){
    for(dislp=0;dislp<size[isl];dislp++,islp++);
    for(int i=0;i<dislp;i++)
       pi[trial][isla][i]=ip[trial][islp-dislp+i];
  }

  for(int j=0; j<isla-1;j++)
    for(int i=1;i<isla;i++)
      if(i>j)
	for(int l=0;l<size[j+isno];l++)
	  for(int k=0;k<size[i+isno];k++)
	    if(isboundary(sample,pi[trial][j][l])==1 || isboundary(sample,pi[trial][i][k])==1)
	      if(dis(pi[trial][j][l],pi[trial][i][k])<distance || dis(pi[trial][j][l],pi[trial][i][k])!=0){
	        distance=dis(pi[trial][j][l],pi[trial][i][k]);
	        r[trial][(2*isla-j-1)*j/2+i-j-1]=distance;
	    }
 
  for(int i=0;i<bars;i++)
    for(int j=0;j<isla*(isla-1)/2;j++)
      if(r[trial][j]>i*interval && r[trial][j]<=(i+1.)*interval)
	counts[trial][i]++;

  double idiff=1000.1;
  double diff=1000.;
  for(;idiff>diff;A=A+0.001){
    idiff=diff;
    for(int i=0;i<bars;i++)
      diff+=abs(F(A,i*interval)-counts[trial][i]);
    diff=diff/bars;
  }
  cout<<A<<endl;
  g[trial][0]=0;
  for(int i=1;i<bars;i++)
    g[trial][i]=counts[trial][i]/F(A,i*interval);
  }

  for(int i=0;i<bars;i++)
    G[i]=0.;

  for(int i=0;i<bars;i++)
    for(int j=0;j<trials;j++)
      G[i]+=g[j][i];

  ofstream mcfile;

  mcfile.open("G_r.dat"); 
  for(int i=1;i<bars;i++)
    mcfile<<i*interval<<' '<<G[i]/trials<<endl;
  mcfile.close();

  mcfile.open("r.dat"); 
  for(int i=0;i<(n_i-isno)*(n_i-isno-1)/2;i++)
    mcfile<<r[0][i]<<endl;
  mcfile.close();

  return 0;
}

int check(int w){
  for(int i=0;i<np;i++){
    if(w==ip[trial][i])
	return 0;
  }
  np++;
  ip[trial][np-1]=w;
  return 1;
}

int getatom(int s[][S],int u){

  int inp=np;

  if(u%S==0){
    if(u/S==0){
      if(s[0][1]==1)
	    check(1);
      if(s[1][0]==1)
	    check(S);
      if(s[0][S-1]==1)
	    check(S-1);
      if(s[S-1][0]==1)
	    check(S*(S-1));
      if(s[1][S-1]==1)
	    check(S+S-1);
      if(s[S-1][1]==1)
	    check(S*(S-1)+1);
    }
    else if(u/S==S-1){
      if(s[0][0]==1)
	    check(0);
      if(s[0][S-1]==1)
	    check(S-1);
      if(s[S-1][S-1]==1)
	    check(S*S-1);
      if(s[S-2][0]==1)
	    check(S*(S-2));
      if(s[S-2][1]==1)
	    check(S*(S-2)+1);
      if(s[S-1][1]==1)
	    check(S*(S-1)+1);
    }
    else{
      if(s[u/S+1][0]==1)
	    check(u+S);
      if(s[u/S+1][S-1]==1)
	    check(u+S+S-1);
      if(s[u/S][S-1]==1)
	    check(u+S-1);
      if(s[u/S-1][0]==1)
	    check(u-S);
      if(s[u/S-1][1]==1)
	    check(u-S+1);
      if(s[u/S][1]==1)
	    check(u+1);
    }
  }

  else if(u%S==S-1){
    if(u/S==0){
      if(s[1][S-1]==1)
	    check(S+S-1);
      if(s[1][S-2]==1)
	    check(S+S-2);
      if(s[0][S-2]==1)
	    check(S-2);
      if(s[S-1][S-1]==1)
	    check(S*S-1);
      if(s[S-1][0]==1)
	    check(S*(S-1));
      if(s[0][0]==1)
	    check(0);
    }
    else if(u/S==S-1){
      if(s[0][S-1]==1)
	    check(S-1);
      if(s[0][S-2]==1)
	    check(S-2);
      if(s[S-1][S-2]==1)
	    check(S*S-2);
      if(s[S-2][S-1]==1)
	    check(S*S-S-1);
      if(s[S-2][0]==1)
	    check(S*(S-2));
      if(s[S-1][0]==1)
	    check(S*(S-1));
    }
    else{
      if(s[u/S+1][S-1]==1)
	    check(u+S);
      if(s[u/S+1][S-2]==1)
	    check(u+S-1);
      if(s[u/S][S-2]==1)
	    check(u-1);
      if(s[u/S-1][S-1]==1)
	    check(u-S);
      if(s[u/S-1][0]==1)
	    check(u-S-S+1);
      if(s[u/S][0]==1)
	    check(u-S+1);
    }
  }

  else{
    if(u/S==0){
      if(s[1][u]==1)
	    check(S+u);
      if(s[1][u-1]==1)
	    check(S+u-1);
      if(s[0][u-1]==1)
	    check(u-1);
      if(s[S-1][u]==1)
	    check(S*(S-1)+u);
      if(s[S-1][u+1]==1)
	    check(S*(S-1)+u+1);
      if(s[0][u+1]==1)
	    check(u+1);
    }
    else if(u/S==S-1){
      if(s[0][u-S*(S-1)]==1)
	    check(u-S*(S-1));
      if(s[0][u-S*(S-1)-1]==1)
	    check(u-S*(S-1)-1);
      if(s[S-1][u-S*(S-1)-1]==1)
	    check(u-1);
      if(s[S-2][u-S*(S-1)]==1)
	    check(u-S);
      if(s[S-2][u-S*(S-1)+1]==1)
	    check(u-S+1);
      if(s[S-1][u-S*(S-1)+1]==1)
	    check(u+1);
    }
    else{
      if(s[u/S+1][u%S]==1)
	    check(u+S);
      if(s[u/S+1][u%S-1]==1)
	    check(u+S-1);
      if(s[u/S][u%S-1]==1)
	    check(u-1);
      if(s[u/S-1][u%S]==1)
	    check(u-S);
      if(s[u/S-1][u%S+1]==1)
	    check(u-S+1);
      if(s[u/S][u%S+1]==1)
	    check(u+1);
    }
  }

  return np-inp;

}

int getsize(int s[][S], int k){
  int nip=getatom(s,k);
  if(nip!=0){
    size[n_i]+=nip;
    int niip=np;
    for(int i=0;i<nip;i++) 
      getsize(s,ip[trial][niip-nip+i]);
  }
  else
    return 0;
}

double dis(int u, int v){

  int a=u%S;
  int b=u/S;
  int c=v%S;
  int d=v/S;
  int x=a-c;
  int y=b-d;
  double sr=0.;

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

  return bl*sqrt(sr);
}

int isboundary(int s[][S],int u){

  int bn=0;
  int i=u%S;
  int j=u/S;

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

  if(bn==6) 
    return 0;
  else
    return 1;

}
