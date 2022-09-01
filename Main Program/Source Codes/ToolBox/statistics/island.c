#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#define S (100)  //Size of the plane

using namespace std;

const int n=S*S;    //Number of atoms approximately(can be more than the exact number, but cannot be less)
int size[n];
int n_i=0;
int ip[n];
int np=0;

int main(){

  int sample[S][S];
  int getsize(int s[][S],int);
  int check(int);

  ifstream fin("0.txt");
    for(int i=0;i<S;i++)
       for(int j=0;j<S;j++) 
         fin>>sample[i][j];
  fin.close();

  for(int pt=0;pt<S*S;pt++){
    if(*(&sample[0][0]+pt)==1)
      if(check(pt)==1){
	size[n_i]++;
        getsize(sample,pt);
        n_i++;
      }
  }

  ofstream mcfile;
  mcfile.open("size.dat");
 
  cout<<"No. of islands: "<<n_i<<endl;
  for(int i=0;i<n_i;i++)
    mcfile<<size[i]<<endl;
  cout<<"No. of points: "<<np<<endl;
  //  for(int j=0;j<np;j++)
  //    cout<<"Point "<<ip[j]<<endl;

  mcfile.close();

  return 0;
}



int check(int w){
  for(int i=0;i<np;i++){
    if(w==ip[i])
	return 0;
  }
  np++;
  ip[np-1]=w;
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
      getsize(s,ip[niip-nip+i]);
  }
  else
    return 0;
}
