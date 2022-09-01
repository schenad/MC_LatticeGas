#include <iostream>
#include<fstream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <ctime>
#define steps (1e6)
#define a (1.)     //Size of 1 unit of plane(Bond length)
#define N (100)  //Size of the plane
#define n (100)   //Number of atoms

using namespace std;

int ran(int x){
  return (rand()%x);
}

int main(){

  int sample[N][N];
  int ran(int);
  int x,y,xp,yp;
  void print(int s[][N]);
  srand((int)time(NULL));

  for(int i=0;i<N;i++)
    for(int j=0;j<N;j++)
      sample[j][i]=0;           //Initialize the plane.

    ifstream fin("0.txt");
    for(int i=0;i<N;i++)
       for(int j=0;j<N;j++) 
         fin>>sample[i][j];
    fin.close();

     for(int j=0;j<N;j++){
      for(int i=0;i<N;i++)
         cout<<sample[j][i]<<"   ";
      cout<<endl<<endl<<setw(2+j*2)<<' ';
     }
     cout<<endl;

  print(sample);

  return 0;

 }

void print(int s[][N]){

  int i,k;

  for(int j=0;j<N;j++){
    for(i=0-j/2,k=0;k<N;k++,i++){
      if(i<0)
        i+=N;
      if(i>=N)
	i-=N;     
    cout<<s[j][i]<<" ";
    }
    cout<<endl;
  }
  cout<<endl;
}
