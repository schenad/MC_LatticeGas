#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#define S (200)
#define trials (10)

using namespace std;

int main(){

  int samples[trials][S][S];
  char imanum[3];
  char filename[30];
  char firstname[20];
  int count=0;
  double coverage;
  cout<<"The prefix of the files to be combined is:"<<endl;

  fstream check;
  while(true){
  cin>>firstname;
  strcpy(filename,firstname);
  strcat(strcat(strcat(filename,"-"),"1"),".txt");
  check.open(filename,ios::in);
  if(!check)
    cout<<"The file is not existed. Please try again. "<<endl<<endl<<"The prefix of the files to be combined is:"<<endl;
  else
    break;
  }

  cout<<"Processing..............."<<endl<<endl<<"File:"<<endl<<endl;

  for(int k=0;k<trials;k++){
    strcpy(filename,firstname);
    sprintf(imanum,"%d",k+1);
    strcat(strcat(strcat(filename,"-"),imanum),".txt");
    cout<<filename<<endl;
    ifstream fin(filename);
    for(int i=0;i<S;i++)
       for(int j=0;j<S;j++)
         fin>>samples[k][i][j];
    fin.close();
  }

  for(int j=0;j<S;j++)            //Check coverage
    for(int i=0;i<S;i++)
      if(samples[0][j][i]==1)
    	  count++;
  coverage=count*100./(S*S);
  cout<<endl<<"(Coverage="<<coverage<<'%'<<')'<<endl;

  if(coverage>50.){
	  for(int trial=0;trial<trials;trial++)
	   for(int j=0;j<S;j++)
	     for(int i=0;i<S;i++){
	       if(samples[trial][j][i]==0)
		 samples[trial][j][i]=1;
	       else
		 samples[trial][j][i]=0;
	     }
   }

  ofstream mcfile;
  mcfile.open("Br.txt");
  for(int trial=0;trial<trials;trial++){
   for(int j=0;j<S;j++){            //Output results
     for(int i=0;i<S;i++)
       mcfile<<samples[trial][j][i]<<" ";
     mcfile<<endl;
     }
    mcfile<<endl;
  }
  mcfile.close();

  cout<<endl<<"have been combined to Br.txt."<<endl;
  return 0;
}
