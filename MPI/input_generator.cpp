#include <bits/stdc++.h>

using namespace std;

int main (){
  
 	

  for(int t=0;t<128;t++){
   string s = "MPIinputs/input" + to_string(t) + ".txt"; 
   const char * fname  = s.c_str();


   freopen (fname,"w",stdout);

	  for(int i=0;i<128;i++){
	    int x = rand() % 128;
	    cout << x << " ";
	  }
  }
  
  return 0;
}
