#include <vector>
#include <map>
#include <math.h>
using namespace std;





///////////////////////////////////////////////////////////////////////////
// Copyright (C) 2018  Bowen Liu                                         //
//                                                                       //
// This program is free software: you can redistribute it and/or modify  //
// it under the terms of the GNU General Public License as published by  //
// the Free Software Foundation, either version 3 of the License, or     //
// (at your option) any later version.                                   //
//                                                                       //
// This program is distributed in the hope that it will be useful,       //
// but WITHOUT ANY WARRANTY; without even the implied warranty of        //
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         //
// GNU General Public License for more details.                          //
//                                                                       //
// You should have received a copy of the GNU General Public License     //
// along with this program.  If not, see <http://www.gnu.org/licenses/>. //
///////////////////////////////////////////////////////////////////////////


/* 
AnDarl.c
$Revision: 1.1 $  $Date: 2014/06/09 07:53:21 $
Original C code by G. and J. Marsaglia
R interface by Adrian Baddeley
https://github.com/cran/goftest
*/
static double adinf(double z) { 
  if(z<2.) return (
      exp(-1.2337141/z)/sqrt(z)
  )*(
      2.00012+(.247105-	
      (.0649821-
      (.0347962-
      (.011672-.00168691*z)
         *z)*z)*z)*z);
  /* max |error| < .000002 for z<2, (p=.90816...) */
  return exp(
    -exp(1.0776-(2.30695-(.43424-(.082433-(.008056 -.0003146*z)
                                    *z)*z)*z)*z));
                                    /* max |error|<.0000008 for 4<z<infinity */
}



static double AD(int n,double z){
  double c,v,x;
  x=adinf(z);
  /* now x=adinf(z). Next, get v=errfix(n,x) and return x+v; */
  if(x>.8) {
    v=(-130.2137+(745.2337-(1705.091-(1950.646-(1116.360-255.7844*x)
                                        *x)*x)*x)*x)/n;
    return x+v;
  }
  c=.01265+.1757/n;
  if(x<c){ 
    v=x/c;
    v=sqrt(v)*(1.-v)*(49*v-102);
    return(x+v*(.0037/(n*n)+.00078/n+.00006)/n);
  }
  v=(x-c)/(.8-c);
  v=-.00022633+(6.54034-(14.6538-(14.458-(8.259-1.91864*v)*v)*v)*v)*v;
  return (x+v*(.04213+.01365/n)/n);
}



static double ADtest(int n, double *x)
{ int i;
  double t,z=0;
  for(i=0;i<n;i++)   {
    t=x[i]*(1.-x[n-1-i]);
    z=z-(i+i+1)*log(t);
  }
  return AD(n,-n+z/n);
}

