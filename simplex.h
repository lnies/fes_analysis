/*
 Copyright (C) 2010 Botao Jia


 This file is an implementation of the downhill simplex optimization algorithm using C++.
 To use BT::Simplex correctly, the followings are needed, inclusively.

 1. f: a function object or a function which takes a vector<class Type> and returns a Type, inclusively.
       Signature, e.g. for double Type, 
            double f(vector<double> x);

 2. init: an inital guess of the fitted parameter values which minmizes the value of f. 
          init must be a vector<D>, where D must be the exactly same type as the vector taken by f.
          init must have the exactly same dimension as the vector taken by f.
	  init must order the parameters, such that the order follows the vector taken by f.
          e.g. f takes vector<double> x, where x[0] represents parameter1; x[1] represents parameter2, etc.
          And init must follow this order exactly, init[0] is the initial guess for parameter1,
          init[1] is the initial guess for parameter2, etc.

3 to 5 are all optional:

3. tol: the termination criteria. 
        It measures the difference of the simplex centroid*(N+1) of consecutive iterations.

4. x:  an initial simplex, which is calculated according to the initial trial parameter values.

5. iterations: maximum iterations.

The return value of BT::Simplex is a vector<Type>, 
which represents the optimized parameter values that minimize f.
The order of the parameter is the same as in the vector<Type> init.


 You can redistribute it and/or modify it at will.

 This program is distributed in the hope that it will be useful, but WITHOUT
 ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 FOR A PARTICULAR PURPOSE.  
*/

//simplex.h
#include <vector>
#include <limits>
#include <algorithm>
#include <functional>
#include <iostream>

using namespace std;

namespace BT{

  template<class D, class OP>
  vector<D> Simplex(OP f,                       // target function
       vector<D> &array,                         // Data array to fit
       vector<D> &proto,                         // Proto array to adjust
       vector<D> &weights,                       // Weight array weight smaples
       int start,                                    // Start sample for comparing vector
       int end,                                      // End sample for comparing vector
			 vector<D> &init,                          // initial guess of the parameters
			 D tol=1E8*numeric_limits<D>::epsilon(),  // termination criteria
			 vector<vector<D> > x =  vector<vector<D> >(), //x: The Simplex
			 int iterations=1E5,                           // Iteration step number
       int debug = 1                                 // Debug output
    ){    

    int N=init.size();                         //space dimension
    const double a=2.0, b=2.0, g=0.5, h=0.5;   //coefficients
                                               //a: reflection  -> xr  
                                               //b: expansion   -> xe 
                                               //g: contraction -> xc
                                               //h: full contraction to x1
    vector<D> xcentroid_old(N,0);   //simplex center * (N+1)
    vector<D> xcentroid_new(N,0);   //simplex center * (N+1)
    vector<D> vf(N+1,0);            //f evaluated at simplex vertexes       
    int x1=0, xn=0, xnp1=0;         //x1:   f(x1) = min { f(x1), f(x2)...f(x_{n+1} }
                                    //xnp1: f(xnp1) = max { f(x1), f(x2)...f(x_{n+1} }
                                    //xn:   f(xn)<f(xnp1) && f(xn)> all other f(x_i)
    int cnt=0; //iteration step number

    vector<D> subst; // substitute vector for the adjust function

    if(x.size()==0) //if no initial simplex is specified
    { //construct the trial simplex
      //based upon the initial guess parameters
      vector<D> del( init );
      transform(del.begin(), del.end(), del.begin(), 
	      bind2nd( divides<D>() , 20) );//'20' is picked 
                                                           //assuming initial trail close to true
    
      for(int i=0; i<N; ++i){
      	vector<D> tmp( init );
      	tmp[i] +=  del[i];
      	x.push_back( tmp );
      }
      x.push_back(init);//x.size()=N+1, x[i].size()=N
      
      //xcentriod
      transform(init.begin(), init.end(), 
		    xcentroid_old.begin(), bind2nd(multiplies<D>(), N+1) );
    }//constructing the simplex finished
    
    //optimization begins
    for(cnt=0; cnt<iterations; ++cnt){

      for(int i=0;i<N+1;++i){
        // subst = array_adjust(proto, x[i], debug);
	      vf[i]= f(array, proto, weights, x[i], start, end, debug); // x[i][1] is the shift parameter
        // vf[i]= f(x[i]);

      }
      
      x1=0; xn=0; xnp1=0;//find index of max, second max, min of vf.
      
      for(int i=0;i<(int)vf.size();++i){
      	if(vf[i]<vf[x1]){
      	  x1=i;
      	}
      	if(vf[i]>vf[xnp1]){
      	  xnp1=i;
      	}
      }
      
      xn=x1;
      
      for(int i=0; i<(int)vf.size();++i){ 
      	if( vf[i]<vf[xnp1] && vf[i]>vf[xn] ) xn=i;
      }
      //x1, xn, xnp1 are found

      vector<D> xg(N, 0);//xg: centroid of the N best vertexes
      for(int i=0; i<(int)x.size(); ++i){
      	if(i!=xnp1)
      	  transform(xg.begin(), xg.end(), x[i].begin(), xg.begin(), plus<D>() );
      }
      transform(xg.begin(), xg.end(), 
		    x[xnp1].begin(), xcentroid_new.begin(), plus<D>());
      transform(xg.begin(), xg.end(), xg.begin(), 
		    bind2nd(divides<D>(), N) );
      //xg found, xcentroid_new updated

      //termination condition
      D diff=0;          //calculate the difference of the simplex centers
                         //see if the difference is less than the termination criteria
      for(int i=0; i<N; ++i)     
	       diff += fabs(xcentroid_old[i]-xcentroid_new[i]);

      if (diff/N < tol) break;              //terminate the optimizer
      else xcentroid_old.swap(xcentroid_new); //update simplex center
      
      //reflection:
      vector<D> xr(N,0); 
      for( int i=0; i<N; ++i)
	       xr[i]=xg[i]+a*(xg[i]-x[xnp1][i]);
      //reflection, xr found
      
      // D fxr=f(xr);//record function at xr
      // subst = array_adjust(proto, xr, debug);
      D fxr=f(array, proto, weights, xr, start, end, debug);//record function at xr
      
      if(vf[x1]<=fxr && fxr<=vf[xn])
    	copy(xr.begin(), xr.end(), x[xnp1].begin() );
          
          //expansion:
          else if(fxr<vf[x1]){
    	vector<D> xe(N,0);
    	for( int i=0; i<N; ++i)
    	  xe[i]=xr[i]+b*(xr[i]-xg[i]);
    	// if( f(xe) < fxr )
      // subst = array_adjust(proto, xe, debug);
      if( f(array, proto, weights, xe, start, end, debug) < fxr )
    	  copy(xe.begin(), xe.end(), x[xnp1].begin() );
    	else
    	  copy(xr.begin(), xr.end(), x[xnp1].begin() );
          }//expansion finished,  xe is not used outside the scope
          
      //contraction:
          else if( fxr > vf[xn] ){
    	vector<D> xc(N,0);
    	for( int i=0; i<N; ++i)
    	  xc[i]=xg[i]+g*(x[xnp1][i]-xg[i]);
    	// if( f(xc) < vf[xnp1] )
      // subst = array_adjust(proto, xc, debug);
      if( f(array, proto, weights, xc, start, end, debug) < vf[xnp1] )
    	  copy(xc.begin(), xc.end(), x[xnp1].begin() );
    	else{
    	  
    	  for( int i=0; i<(int)x.size(); ++i ){
    	    if( i!=x1 ){ 
    	      for(int j=0; j<N; ++j) 
    		x[i][j] = x[x1][j] + h * ( x[i][j]-x[x1][j] );
    	    }
    	  }
    	  
    	}
      }//contraction finished, xc is not used outside the scope
      // subst = array_adjust(proto, x[x1], debug);
      if (debug == 2) printf("%3.5f %3.5f %3.5f\n", x[x1][0], x[x1][1], f(array, proto, weights, x[x1], start, end, debug)); 
    }//optimization is finished

    if(cnt==iterations){//max number of iteration achieves before tol is satisfied
      cout<<"Iteration limit achieves, result may not be optimal"<<endl;
    }
    if (debug == 2) printf("+++++ %3.3f %3.3f +++++\n", x[x1][0], x[x1][1]);
    return x[x1];
  }
    
}//BT::
