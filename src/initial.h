#include <Rcpp.h>
#include <R.h>
using namespace Rcpp;
using namespace std;
#define NODE_TERMINAL 1
#define NODE_TOSPLIT  2
#define NODE_INTERIOR 0
#define RAND_PREDICT 99



// Here are some basic functions used in the procedure of constructing a KS-test-tree.


//' @title A order-located function
//' @description A function returns the index (starting with 0) of the i-th smallest element in the vector in its i-th element
//' @param vec a numeric vector
//' @return an integer vector with i-th element equaling the index (starting with 0) of the i-the smallest element in the input vector
//' @examples
//' \dontrun{
//' x<-c(1,3,2,6)
//' ordered(x)}
//' @export
//[[Rcpp::export]]
IntegerVector ordered(NumericVector vec){
  int n = vec.size();
  IntegerVector vec_order(n);
  std::iota(vec_order.begin(),vec_order.end(),0); //Initializing
  stable_sort(vec_order.begin(),vec_order.end(), [&vec](int i,int j){return vec[i]<vec[j];}); //该值是第几小的
  return vec_order;
}



// calculate the criteria for split
List ccal_r(NumericMatrix x, IntegerVector y,int s_min){
  
  int n = y.length();
  int n_1 = sum(y);
  int p=x.ncol();
  int dim=-1;
  int pos=-1;
  int ntie=0;
  double impurity=DBL_MAX;
  
  NumericMatrix results(p, n-1);
  
  
  double x_now;
  
  
  
  for(int i=0;i<p;++i){
    
    IntegerVector ind = ordered(x(_,i)); 
    
    double sum_y = 0;
    
    
    for(int t=0;t<n-1;++t){
      
      x_now=x(ind[t],i);
      int t0=t;
      
      //find the split point
      for(int k=t+1;k<n;++k){
        
        if(x(ind[k],i)==x_now){
          t=k;
        }else{
          break;
        }
      }
      if(t==n-1){ 
        
        continue;}
      
      for(int k=t0;k<t+1;++k){
        sum_y +=y[ind[t]];
      }
      
      double aa=sum_y/(t+1); 
      double bb=(n_1-sum_y)/(n-t-1); 
      
      for(int j=0;j<t+1;++j){
        results(i,t)+=(y[ind[j]]-aa)*(y[ind[j]]-aa);
      }
      for(int j=t+1;j<n;++j){
        results(i,t)+=(y[ind[j]]-bb)*(y[ind[j]]-bb);
      }
      
      if (impurity > results(i,t)){
        
        impurity =results(i,t);
        dim=i;
        pos=t;
      }
      if(impurity==results(i,t)){
        ntie++;
        
        if(unif_rand()<1.0/ntie){
          dim = i;
          pos = t;
        }
      }
    }
  }
  
  
  if(dim==-1){
    //no split
    return List::create(Named("dim")=dim);
  }
  
  IntegerVector ind = ordered(x(_,dim));
  IntegerVector child(n);
  IntegerVector status(3);
  
  double best_split=(x(ind[pos],dim)+x(ind[pos+1],dim))/2.0;
  int n1l=0;
  int n1r=0;
  int nr=0;
  for(int i=0;i<n;++i){
    if(x(i,dim)>best_split){
      child[i]=1;//right child
      n1r+=y[i];
      nr++;
    }else{
      n1l+=y[i];
      child[i]=0;//left child
    }
  }
  
  int n0r=nr-n1r;
  int n0l=n-nr-n1l;
  
  
  // determine the status of children
  status[0] = NODE_INTERIOR;
  if(n0l*n1l>0){
    if(n0l+n1l>s_min){
      status[1]=NODE_TOSPLIT;
    }else{
      status[1]=NODE_TERMINAL;
    }
  }else{
    status[1]=NODE_TERMINAL;
  }
  
  if(n0r*n1r>0){
    if(n0r+n1r>s_min){
      status[2]=NODE_TOSPLIT;
    }else{
      status[2]=NODE_TERMINAL;
    }
  }else{
    status[2]=NODE_TERMINAL;
  }
  
  
  return List::create(Named("dim")=dim, Named("nl1")=n1l, Named("nl0")=n0l,
                            Named("nr1")=n1r, Named("nr0")=n0r, Named("split")=best_split,
                                  Named("status")=status,Named("child")=child);
  
}