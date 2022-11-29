#include <Rcpp.h>
#include <R.h>
#include "initial.h"
using namespace Rcpp;
using namespace std;
#define NODE_TERMINAL 1
#define NODE_TOSPLIT  2
#define NODE_INTERIOR 0
#define RAND_PREDICT 99



List cnode_split(NumericMatrix x_split, IntegerVector y_split, List node_info, IntegerMatrix where,
                 IntegerVector dim_for_split, int q, int s_min, bool random_split=true){
 
  
  IntegerVector node=node_info["node"];
  IntegerVector parent=node_info["parent"];
  IntegerVector left_child=node_info["left_child"];
  IntegerVector right_child=node_info["right_child"];
  IntegerVector dim=node_info["split_dim"];
  NumericVector split=node_info["split_pos"];
  IntegerVector num_0=node_info["num_0"];
  IntegerVector num_1=node_info["num_1"];
  IntegerVector node_status=node_info["node_status"];
  IntegerVector label=node_info["label"];
  
  IntegerVector node_split;//nodes to split
  int count=0;
  
  IntegerVector dimtemp=seq(0,x_split.ncol()-1); 
  int ncur=node.length();
  
  for(int i=0;i<ncur;++i){
    if(node_status[i]==NODE_TOSPLIT){
      if(random_split) dim_for_split=sample(dimtemp,q,false);
      
      NumericMatrix x_in(sum(where(1,_)==i),q);
      IntegerVector y_in;
      int l=0;
      
      IntegerVector sam_now;
      
      for(int k=0;k<where.ncol();++k){
       
        if(where(1,k)==i){
          for(int m=0;m<q;++m){
            x_in(l,m)=x_split(k,dim_for_split[m]);
          }
          
          sam_now.push_back(k);
          y_in.push_back(y_split[k]);
          l++;
        }
      }
      
      
      List rs = ccal_r(x_in,y_in,s_min);
      
      int dim_now = rs["dim"];
      if(dim_now==-1){
        
        node_status[i]=NODE_TERMINAL;
      }else{
        IntegerVector statuss=rs["status"];
        
        
        IntegerVector child= rs["child"];
        node_status[i]=statuss[0];
        dim[i] = dim_for_split[dim_now];
        split[i]=rs["split"];
        left_child[i] = node.length();
        right_child[i] = node.length()+1;
        node.push_back(node.length());
        node.push_back(node.length());
        parent.push_back(i);
        parent.push_back(i);
        left_child.push_back(-1);
        left_child.push_back(-1);
        right_child.push_back(-1);
        right_child.push_back(-1);
        dim.push_back(-1);
        dim.push_back(-1);
        split.push_back(0);
        split.push_back(0);
        node_status.push_back(statuss[1]);
        node_status.push_back(statuss[2]);
        num_0.push_back(rs["nl0"]);
        num_0.push_back(rs["nr0"]);
        num_1.push_back(rs["nl1"]);
        num_1.push_back(rs["nr1"]);

        
        if(num_0[node.length()-2]==num_1[node.length()-2]){
          label.push_back(RAND_PREDICT);
        }else{
          label.push_back(int(num_0[node.length()-2]<num_1[node.length()-2]));
        }
        
        
        if(num_0[node.length()-1]==num_1[node.length()-1]){
          label.push_back(RAND_PREDICT);
        }else{
          label.push_back(int(num_0[node.length()-1]<num_1[node.length()-1]));
        }
        
        //move data
        
        
        for(int s=0;s<sam_now.length();s++){
          if(child[s]) where(1,sam_now[s])=right_child[i];
          else where(1,sam_now[s])=left_child[i];
        }
        
        count+=1;
      }
    }
    
  }

  node_info["node"]=node;
  node_info["parent"]=parent;
  node_info["left_child"]=left_child;
  node_info["right_child"]=right_child;
  node_info["split_dim"]=dim;
  node_info["node_status"]=node_status;
  node_info["split_pos"]=split;
  node_info["num_0"]=num_0;
  node_info["num_1"]=num_1;
  node_info["label"]=label;
  return List::create(Named("node_info")=node_info,Named("where")=where,Named("count")=count);
  
}



//' @title Predict the label of a sample using our tree
//' @description Predict the label of a sample using our tree
//' @param x_test A vector, the sample 
//' @param node_info A tree producted by tree_create
//' @return The predicted label of x_test
//' @examples
//' \dontrun{
//' data(iris)
//' iris_data <- as.matrix(iris[1:100,-5])
//' iris_label <- as.numeric(iris[1:100,5])-1
//' set.seed(15)
//' ind <- 1:100
//' ind_tr <- sample(ind,60)
//' ind_test <- ind[-ind_tr]
//' t <- tree_create(x=iris_data,y=iris_label,q=4,random_split=F)
//' pre <- sapply(1:40,tree_predict(t,iris_data[ind_test[i],]))
//' }
//' @export
//[[Rcpp::export]] 
int tree_predict(List node_info, NumericVector x_test){
  IntegerVector node_status = node_info["node_status"];
  IntegerVector dim_ori = node_info["split_dim"];
  NumericVector split= node_info["split_pos"];
  IntegerVector left_child = node_info["left_child"];
  IntegerVector right_child = node_info["right_child"];
  IntegerVector label=node_info["label"];
  int node_now = 0;
  IntegerVector dim = clone(dim_ori)-1;
  while(node_status[node_now]!=NODE_TERMINAL){
    if(x_test[dim[node_now]]<=split[node_now]){
      node_now = left_child[node_now];
    }else{
      node_now = right_child[node_now];
    }
  }
  
  int ans;
  if(label[node_now]==RAND_PREDICT){
    ans=unif_rand();
  }else{
    ans=label[node_now];
  }
  
  return ans;
}


//' @title A random-split method to build a tree for binary classification
//' @description We use recursive partition to build a tree for binary classification problem. We take the mean-square error instead of gini as the split criteria. We can choose to build a random-split tree or a deterministic tree by setting the parameter \code{q} and \code{random_split}. This function is the base for a random forest.
//' @param x Numeric matrix, with \code{n} rows and \code{p} columns.
//' @param y The binary label, either 0 or 1.
//' @param q The number of dimensions of \code{x} to be considered for split at each node. (Default q=\code{sqrt(p)})
//' @param s_min If the number of sample on the node is less than s_min (default 1), latter split will be ended.
//' @param random_split \code{T} or \code{F}. \code{T} means that we would like to randomly choose \code{q} dimensions of \code{x} to split at each node. \code{F} means we just split in q dimensions randomly chosen at the beginning.
//' @return A tree saved as a list contains information
//' @examples
//' \dontrun{
//' data(iris)
//' iris_data <- as.matrix(iris[1:100,-5])
//' iris_label <- as.numeric(iris[1:100,5])-1
//' set.seed(15)
//' ind <- 1:100
//' ind_tr <- sample(ind,60)
//' ind_test <- ind[-ind_tr]
//' t <- tree_create(x=iris_data,y=iris_label,q=4,random_split=F)
//' pre <- sapply(1:40,tree_predict(t,iris_data[ind_test[i],]))
//' }
//' @export
//[[Rcpp::export]]
List tree_create(NumericMatrix x, IntegerVector y, int q=NA_LOGICAL, 
                 int s_min=1, bool random_split=true){
  //the split will be ended if the sample of the child is less than s_min
  // x a n*p matrix, containing samples

  int nvar=x.ncol();
  if(q==NA_LOGICAL) q=floor(sqrt(nvar)); 
  
  int nsplit=y.length();
  IntegerMatrix where(2,nsplit);  //recording where leaf the i-th sample lie in
 
  for(int i=0;i<nsplit;i++){
    where(0,i)=i;
    where(1,i)=0;
  }
  
  
  //initialing
  IntegerVector node={0};
  IntegerVector parent={0};
  IntegerVector left_child={-1};
  IntegerVector right_child={-1};
  IntegerVector dim={-1};
  IntegerVector status={NODE_TOSPLIT};
  NumericVector split = {0};
  IntegerVector num_1 = {sum(y)};
  IntegerVector num_0 = {nsplit-num_1[0]};
  
  int n1s = sum(y);
  int n0s = y.length()-n1s;
  int label0;
  if(n0s==n1s){
    label0 = RAND_PREDICT;
  }else{
    label0 = (n1s>n0s);
  }
  IntegerVector label = {label0};
  
  List node_info = List::create(Named("node")=node,Named("parent")=parent,Named("left_child")=left_child,
                                      Named("right_child")=right_child,Named("split_dim")=dim,
                                      Named("split_pos")=split,Named("node_status")=status,Named("num_0")=num_0,
                                      Named("num_1")=num_1,Named("label")=label);
 
  IntegerVector dimtemp = seq(0,nvar-1);
  IntegerVector dim_for_split=sample(dimtemp, q, false);
  
  int count=1;
  while(count){
    
    List tree_now = cnode_split(x,y,node_info,where,dim_for_split,q,s_min,random_split);
    node_info = tree_now["node_info"];
    IntegerMatrix C = tree_now["where"];
    where = clone(C);
    count = tree_now["count"];
    node=node_info["node"];
    
  }
 
  IntegerVector dim_now = node_info["split_dim"];
  dim_now = dim_now+1;
  node_info["split_dim"] = dim_now;
  
  return node_info;
}