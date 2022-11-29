## ----eval=FALSE---------------------------------------------------------------
#  List tree_create(NumericMatrix x, IntegerVector y, int q=NA_LOGICAL,
#                   int s_min=1, bool random_split=true){
#    # the split will be ended if the sample of the child is less than s_min
#    #  x a n*p matrix, containing samples
#  
#    int nvar=x.ncol();
#    if(q==NA_LOGICAL) q=floor(sqrt(nvar));
#  
#    int nsplit=y.length();
#    IntegerMatrix where(2,nsplit);  //recording where leaf the i-th sample lie in
#  
#    for(int i=0;i<nsplit;i++){
#      where(0,i)=i;
#      where(1,i)=0;
#    }
#  
#    # initialing
#    IntegerVector node={0};
#    IntegerVector parent={0};
#    IntegerVector left_child={-1};
#    IntegerVector right_child={-1};
#    IntegerVector dim={-1};
#    IntegerVector status={NODE_TOSPLIT};
#    NumericVector split = {0};
#    IntegerVector num_1 = {sum(y)};
#    IntegerVector num_0 = {nsplit-num_1[0]};
#  
#    int n1s = sum(y);
#    int n0s = y.length()-n1s;
#    int label0;
#    if(n0s==n1s){
#      label0 = RAND_PREDICT;
#    }else{
#      label0 = (n1s>n0s);
#    }
#    IntegerVector label = {label0};
#  
#    List node_info = List::create(Named("node")=node,Named("parent")=parent,
#                      Named("left_child")=left_child,Named("right_child")=right_child,
#                      Named("split_dim")=dim,Named("split_pos")=split,
#                      Named("node_status")=status,Named("num_0")=num_0,
#                      Named("num_1")=num_1,Named("label")=label);
#  
#    IntegerVector dimtemp = seq(0,nvar-1);
#    IntegerVector dim_for_split=sample(dimtemp, q, false);
#  
#    int count=1;
#    while(count){
#  
#      List tree_now = cnode_split(x,y,node_info,where,dim_for_split,q,s_min,random_split);
#      node_info = tree_now["node_info"];
#      IntegerMatrix C = tree_now["where"];
#      where = clone(C);
#      count = tree_now["count"];
#      node=node_info["node"];
#  
#    }
#  
#    IntegerVector dim_now = node_info["split_dim"];
#    dim_now = dim_now+1;
#    node_info["split_dim"] = dim_now;
#  
#    return(node_info);
#  }

## ----eval=F-------------------------------------------------------------------
#  int tree_predict(List node_info, NumericVector x_test){
#    IntegerVector node_status = node_info["isleaf"];
#    IntegerVector dim_ori = node_info["split_dim"];
#    NumericVector split= node_info["split_pos"];
#    IntegerVector left_child = node_info["left_child"];
#    IntegerVector right_child = node_info["right_child"];
#    IntegerVector label=node_info["label"];
#    int node_now = 0;
#    IntegerVector dim = clone(dim_ori)-1;
#    while(node_status[node_now]!=NODE_TERMINAL){
#      if(x_test[dim[node_now]]<=split[node_now]){
#        node_now = left_child[node_now];
#      }else{
#        node_now = right_child[node_now];
#      }
#    }
#  
#    int ans;
#    if(label[node_now]==RAND_PREDICT){
#      ans=unif_rand();
#    }else{
#      ans=label[node_now];
#    }
#  
#    return (ans);
#  }
#  

## ----eval=TRUE----------------------------------------------------------------
data(iris)
iris_data <- as.matrix(iris[,-5])
iris_label <- as.numeric(iris[,5])
iris_label[iris_label!=1]<-0 # transformed to binary classification
n <- length(iris_label)
iris_ind <- 1:n

# separate training and testing set
set.seed(22003)
tr_ind <- sample(iris_ind,100) # training set
te_ind <- iris_ind[-tr_ind]  # testing set
p <- ncol(iris_data)

# build the tree
t1 <- tree_create(x=iris_data[tr_ind,],y=iris_label[tr_ind],q=p,random_split=F)

# show the tree
tm <- matrix(unlist(t1),ncol=10)
colnames(tm) <- names(t1)
knitr::kable(tm)

## -----------------------------------------------------------------------------
test_pre <- sapply(1:50,function(i) tree_predict(t1,iris_data[te_ind[i],]))
cat("The test accuracy is\n")
mean(test_pre==iris_label[te_ind])

## -----------------------------------------------------------------------------
# Generate simulation data
n <- 50
p <- 400
set.seed(155)
x <- matrix(rnorm(n*p),nrow=n)
y <- as.numeric(x[,1]+x[,10]>0)


# create a random-split tree
set.seed(18)
t2 <- tree_create(x=x,y=y)

# show the tree
tm <- matrix(unlist(t2),ncol=10)
colnames(tm) <- names(t2)
knitr::kable(tm)

