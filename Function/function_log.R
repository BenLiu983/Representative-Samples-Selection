# Function 1: create a function to select samples from logistic regression model
am_log <- function(DATA,lambda,p1_low,select_row) 
{
  # DATA: the original dataset
  # lambda: a threshold to judge whether a sample is representative 
  # p1_low: a percentage, if the percentage of final representative samples < p1_low,
        # exclude these final samples from the original dataset, then restart the sample selection
  # select_row: the order number of rows selected from the original dataset "DATA"
  
  N=nrow(DATA)
  p= ncol(DATA)
  p1_hat =0
  X=DATA[,-p]
  Y=DATA[,p]
  
  # if no given sample, choose half of dataset randomly
  if (sum(select_row)==0)
  {
    select_row = sample(N, floor(N/2))
    x_step = DATA[select_row,] 
  }
  
  else
  {
    x_step = DATA[select_row,]
  }
  
  # if the percentage of the selected representative samples > p1_low, stop1 = 1, end the outer loop
  stop1=0
  
  # if the updated data remain the same, stop2 = 1, end the inner loop 
  stop2=0
  
  # max times of iteration 
  max_iter=100
  i=0
  
  while (stop1==0) 
  {
    while (stop2 == 0 && i <= max_iter)
    {
      # compute beta from the selected samples
      #df = data.frame(y=x_step[,p],x=x_step[,1:(p-1)])
      #beta=glm( x_step[,p]~x_step[,1:(p-1)],data=df,family="binomial")
      
      # ridge 
      x=x_step[,-p]
      y=x_step[,p]
      
      #na_index <- is.na(Y)   # in case Y is na
      cv=cv.glmnet(x, y, alpha=0,nfolds=3, family="binomial")
      betacv = coef(cv, s = "lambda.min")        # or "lambda.1se"
      
      # compute probability of each row of the original dataset
      x_all_plus1 = as.matrix(cbind(rep(1,nrow(DATA)),X))
      
      z_all =  x_all_plus1 %*% betacv
      pr_all = 1/(1+exp(-z_all))  
      
      devia = as.vector(2*log((Y/pr_all)^Y) + 2*log(((1-Y)/(1-pr_all))^(1-Y)))
      
      #compute the density of each row of the original dataset
      #den =as.vector(pi_star * (pr_all^(DATA[,p])) * (rep(1,N)-pr_all)^(rep(1,N)-DATA[,p]) )
      
      # select the rows, whose deviance is smaller than threshold
      num =which( devia<lambda, arr.ind = T)
      
      # update dataset
      x_stepnew=DATA[num,]
      
      # if the every row num of new data = row num of data before updating, stop the inner while
      if ((length(select_row)==length(num)) && all(select_row == num))  
        
      {  
        p1_hat = nrow(x_stepnew)/N
        stop2=1
      }
      else
      {
        # show the dimension of each step and upate x_step
        #print(dim(x_step))
        select_row=num
        x_step=DATA[select_row,]
      }
      
      i=i+1
    }
    
    # if selected representative samples's percentage > percentage set up, stop the whole outer while loop
    if ((p1_hat>p1_low) || (length(num)==0))
    {stop1=1}
    
    # else, exclude the selected samples, use the rest dataset to start over
    else
    {
      x_step= DATA[-num,]
      dataframe=as.data.frame(x_step)
      N1=nrow(x_step)
      x_step= x_step[sample(nrow(x_step), N1/2), ]
      stop2=0
      i=0
    }
  }
  
  # save some indexes
  assign("final_row",select_row,envir = globalenv())
  assign("final_beta",
         betacv,envir = globalenv())
  
  return (list(select_row, betacv))
}


# Function 2: create a function to select samples from logistic regression model with backtrack
fam_log <- function(DATA,lambda,p1_low,select_row) 
{
  # DATA: the original dataset
  # lambda: a threshold to judge whether a sample is representative 
  # p1_low: a percentage, if the percentage of final representative samples < p1_low,
        # exclude these final samples from the original dataset, then restart the sample selection
  # select_row: the order number of rows selected from the original dataset "DATA"
  
  N=nrow(DATA)
  p= ncol(DATA)
  DATA=as.matrix(DATA)
  X = DATA[,1:(p-1)]
  Y = DATA[,p]
  X_1 = cbind(rep(1,N),X)
  beta = matrix(rep(0,p),p,1)  
  p1_hat =0
  
  # if no given sample, choose half of dataset randomly
  if (sum(select_row)==0)
  {
    select_row = sample(N, floor(N/2))
    x_step = DATA[select_row,] 
  }
  
  else
  {
    x_step = DATA[select_row,]
  }
  
  # if the percentage of selected representative samples > p1_low, stop1 = 1, end the outer loop
  stop1=0
  
  # if the updated dataset remain the same, stop2 = 1, end the inner loop 
  stop2=0
  
  # max times of iteration 
  max_iter=100
  i=0
  
  while (stop1==0) 
  {
    while (stop2 == 0 && i <= max_iter)
    {
      # important parameter, depending on the gradient
      step = 3*10^(-4)
      
      maxiter2 = 100
      a = 0.2
      b = 0.9
      
      for (iter in 1:maxiter2)
      {
        beta_copy = beta
        
        # minimize the negative maximum likelihood
        X1 = x_step[,1:(p-1)]
        Y1 = x_step[,p]
        X1_1 = cbind(rep(1,nrow(x_step)),X1)
        
        # calculate the objective function
        q = sum(na.omit(Y1*log(1 + exp(-(X1_1 %*% beta))) + (1-Y1)*log(1 + exp(X1_1 %*% beta))) )
        
        # calculate gradient and updata beta
        grad = t(X1_1) %*% ( (1+exp(-(X1_1 %*% beta)))^(-1)  - Y1  )
        
        beta_new= beta - step* grad
        
        # updated objective
        q_new = sum(na.omit(Y1*log(1 + exp(-(X1_1 %*% beta_new))) + (1-Y1)*log(1 + exp(X1_1 %*% beta_new)))) 
        
        q_object = q - a*step*t(grad) %*% grad
        
        if (q_new > q_object)
        {
          step = step*b
        }
        
        else
        {
          beta = beta_new
          break
        }
      }
      
      # compute probability of each row of the original dataset
      x_all_plus1 = as.matrix(cbind(rep(1,nrow(DATA)),X))
      
      z_all =  x_all_plus1 %*% beta
      pr_all = 1/(1+exp(-z_all))  
      

      #compute the deviance of each row
      devia = as.vector(2*log((Y/pr_all)^Y) + 2*log(((1-Y)/(1-pr_all))^(1-Y)))

      # select the rows, whose deviance is smaller than threshold
      num =which(devia<lambda, arr.ind = T)
      
      # update dataset
      x_stepnew=DATA[num,]
      
      # if the every row num of new dataset = row num of data before updating, stop the inner while
      if ((length(select_row)==length(num)) && all(select_row == num))  
        
      {  
        p1_hat = nrow(x_stepnew)/N
        stop2=1
      }
      
      else
      {
        # show the dimension of each step and upate x_step
        #print(dim(x_step))
        select_row=num
        x_step=DATA[select_row,]
      }
      
      i=i+1
    }
    
    # if selected representative samples' percentage > percentage set up, stop the whole outer while loop
    if ((p1_hat>p1_low) || (length(num)==0))
    {stop1=1}
    
    # else, exclude the selected samples, use the rest dataset to start over
    else
    {
      x_step= DATA[-num,]
      dataframe=as.data.frame(x_step)
      N1=nrow(x_step)
      x_step= x_step[sample(nrow(x_step), N1/2), ]
      stop2=0
      i=0
    }
  }
  
  # save some indexes
  assign("final_row_bt",select_row,envir = globalenv())
  assign("final_beta_bt", beta,envir = globalenv())
  
  return (list(select_row, beta))
}



# Function 3: create a function to create static plot for each step
auto_plot_log <- function(DATA,lambda,p1_low,select_row) 
{
  # DATA: the original dataset
  # lambda: a threshold to judge whether a sample is representative 
  # p1_low: a percentage, if the percentage of final representative samples < p1_low,
           # exclude these final samples from the original dataset, then restart the sample selection
  # select_row: the order number of rows selected from the original dataset "DATA"
  
  N= nrow(DATA)
  p= ncol(DATA)
  
  # if no given samples, choose half of dataset randomly
  if (sum(select_row)==0)
  {
    select_row = sample(N, floor(N/2))
    x_step = DATA[select_row,] 
    dataframe = as.data.frame(x_step)
    
    # prepare for the ggplot
    colnames(dataframe) <- c("x","y","type")
    dataframe$type=as.factor(dataframe$type)
  }
  
  else
  {
    x_step = DATA[select_row,] 
    dataframe = as.data.frame(x_step)
    colnames(dataframe) <- c("x","y","type")
    dataframe$type=as.factor(dataframe$type)
  }
  
  # if the percentage of selected representative samples > p1_low, stop1 = 1, end the outer loop
  stop1=0
  
  # if the updated dataset remain the same, stop2 = 1, end the inner loop 
  stop2=0
  
  # max times of iteration 
  max_iter=100
  i=0
  
  while (stop1==0) 
  {
    while (stop2 == 0 && i <= max_iter)
    {
      # prepare for the plot when i >=1
      dataframe=as.data.frame(x_step)
      colnames(dataframe) <- c("x","y","type")
      dataframe$type=as.factor(dataframe$type)
      
      # compute beta from the selected samples 
      df = data.frame(y=x_step[,2],x=x_step[,1])
      beta=glm(x_step[,2]~x_step[,1],data=df,family="binomial")
      
      # compute the probability for each row of selected samples
      inter=beta$coefficients[1]
      slope=beta$coefficients[2]
      z = inter + slope*x_step[,1] 
      pr1 = 1/(1+exp(-z)) 
      
      library(ggplot2)
      cols <- c("1" = "red", "0" = "blue")
      
      # diaplay plot
      print(ggplot(dataframe, 
                   aes(x=x,  
                       y=y,  
                       fill =type))  +   
              xlim(-15,15)+ylim(0,1)+
              geom_point() +
              geom_point(aes(x=x, y=pr1), colour="gold")+
              geom_point(alpha = 1,shape=21, size =3)+ 
              scale_fill_manual(values = cols)+
              labs(title= paste("step" ,as.character(i)),envir = globalenv()))
      
      # select the rows, whose deviance is smaller than the threshold
      x_all_plus1 = as.matrix(cbind(rep(1,nrow(DATA)),DATA[,1]))
      
      z_all =  x_all_plus1 %*% beta$coefficients
      pr_all = 1/(1+exp(-z_all))  
      devia = as.vector(2*log((DATA[,2]/pr_all)^DATA[,2]) + 2*log(((1-DATA[,2])/(1-pr_all))^(1-DATA[,2])))
      
      num = which(devia<lambda, arr.ind = T) 
      
      # update dataset
      x_stepnew=DATA[num,]
      
      # if every row num of new dataset = row num of dataset before updating, stop the inner while
      if ((length(select_row)==length(num)) && all(select_row == num))  
      {  
        p1_hat = nrow(x_stepnew)/N
        stop2=1
      }
      else
      {
        #print(dim(x_step))
        select_row=num
        x_step=DATA[select_row,]
      }
      
      # save dataframe and beta of each step for gif
      assign(paste("df", i+1, sep = ""),
             dataframe,envir = globalenv() )
      assign(paste("beta", i+1, sep = ""),
             beta,envir = globalenv() )
      i=i+1
    }
    
    # if selected representative samples' percentage > percentage set up, stop the whole outer while
    if (p1_hat>p1_low)
    {stop1=1}
    
    # else, exclude the selected samples, use the rest dataset to start over
    else
    {
      x_step= DATA[-num,]
      N1=nrow(x_step)
      x_step= x_step[sample(nrow(x_step), N1/2), ]
      dataframe=as.data.frame(x_step)
      colnames(dataframe) <- c("x","y","type")
      dataframe$type=as.factor(dataframe$type)
      stop2=0
      i=0
    }
  }
  # save index
  assign("final_row",select_row,envir = globalenv())
  assign("final_beta",beta$coefficients,envir = globalenv())
}


# Function 4 : gif 

dynamic_plot_log <- function(dataframe_list,beta_list) 
{
  df_list1=dataframe_list
  b_list1=beta_list
  cols <- c("1" = "red", "0" = "blue")
  library(animation)
  library(ggplot2)
  i=0
  
  saveGIF({
    for(i in seq_along(df_list1)){
      # compute the probability of each row of selected samples
      z = b_list1[[i]]$coefficients[1] + b_list1[[i]]$coefficients[2] * df_list1[[i]][,1] 
      pr1 = 1/(1+exp(-z)) 
      
      a <- ggplot(df_list1[[i]], 
                  aes(x = x,  y = y, fill = type))  +   
        xlim(-15,15)+ylim(0,1)+
        geom_point() +
        geom_point(aes(x=x, y=pr1), colour="gold")+
        geom_point(alpha = 1,shape=21, size =3)+ 
        scale_fill_manual(values = cols) +
        labs(title= paste("step" ,as.character(i)),envir = globalenv())
      print(a)
    }
  }, interval = 1, movie.name="test.gif")
}


# Function 5:  Create a Compound Symmetry sigma
sigma_compound <- function(p,rou) 
{
  sig1= diag(rou,p,p)
  sig2= matrix(rou,p,p)
  sig=sig1+sig2
  return(sig)
}


# Function 6: Newton_log
newton_log <- function(X, Y, beta)
{
  n=nrow(X)

  for (i in 1:100)
  {
    beta0=beta
    prob = exp(X %*% beta)/(1+exp(X %*% beta))
    
    # gradient
    grad = t(X) %*% (prob-Y)
    
    # Hessian Matrix
    M = exp(X %*% beta) / (1+(exp(X %*% beta))^2)
    H = t(X) %*% diag(M[1],nrow(X),nrow(X)) %*% X
    #update beta
    beta = beta -solve(H) %*% grad
    # print(beta)
    
      if (norm(beta-beta0,type="2")<=0.01)
        break
    
  }
  beta_final=beta
  
  return (beta_final) 

}


# Function 7: Newton_log_ridge
newton_log_ridge <- function(X, Y, beta,lambda)
{
  n=nrow(X)
  
  for (i in 1:100)
  {
    beta0=beta
    prob = exp(X %*% beta)/(1+exp(X %*% beta))
    
    # gradient
    grad = t(X) %*% (prob-Y) + 2*lambda*beta
    
    # Hessian Matrix
    M = exp(X %*% beta) / (1+(exp(X %*% beta))^2)
    H = t(X) %*% diag(M[1],nrow(X),nrow(X)) %*% X+ 2*lambda*diag(1,ncol(X),ncol(X))
    #update beta
    beta = beta -solve(H) %*% grad
    # print(beta)
    
    if (norm(beta-beta0,type="2")<=0.01)
      break
    
  }
  beta_final=beta
  
  return (beta_final) 
  
}

# Function 8: backtrack 
backtrack_log<-function(data)
{
  data=as.matrix(data)
  n = nrow(data)
  p = ncol(data)
  X = data[,1:(p-1)]
  Y = data[,p]
  X_1 = cbind(rep(1,n),X)
  beta = matrix(rep(0,p),p,1)  
  
  # important parameter, depending on the gradient
  step = 3*10^(-2)
  maxiter = 2000
  a = 0.2
  b = 0.9
  e= 10^(-5)
  
  for (iter in 1:maxiter)
  {
    beta_copy = beta
    
    # minimize the negative maximum likelihood
    q = sum(na.omit(Y*log(1 + exp(-(X_1 %*% beta))) + (1-Y)*log(1 + exp(X_1 %*% beta))) )
    
    grad = t(X_1) %*% ( (1+exp(-(X_1 %*% beta)))^(-1) - Y)
    
    beta_new= beta - step* grad
    
    # updated objective
    q_new = sum(na.omit(Y*log(1 + exp(-(X_1 %*% beta_new))) + (1-Y)*log(1 + exp(X_1 %*% beta_new)))) 
    q_object = q - a*step*t(grad) %*% grad
    
    if(q_new > q_object)
    {
      step = step*b
    }
    else
    {
      beta = beta_new
      if (norm(beta_new-beta_copy, type="2")<e)
      {break}
    }
  }
  
  return(list(beta,iter))
}


