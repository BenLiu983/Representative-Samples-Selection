# Function 0: create a function to select samples from Poisson regression model
am_poi <- function(DATA,lambda,p1_low,select_row) 
{
  # DATA: the original dataset
  # lambda: a threshold to judge whether a sample is representative 
  # p1_low: a percentage, if the percentage of final representative samples < p1_low,
         # exclude these final samples from the original dataset, then restart the sample selection
  # select_row: the order number of rows selected from the original dataset "DATA"
  
  N=nrow(DATA)
  p= ncol(DATA)
  X = DATA[,-p]
  Y = DATA[,p]
  
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
      # compute beta from the selected represenatative samples
      data <- data.frame(y=x_step[,p], x=x_step[,-p])             
      beta <- glm(x_step[,p] ~ x_step[,-p], data=data, family=poisson(link="log"))
      
      # compute predicted y 
      DATA1= cbind(rep(1,N),DATA)
      y_pre = exp((DATA1[,-(p+1)])%*% beta$coefficients)
      
      #deviance
      devia = as.vector(2*(y_pre-Y)+2*Y*log(Y/y_pre))
      
      # select the rows, whose deviance is smaller than the threshold
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
  assign("final_row",select_row,envir = globalenv())
  assign("final_beta",beta,envir = globalenv())
  return (list(select_row, beta))
}


# Function 2: A function to select representative samples from poisson regression model with backtrack
fam_poi <- function(DATA,lambda,p1_low,select_row) 
{
  # DATA: the original dataset
  # lambda: a threshold to judge whether a sample is representative 
  # p1_low: a percentage, if the percentage of final representative samples < p1_low,
        # exclude these final samples from the original dataset, then restart the sample selection
  # select_row: the order number of rows selected from the original dataset "DATA"
  
  N=nrow(DATA)
  p= ncol(DATA)
  beta = matrix(rep(0,p),p,1) 
  p1_hat = 0

  # if no given samples, choose half of dataset randomly
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
  
  # if the updated data remain the same, stop2 = 1, end the inner loop 
  stop2=0
  
  # max times of iteration 
  max_iter=100
  i=0
  
  while (stop1==0) 
  {
    while (stop2 == 0 && i <= max_iter)
    {

      # backtrack
      X = x_step[,1:(p-1)]
      Y = x_step[,p]
      X_1 = as.matrix(cbind(rep(1,nrow(x_step)),X))
      
      # most important parameter, depending on the gradient
      step = 5*10^(-7)
      maxiter = 100
      
      a = 0.2
      b = 0.9
      
      for (iter in 1:maxiter)
      {
        beta_copy = beta
        
        # compute the objective function
        q = sum(exp(X_1 %*% beta)) - Y %*% (X_1 %*% beta) 
        
        # compute gradient and update beta
        grad = t(X_1) %*% (exp(X_1 %*% beta)-Y) 
        beta_new= beta - step* grad
        
        # updated objective function
        q_new = sum(exp(X_1 %*% beta_new)) - Y %*% (X_1 %*% beta_new) 
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
      
      # choose the rows, whoses deviance is smaller than the threshold
      DATA1= cbind(rep(1,N),DATA)

      y_pre = exp((DATA1[,-(p+1)])%*% beta)
      devia = as.vector(2*(y_pre-DATA[,p])+2*DATA[,p]*log(DATA[,p]/y_pre))
      
      num =which(devia < lambda, arr.ind = T)
  
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
        # print(dim(x_step))
        select_row=num
        x_step=DATA[select_row,]
        }
      
      i=i+1
    }
    
    # if selected samples' percentage > percentage set up, or select_row = 0, stop the outer while loop
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
  assign("final_beta_bt",beta,envir = globalenv())
  return (list(select_row, beta))
}

# Function 3: create a function to create static plot for each step
auto_plot_poi <- function(DATA,lambda,p1_low,select_row) 
{
  # DATA: the original dataset
  # lambda: a threshold to judge whether a sample is representative 
  # p1_low: a percentage, if the percentage of final representative samples < p1_low,
        # exclude these final samples from the original dataset, then restart the sample selection
  # select_row: the order number of rows selected from the original dataset "DATA"
  
  N=nrow(DATA)
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
      # prepare for the plot after each update
      dataframe=as.data.frame(x_step)
      colnames(dataframe) <- c("x","y","type")
      dataframe$type=as.factor(dataframe$type)
      
      # compute beta from the selected samples 
      data <- data.frame(y=x_step[,2], x=x_step[,1])             
      beta <- glm(x_step[,2] ~ x_step[,1], data=data, family=poisson(link="log"))
      
      # compute the probability for each row of selected samples
      inter=beta$coefficients[1]
      slope=beta$coefficients[2]
      beta$model$fitted <- predict(beta, type = "response")
      
      library(ggplot2)
      cols <- c("1" = "red", "0" = "blue")
      
      # diaplay plot
      print(ggplot(dataframe, 
                   aes(x=x,  
                       y=y,  
                       fill =type))  +   
              xlim(-5,5)+ylim(-4,3100)+
              geom_point(alpha = 1,shape=21, size =3)+ 
              scale_fill_manual(values = cols)+
              geom_smooth(aes(y=fitted(beta)))+
              labs(title= paste("step" ,as.character(i)),envir = globalenv()))

      fitt=fitted(beta)
      
      # compute predicted y 
      DATA1= cbind(rep(1,N),DATA[,1:(p-2)])
      y_pre = exp(DATA1 %*% beta$coefficients)
      
      # select the rows, whose deviance is smaller than the threshold
      devia = as.vector(2*(y_pre-DATA[,p-1])+2*DATA[,p-1]*log(DATA[,p-1]/y_pre))
      num = which(devia<lambda, arr.ind = T) 
      
      # update dataset
      x_stepnew=DATA[num,]
      
      # if every row num of new dataset = row num of data before updating, stop the inner while
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
      assign(paste("fit_value", i+1, sep = ""),
             fitted(beta),envir = globalenv() )
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
  return (list(select_row, beta)) 
}


# Function 4 : gif 

dynamic_plot_poi <- function(dataframe_list,fit_list) 
{
  df_list1=dataframe_list
  fit_list1=fit_list
  cols <- c("1" = "red", "0" = "blue")
  library(animation)
  library(ggplot2)
  i=0
  
  saveGIF({
    for(i in seq_along(df_list1)){
     # compute the probability of each row of selected samples
      a <- ggplot(df_list1[[i]], 
                  aes(x = x,  y = y, fill = type))  +   
        xlim(-5,5)+ylim(-4,1100)+
        geom_point() +
        geom_point(alpha = 1,shape=21, size =3)+ 
        scale_fill_manual(values = cols) +
        geom_smooth(aes(y=fit_list1[[i]]))+
        labs(title= paste("step" ,as.character(i)),envir = globalenv())
      print(a)
    }
  }, interval = 1, movie.name="test.gif")
}

# Function 5: Create a Autoregressive correlation sigma
sigma_auto <- function(p,rou) 
{
  sig= matrix(1,p,p)
  for (i in 1:p-1)
  {
    a=(1+i):p
    b=1:(p-i)
    w1=(b-1)*p+a
    w2=(a-1)*p+b   
    sig[w1]=rou^(i)
    sig[w2]=rou^(i)
  }
  
  return(sig)
}


# Function 6 : Create a compound symmetry sigma
sigma_compound <- function(p,rou) 
{
  sig1= diag(rou,p,p)
  sig2= matrix(rou,p,p)
  sig=sig1+sig2
  return(sig)
}


# Function 7 : Backtrack
backtrack_poi<-function(data)
{
  n = nrow(data)
  p = ncol(data)
  X = data[,1:(p-1)]
  Y = data[,p]
  X_1 = cbind(rep(1,n),X)
  theta = matrix(rep(0,p),p,1)  
  
  # most important parameter, depending on the gradient
  step = 9*10^(-7)
  maxiter = 100
  a = 0.2
  b = 0.9
  
  for (iter in 1:maxiter)
  {
    theta_copy = theta
    # minimize the objective function
    q = sum(exp(X_1 %*% theta)) - Y %*% (X_1 %*% theta) 
    
    grad = t(X_1) %*% (exp(X_1 %*% theta)-Y) 
    theta_new= theta - step* grad
    
    # updated objective
    q_new = sum(exp(X_1 %*% theta_new)) - Y %*% (X_1 %*% theta_new) 
    q_object = q - a*step*t(grad) %*% grad
    
    if (q_new > q_object)
    {
      step = step*b
    }
    else
    {
      theta = theta_new
      break
    }
  }
  
  return(list(theta,iter))
}

  

 
