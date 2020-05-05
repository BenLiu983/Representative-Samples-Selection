# Function 1: create a function to select samples
am_normal <- function(DATA,sig,lambda,p1_low,select_row) 
{
  # DATA: the original dataset
  # sig: the covariance matrix of the dataset, if it is unknown, input matrix(0,1,1)
  # lambda: a threshold to judge whether a sample is representative 
  # p1_low: a percentage, if the percentage of final representative samples < p1_low,
       # exclude these final samples from the original dataset, then restart the sample selection
  # select_row: the order number of rows selected from the original dataset "DATA"
  
  N= nrow(DATA)
  p= ncol(DATA)
  p1_hat = 0
  
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
  
  # if the updated selected samples remain the same, stop2 = 1, end the inner loop 
  stop2=0
  
  # max times of iteration 
  max_iter=100
  i=0
  
  while (stop1==0) 
  {
    while (stop2 == 0 && i <= max_iter)
    {
      # mu
      mu_hat = colMeans(x_step)  
      
      # sigma
      if (norm(sig)==0)
      {
        number_sig = nrow(x_step) 
        mu_matrix = matrix(rep(mu_hat,number_sig),number_sig,p,byrow=TRUE)   
        sigma_hat = 1/number_sig * t(x_step-mu_matrix) %*% (x_step-mu_matrix)
      }
      
      else
      {
        sigma_hat = sig
      }
    
      # select the order number of samples, whose deviance is smaller than the threshold
      
      mu_rep= matrix(rep(mu_hat,N),N, p,byrow=TRUE)
                        
      devia  =  (DATA-mu_rep) %*% solve(sigma_hat) %*% t(DATA-mu_rep)
      devia_vec = diag(devia)
      
      num =which(devia_vec<lambda, arr.ind = T)
      
      # update dataset
      x_stepnew=DATA[num,]
      
      # if the each row umn of new dataset = row num of dataset before, stop the inner while
      if ((length(select_row)==length(num)) && all(select_row == num))
        
      {  
        p1_hat =nrow(x_stepnew)/N
        stop2=1
      }
      
      else
      {
        # print(dim(x_step))
        select_row=num
        x_step=DATA[select_row,]
      }

      i=i+1
    }
    
    # if percentage of representative samples > percentage set up, stop the whole outer while loop
    if ((p1_hat>p1_low) || (length(num)==0))
    {stop1=1}
    
    # else, exclude the samples selected, use the rest samples to start over
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
  
  # save some final indexes
  assign("final_row",
         select_row,envir = globalenv())
  assign("final_mu",
         mu_hat,envir = globalenv())
  assign("final_sigma",
         sigma_hat,envir = globalenv())
  
  return (list(select_row, mu_hat,sigma_hat)) 
}


# Function 2: create a function to select samples, backtrack
fam_normal <- function(DATA,sig,lambda,p1_low,select_row) 
{
  # DATA: the original dataset
  # sig: the covariance matrix of the dataset, if it is unknown, input matrix(0,1,1)
  # lambda: a threshold to judge whether a sample is representative 
  # p1_low: a percentage, if the percentage of final representative samples < p1_low,
          # exclude these final samples from the original dataset, then restart the sample selection
  # select_row: the order number of rows selected from the original dataset "DATA"
  
  N=nrow(DATA)
  p= ncol(DATA)
  sigma_hat = sig 
  theta = matrix(rep(0,p),p,1)  
  
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
  
  # if the percentage of final representative samples > p1_low, stop1 = 1, end the outer loop
  stop1=0
  
  # if the updated samples remain the same, stop2 = 1, end the inner loop 
  stop2=0
  
  # max times of iteration 
  max_iter=100
  i=0
  
  while (stop1==0) 
  {
    while (stop2 == 0 && i <= max_iter)
    {
      X = as.matrix(x_step)
      
      # most important parameter, depend on gradient
      step = 2*10^(-3)
      
      maxiter = 100
      a = 0.3
      b = 0.7
      
      for (iter in 1:maxiter)
      {
        # calculate the objective function
        theta_copy = theta
        theta_matrix = matrix(rep(theta,nrow(x_step)),nrow(x_step),p,byrow = TRUE)
        q = 1/2 * sum(diag((X - theta_matrix) %*% solve(sigma_hat) %*% t(X - theta_matrix)))
        
        # calculate the gradient and update theta
        grad = solve(sigma_hat) %*% as.matrix(colSums(theta_matrix - X))
        theta_new = theta - step* grad
        
        # calculate new objective function
        theta_matrix_new = matrix(rep(theta_new,nrow(x_step)),nrow(x_step),p,byrow = TRUE)
        q_new = 1/2 * sum(diag((X - theta_matrix_new) %*% solve(sigma_hat) %*% t(X - theta_matrix_new)))
        
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
      
      # select the order number of samples, whose deviance is smaller than a threshold

      mu_rep= matrix(rep(theta,N),N, p,byrow=TRUE)
      
      devia = (DATA-mu_rep) %*% solve(sigma_hat) %*% t(DATA-mu_rep)
      devia_vec = diag(devia)
      
      num = which(devia_vec<lambda, arr.ind = T)
      
      # update data
      x_stepnew=DATA[num,]
      
      # if the each row num of new samples = row num of samples before, stop the inner while
      if ((length(select_row)==length(num)) && all(select_row == num))
        
      {  
        p1_hat =nrow(x_stepnew)/N
        stop2=1
      }
      
      else
      {
        # print(dim(x_step))
        select_row=num
        x_step=DATA[select_row,]
      }
      
      i=i+1
    }
    
    # if percentage of representative samples > percentage set up, stop the whole outer while loop
    if ((p1_hat>p1_low) || (length(num)==0))
    {stop1=1}
    
    # else, exclude the data selected, use the rest dataset to start over
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
  assign("final_row_bt",
         select_row,envir = globalenv())
  assign("final_mu_bt",
         theta,envir = globalenv())
  
  return (list(select_row, theta)) 
}


# Function 3: create a function to create static plot for each step
am_plot <- function(DATA,sig,lambda,p1_low,select_row) 
{
  # DATA: the original dataset
  # sig: the covariance matrix of the dataset, if it is unknown, input matrix(0,1,1)
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
    colnames(dataframe) <- c("x1","x2","type")
    dataframe$type=as.factor(dataframe$type)
  }
  
  else
  {
    x_step = DATA[select_row,] 
    dataframe = as.data.frame(x_step)
    colnames(dataframe) <- c("x1","x2","type")
    dataframe$type=as.factor(dataframe$type)
  }
  
  # if the percentage of final representative samples > p1_low, stop1 = 1, end the outer loop
  stop1=0
  
  # if the updated samples remain the same, stop2 = 1, end the inner loop 
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
      colnames(dataframe) <- c("x1","x2","type")
      dataframe$type=as.factor(dataframe$type)
      
      cols <- c("1" = "red", "0" = "blue")
      library(ggplot2)
      
      # diaplay plot
      print(ggplot(dataframe, 
                   aes(x=x1,  
                       y=x2,  
                       fill =type))  +   
              geom_point(alpha = 1,shape=21, size =3)+ 
              scale_fill_manual(values = cols)+
              xlim(-20,10)+ylim(-10,20)+labs(title= paste("step" ,as.character(i))),envir = globalenv())
      
      # mu
      mu_hat = colMeans(x_step[,-p])
      
      # sigma
      if (norm(sig)==0)
      {
        number_sig=nrow(x_step)
        mu_matrix = matrix(rep(mu_hat,number_sig),number_sig,p-1,byrow=TRUE)
        sigma_hat = 1/number_sig * t(x_step[,-p]-mu_matrix) %*% (x_step[,-p]-mu_matrix)
      }
      
      else
      {
        sigma_hat = sig
      }
      
      # select the rows, whose deviance is smaller than the threshold
      mu_rep = matrix(rep(mu_hat,N),N, p-1,byrow=TRUE)
      devia  =  (DATA[,-p]-mu_rep) %*% solve(sigma_hat) %*% t(DATA[,-p]-mu_rep)
      devia_vec = diag(devia)
      num = which(devia_vec<lambda, arr.ind = T)
      
      # update dataset
      x_stepnew=DATA[num,]
      
      # if the row num of new dataset = row num of data before updating, stop the inner while
      if   ((nrow(x_step)==nrow(x_stepnew)) && all(x_step == x_stepnew) )
      {  
        p1_hat = nrow(x_stepnew)/N
        stop2=1
      }
      
      else
      {
        #print(dim(x_step))
        x_step=x_stepnew
      }
      
      # save dataframe of each step for gif, "df" + 100 is for the dynamic gif
      assign(paste("df", i+100, sep = ""),
             dataframe,envir = globalenv() )
      assign(paste("mu_hat", i+1, sep = ""),
             mu_hat,envir = globalenv() )
      assign(paste("sigma_hat", i+1, sep = ""),
             sigma_hat,envir = globalenv() )
      assign(paste("deviation", i+1, sep = ""),
             norm(mu_hat-c(2,6),type="2")  / norm(c(2,6), type="2") ,envir = globalenv() )
      i=i+1
    }
    
    # if selected representative samples' percentage > percentage set up, stop the whole outer while
    if (p1_hat>p1_low)
    {stop1=1}
    
    # else, exclude the selected samples, use the rest dataset to start over
    else
    {
      x_step= DATA[-num,]
      dataframe=as.data.frame(x_step)
      N1=nrow(x_step)
      x_step= x_step[sample(nrow(x_step), N1/2), ]
      dataframe=as.data.frame(x_step)
      colnames(dataframe) <- c("x1","x2","type")
      dataframe$type=as.factor(dataframe$type)
      stop2=0
      i=0
    }
  }
}


# Function 4 : gif 
dynamic_plot <- function(dataframe_list) 
{
  df_list1=dataframe_list
  cols <- c("1" = "red", "0" = "blue")
  library(animation)
  library(ggplot2)
  i=0
  saveGIF({
    for(i in seq_along(df_list1)){
      a <- ggplot(df_list1[[i]], 
                  aes(x = x1,  y = x2, fill = type))  +   
        geom_point(alpha = 1,shape = 21, size = 3)+ 
        scale_fill_manual(values = cols) +
        xlim(-20, 10) + ylim(-10, 20)+labs(title= paste("step" ,as.character(i)),envir = globalenv())
      print(a)
    }
  }, interval = 1, movie.name="test.gif")
}


# Function 5 : Create a  compound symmetry sigma
sigma_compound <- function(p,rou) 
{
  sig1= diag(rou,p,p)
  sig2= matrix(rou,p,p)
  sig=sig1+sig2
  return(sig)
}


# Function 6: Create a Autoregressive correlation sigma
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







