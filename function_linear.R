# Function 1: create a function to select representative samples from linear model
am_linear <- function(DATA,sig,lambda,p1_low,select_row) 
{
  # DATA: the original dataset
  # sig: the covariance matrix of the covariates, if it is unknown, input 0
  # lambda: a threshold to judge whether a sample is representative 
  # p1_low: a percentage, if the percentage of final representative samples < p1_low,
         # exclude these final samples from the original dataset, then restart the sample selection
  # select_row: the order number of rows selected from the original dataset "DATA"
 
  N=nrow(DATA)
  p= ncol(DATA)
  DATA_1 = as.matrix(cbind(rep(1,N),DATA[,-p]))
  
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
      n1=nrow(x_step)
      
      # add a column 1 to compute the intersept, the last column of x_step is "y" 
      x_step_1 = as.matrix(cbind(rep(1,n1),x_step[,-p]))

      beta = (solve(t(x_step_1) %*% x_step_1))%*% t(x_step_1) %*% x_step[,p]

      # compute expected "y" for each row of orginal dataset
      x_step_all = as.matrix(cbind(rep(1,N),DATA[,1:(p-1)]))
      mu_all =  x_step_all %*% beta 
      
      # sigma
      if (sig==0)
      {
        sigma_hat=var(x_step[,p] - x_step_1 %*% beta)*(n1-1)/n1
      }
      
      
      else
      {
        sigma_hat = sig
      }
      
      # select the rows, whose deviance is smaller than a threshold
      # den = pi_star * dnorm(x=DATA[,p], mean=mu_all,sd=sqrt(sigma_hat))
      devia =  (DATA[,p]-mu_all)  %*% sigma_hat %*% t(DATA[,p]-mu_all)
      devia_vec = diag(devia)
      num =which(devia_vec < lambda, arr.ind = T)
      
      # update data
      x_stepnew=DATA[num,]
      
      # if the every row num of new dataset = row num of data before updating, stop the inner while
      if (  (length(select_row) == length(num))    && all(select_row == num))  

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
    if (p1_hat>p1_low)
    {stop1=1}
    
    # else, exclude the samples selected, use the rest dataset to start over
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
  assign("final_row",
         select_row,envir = globalenv())
  assign("final_beta",
         beta,envir = globalenv())
  assign("final_sigma",
         sigma_hat,envir = globalenv())
  
  return (list(select_row, beta)) 
}


# Function 2: create a function to select representative samples from linear model with backtrack
fam_linear <- function(DATA,sig,lambda,p1_low,select_row) 
{
  # DATA: the original dataset
  # sig: the covariance matrix of the covariates, if it is unknown, input 0
  # lambda: a threshold to judge whether a sample is representative 
  # p1_low: a percentage, if the percentage of final representative samples < p1_low,
        # exclude these final samples from the original dataset, then restart the sample selection
  # select_row: the order number of rows selected from the original dataset "DATA"
  
  N=nrow(DATA)
  p= ncol(DATA)
  p1_hat = 0
  
  DATA_1 = as.matrix(cbind(rep(1,N),DATA[,-p]))
  #beta = (solve(t(DATA_1) %*% DATA_1))%*% t(DATA_1) %*% DATA[,p]
  betav = rep(0,200)
  beta = rep(0,p)
  
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
  
  # if the updated data remain the same, stop2 = 1, end the inner loop 
  stop2=0
  
  # max times of iteration 
  max_iter=200
  i=0
  
  while (stop1==0) 
  {
    while (stop2 == 0 && i <= max_iter)
    {
      n1=nrow(x_step)
      

      # sigma
      if (sig==0)
      {
        sigma=var(x_step[,p] - x_step_1 %*% beta)*(n1-1)/n1
      }
      
      else
      {
        sigma = sig
      }
      
      # backtrack
      X = x_step[,1:(p-1)]
      Y = x_step[,p]
      X_1 = as.matrix(cbind(rep(1,nrow(x_step)),X))
    
      #  important 
      step = 1*10^(-4)
      
      maxiter2 = 100
      a = 0.4
      b = 0.7
      
      for (iter in 1:maxiter2)
      {
        beta_copy = beta
        
        # calculate the objective function
        q = 1/(2*sigma^2) * t(Y - (X_1 %*% beta)) %*% (Y - (X_1 %*% beta))
        
        # calculate the gradient and update theta
        grad = 1/sigma^2 * t(X_1) %*% (X_1 %*% beta-Y) 
        beta_new = beta - step* grad
        
        # calculate new objective function
        q_new = 1/(2*sigma^2) * t(Y - (X_1 %*% beta_new)) %*% (Y - (X_1 %*% beta_new))
        
        q_object = q - a*step*t(grad) %*% grad   # step
        
        if (q_new > q_object)    # 
        {
          step = step*b
        }
        
        else
        {
          beta = beta_new
          break
        }
      }

      # compute expected "y" for each row of orginal dataset
      x_step_all = as.matrix(cbind(rep(1,N),x_step00[,1:(p-1)]))
      mu_all =  x_step_all %*% beta 
      
      # select the rows whose denviance is smaller than the threshold
      
      devia =  (DATA[,p]-mu_all)  %*% sigma %*% t(DATA[,p]-mu_all)
      devia_vec = diag(devia)
      num =which(devia_vec < lambda, arr.ind = T)
      
      # update dataset
      x_stepnew=DATA[num,]
      
      # if the every row num of new dataset = row num of data before updating, stop the inner while
      if (  (length(select_row) == length(num))    && all(select_row == num)) 
        
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
  assign("final_row_bt",
         select_row,envir = globalenv())
  assign("final_beta_bt",
         beta,envir = globalenv())
  
  return (list(select_row, beta)) 
}


# Function 3: create a function to create static plot for each step
plot_linear <- function(DATA,sig,lambda,p1_low,select_row) 
{
  # DATA: the original dataset
  # sig: the covariance matrix of the covariates, if it is unknown, input 0
  # lambda: a threshold to judge whether a sample is representative 
  # p1_low: a percentage, if the percentage of final representative samples < p1_low,
        # exclude these final samples from the original dataset, then restart the sample selection
  # select_row: the order number of rows selected from the original dataset "DATA"
  
  N=nrow(DATA)
  p= ncol(DATA)
  
  # if no given sample, choose half of dataset randomly
  if (sum(select_row)==0)
  {
    select_row = sample(N, floor(N/2))
    x_step = DATA[select_row,] 
    dataframe = as.data.frame(x_step)
    
    # prepare for the plot
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
      # prepare for the plot after each update
      dataframe=as.data.frame(x_step)
      colnames(dataframe) <- c("x","y","type")
      dataframe$type=as.factor(dataframe$type)
      
      n1=nrow(x_step)
      
      # add a column 1 to compute the intercept, and exclude the last 2 columns ("y" and "type")
      x_step_1 = as.matrix(cbind(rep(1,n1),x_step[,-c((p-1):p)]))
      beta = (solve(t(x_step_1) %*% x_step_1))%*% t(x_step_1) %*% x_step[,p-1]
      
      cols <- c("1" = "red", "0" = "blue")
      library(ggplot2)
      
      # save plot 
      assign(paste("plot", i+1, sep = ""),
             ggplot(dataframe, 
                    aes(x=x,  
                        y=y,  
                        fill =type))  +   
               geom_point(alpha = 1,shape=21, size =3)+ 
               scale_fill_manual(values = cols)+
               xlim(-50,50)+ylim(-100,500)+geom_abline(intercept = beta[1,1], slope = beta[2,1]) +
               labs(title= paste("step" ,as.character(i))),envir = globalenv())
      
      # diaplay plot
      print(ggplot(dataframe, 
                   aes(x=x,  
                       y=y,  
                       fill =type))  +   
              geom_point(alpha = 1,shape=21, size =3)+ 
              scale_fill_manual(values = cols)+
              xlim(-50,50)+ylim(-100,500)+geom_abline(intercept = beta[1,1], slope = beta[2,1]) +
              labs(title= paste("step" ,as.character(i))),envir = globalenv())
      
      # compute expected "y" for each row of original data
      x_step_all = as.matrix(cbind(rep(1,N),x_step01[,-c(2:3)]))
      mu_all =  x_step_all %*% beta 
      
      # sigma
      if (sig==0)
      {
        sigma_hat=var(x_step[,p-1] - x_step_1 %*% beta)*(n1-1)/n1
      }
      
      else
      {
        sigma_hat = sig
      }
      
      # select the rows, whose deviance is smaller than the threshold
      devia =  (DATA[,p-1]-mu_all)  %*% sigma_hat %*% t(DATA[,p-1]-mu_all)
      devia_vec = diag(devia)
      num =which(devia_vec < lambda, arr.ind = T)
      
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
      assign(paste("df", i+100, sep = ""),
             dataframe,envir = globalenv() )
      assign(paste("beta", i+100, sep = ""),
             beta,envir = globalenv() )
      
      i=i+1
    }
    
    # if selected representative samples' percentage > percentage set up, stop the whole outer while
    if (p1_hat>p1_low)
    { 
      stop1=1
    }
    
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
}

# Function 4 : gif 
dynamic_plot_linear <- function(dataframe_list,beta_list) 
{
  df_list1=dataframe_list
  b_list1=beta_list
  cols <- c("1" = "red", "0" = "blue")
  library(animation)
  library(ggplot2)
  i=0
  saveGIF({
    for(i in seq_along(df_list1)){
      a <- ggplot(df_list1[[i]], 
                  aes(x = x,  y = y, fill = type))  +   
        geom_point(alpha = 1,shape = 21, size = 3)+ 
        scale_fill_manual(values = cols) +
        xlim(-80, 80) + ylim(-140, 400)+geom_abline(intercept =  b_list1[[i]][1], slope = b_list1[[i]][2]) +
        labs(title= paste("step" ,as.character(i)),envir = globalenv())
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


