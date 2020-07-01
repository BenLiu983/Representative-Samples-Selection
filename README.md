# Representative samples selection in an implicit mixture model with the Approximation Maximization algorithm

Technological innovations have made a profound impact on knowledge discovery. Extracting useful samples from massive dataset is essential in many modern
scientific areas. In this project, we developed an Approximation Maximization (AM) algorithm to select representative samples in an implicit mixture model. 

# Methodology 
The detailed theoretical derivation of AM algorithm can be refer to the project report.

![plot1](https://user-images.githubusercontent.com/64850893/86276549-16837a00-bba3-11ea-8430-abe074154938.jpg)

FAM is a more general version of AM. 
![plot4](https://user-images.githubusercontent.com/64850893/86276955-da044e00-bba3-11ea-801e-a5ec2cb2263d.jpg)

# Experiment settings
In this section, I experiment our AM algorithm under four classes of implicit mixture models, where the explicit models are Gaussian distribution, linear regression, logistic regression and Poisson regression, respectively. The last three models can be summarized as implicit general linear mixture models. Regarding to implicit component models, I design some totally different models. Based on the above four models, I compare the estimating and selecting performance of our AM and FAM with three methods. In the first method, all samples are directly assumed to come from the explicit model, and the parameters of interest are estimated by maximum likelihood estimation. The first method just focuses on parameter estimation and is called "MLE". The second method considers using EM algorithm based on mixture models to conduct estimation and representative samples selection. For ease of implementation, the implicit model is assumed to have same distribution function as the explicit model, but their parameters are different. Among all subgroups after clustering, I choose the subgroup with largest mean likelihood as the representative samples set. This method is denoted as "EM". In the third method, the whole dataset will be clustered into two subsets by K-means, and the parameters will be estimated by MLE under the explicit model assumption. The representative samples set is determined by the same way as "EM". This estimating and selection procedure is referred to as "K-means".

In different models, I need to generate some multiple dimensional variables,such as y in implicit multiple Gaussian mixture model or the covariate vector x in implicit general linear mixture model. In this report, I consider three different covariance structures for y and x. The first structure is that the variables are independent with each other. The second structure allows that variables are auto-regressive correlated, which means the correlation of adjacent variables are higher than that of the distant ones. The final structure is compound symmetry, where the correlation between two variables remains the same.

I use the deviation (DEV) of true parameter to evaluate the estimation performance. In addition, I evaluate the performance of representative samples selection in terms of averaging positive selection rate (PSR) and false discovery rate (FDR). Further to that, the final number (FN) of selected representative samples by each algorithm is also given, which is a more intuitive index to show the selection performance.

# Case 1.1 Implicit normal mixture model (2 variables)

![plot2](https://user-images.githubusercontent.com/64850893/86276183-73326500-bba2-11ea-8fd8-c126cde82109.jpg)

From the above plot, it is obvious that more and more red points are selected into representative samples set as iteration grows. On the other side, the blue noise samples are ruled out gradually by our AM algorithm, which means AM algorithm can effectively conduct the representative samples selection.


# Case 1.2 Implicit normal mixture model (multiple variables)

![plot3](https://user-images.githubusercontent.com/64850893/86276657-3fa40a80-bba3-11ea-85c2-fddd09cad604.jpg)

From the above table, it is noticeable that more complex covariance structure in this experiment does not lead to obvious reduction of accuracy for AM and FAM. By comparing AM with FAM, AM is computationally more efficient than FAM, which is due to a close form of true parameter in this model.


# Case 2.1 Implicit linear mixture model (2 variables)

![plot5](https://user-images.githubusercontent.com/64850893/86277218-54cd6900-bba4-11ea-867d-503cb6a8a3e6.jpg)

According to the figure, AM algorithm successfully selects almost all the red representative samples and few blue noise samples after 7 iterations.


# Case 2.2 Implicit linear mixture model (multiple variables)

![plot6](https://user-images.githubusercontent.com/64850893/86277385-9827d780-bba4-11ea-86e2-7e4b22c3286f.jpg)

In this experiment, the performance of all 5 methods is similar in 3 setups. AM obviously outperforms other methods in terms of DEV and selection indexes. Due to the similar reasons, "MLE", "EM" and "K-means" fail to obtain accurate estimate and representative samples set. As for FAM, although it remains good selection performance in this experiment, it is noteworthy that its estimation accuracy is lower than that in previous experiments.


























This approach shows potential to improve the accuracy of estimated parameters of the explicit model, which helps us extract the representative samples out of a complicate
dataset. Experiments are conducted in 4 main cases, where the explicit models are Gaussian distribution, linear regression, logistics regression and Poisson regression, respectively. Under 3 different correlation structures among variables, the AM algorithm is robust and outperforms other method in terms of estimation accuracy and selection
consistency in most experiments.



