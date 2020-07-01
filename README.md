# Representative samples selection in an implicit mixture model with the Approximation Maximization algorithm

Technological innovations have made a profound impact on knowledge discovery. Extracting useful samples from massive dataset is essential in many modern
scientific areas. In this project, we developed an Approximation Maximization (AM) algorithm to select representative samples in an implicit mixture model. 

# Methodology 
The detailed theoretical derivation of AM algorithm can be referred to the project report.

![plot1](https://user-images.githubusercontent.com/64850893/86276549-16837a00-bba3-11ea-8430-abe074154938.jpg)

FAM is a more general version of AM. 
![plot4](https://user-images.githubusercontent.com/64850893/86276955-da044e00-bba3-11ea-801e-a5ec2cb2263d.jpg)

# Experiment settings
I experiment our AM algorithm under four classes of implicit mixture models, where the explicit models are Gaussian distribution, linear regression, logistic regression and Poisson regression, respectively. Regarding to implicit component models, I designe some totally different models. Based on the above four models, I compare the estimating and selecting performance of our AM and FAM with three methods. In the first method, all samples are directly assumed to come from the explicit model, and the parameters of interest are estimated by maximum likelihood estimation. The first method is called "MLE". The second method considers using EM algorithm based on mixture models to conduct estimation and representative samples selection, which is denoted as "EM". In the third method, the whole dataset will be clustered into two subsets by K-means, and the parameters will be estimated by MLE under the explicit model assumption. This estimating and selection procedure is referred to as "K-means".

In different models, I need to generate some multiple dimensional variables, such as y in implicit multiple Gaussian mixture model or the covariate vector x in implicit general linear mixture model. In this report, I consider three different covariance structures for y and x. The first structure is that the variables are independent with each other. The second structure allows that variables are auto-regressive correlated, which means the correlation of adjacent variables are higher than that of the distant ones. The final structure is compound symmetry, where the correlation between two variables remains the same.

I use the deviation (DEV) of true parameter to evaluate the estimation performance. In addition, I evaluate the performance of representative samples selection in terms of averaging positive selection rate (PSR) and false discovery rate (FDR). Further to that, the final number (FN) of selected representative samples by each algorithm is also given, which is a more intuitive index to show the selection performance.

# Case 1.1 Implicit normal mixture model (2 variables)

![plot2](https://user-images.githubusercontent.com/64850893/86276183-73326500-bba2-11ea-8fd8-c126cde82109.jpg)

From the above plot, it is obvious that more and more red points are selected into representative samples set as iteration grows. On the other side, the blue noise samples are ruled out gradually by our AM algorithm, which means AM algorithm can effectively conduct the representative samples selection.


# Case 1.2 Implicit normal mixture model (multiple variables)

![plot3](https://user-images.githubusercontent.com/64850893/86276657-3fa40a80-bba3-11ea-85c2-fddd09cad604.jpg)

From the above table, it is noticeable that more complex covariance structure in this experiment does not lead to obvious reduction of accuracy for AM and FAM. By comparing AM with FAM, AM is computationally more efficient than FAM, which is due to a close form of true parameter in this model.


# Case 2.1 Implicit linear mixture model (1 variable)

![plot5](https://user-images.githubusercontent.com/64850893/86277218-54cd6900-bba4-11ea-867d-503cb6a8a3e6.jpg)

According to the figure, AM algorithm successfully selects almost all the red representative samples and few blue noise samples after 7 iterations.


# Case 2.2 Implicit linear mixture model (multiple variables)

![plot6](https://user-images.githubusercontent.com/64850893/86277385-9827d780-bba4-11ea-86e2-7e4b22c3286f.jpg)

In this experiment, the performance of all 5 methods is similar in 3 setups. AM obviously outperforms other methods in terms of DEV and selection indexes. Due to the similar reasons, "MLE", "EM" and "K-means" fail to obtain accurate estimate and representative samples set. As for FAM, although it remains good selection performance in this experiment, it is noteworthy that its estimation accuracy is lower than that in previous experiments.


# Case 3.1 Implicit logistic mixture model (1 variable)

![plot7](https://user-images.githubusercontent.com/64850893/86277737-39af2900-bba5-11ea-80f6-e650a45a2c47.jpg)

In summary, the AM algorithm successfully selects a majority of the representative samples with only few noise samples after 4 iterations.


# Case 3.2 Implicit logistic mixture model (multiples variables)

![plot8](https://user-images.githubusercontent.com/64850893/86277895-7ed35b00-bba5-11ea-9c4f-0f28972eb7d1.jpg)

The estimation and selection results of five methods are presented in the above table. To start with, "EM", AM and FAM perform well in selecting representative samples in terms of relatively high PSR and low FDR. From the perspective of DEV, "EM" performs more accurate than AM and FAM. The result of "EM" is different from that in previous models, which is due to that the implicit model assumption of "EM" matches the true model setup. Even so, our AM and FAM still obtain satisfactory results without the assumption on implicit component model. It is noteworthy that, unlike the experiments in Gaussian or linear mixture models, the PSR or FDR in this experiment does not reach the perfect level 100% or 0.


# Case 4.1 Implicit Poisson mixture model (1 variable)

![plot9](https://user-images.githubusercontent.com/64850893/86278141-e1c4f200-bba5-11ea-895e-17d92ee26e12.jpg)
In this case, most of the representative samples with few noise samples are selected by AM after 5 iterations.



# Case 4.2 Implicit Poisson mixture model (multiple variables)
![plot10](https://user-images.githubusercontent.com/64850893/86278315-310b2280-bba6-11ea-9a2d-8ebce871ea62.jpg)
In this experiment, AM and FAM perform better than other methods when estimation and selection are together considered. It is noteworthy that, by the above data generation methods, some noise samples are close to true Poisson regression model, which means a small amount of noise samples will be chosen as representative samples.




# Conclusion
In this report, we developed an Approximation Maximization (AM) algorithm to select representative samples in an implicit mixture model. This approach shows potential to improve the accuracy of estimated parameters of the explicit model, which helps us extract the representative samples out of a complicate dataset. Beside, we further designed a more general first-order AM algorithm. The corresponding estimating and selecting procedure is compatible with a broad range of explicit models. I implemented experiments in 4 cases, where the explicit models are Gaussian distribution, linear regression, logistics regression and Poisson regression, respectively. Under 3 different correlation structures among variables, the AM algorithm is robust and outperforms other methodS in terms of estimation accuracy and selection consistency in most experiments.


# Future work

In our future investigation, we will attempt to study the specific influence of initial parameter input in our algorithm. Besides, we will design a general and
data-driven threshold rule. In addition, the in-depth theoretical guarantees of AM is another important point in our future work.



