# Representative Samples Selection in An Implicit Mixture Model with The Approximation Maximization Algorithm

Technological innovations have made a profound impact on knowledge discovery. Extracting useful samples from massive dataset is essential in many modern
scientific areas. In this project, I developed an Approximation Maximization (AM) algorithm to select representative samples in an implicit mixture model. 

# Methodology 
The detailed theoretical derivation of AM algorithm can be referred as follows.

![plot1](https://user-images.githubusercontent.com/64850893/86276549-16837a00-bba3-11ea-8430-abe074154938.jpg)

FAM is a more general version of AM. 
![plot4](https://user-images.githubusercontent.com/64850893/86276955-da044e00-bba3-11ea-801e-a5ec2cb2263d.jpg)

# Experiment settings

Data: Four classes of implicit mixture datasets, where the explicit models are normal distribution, linear regression, logistic regression, and Poisson regression. 

Models: AM, FAM (a more general version of AM), EM, MLE, and K-means.

Covariance structures: Independent with each other (S1), auto-regressive correlated (S2), and compound symmetry(S3).

Metrics:  Deviation of true parameters (DEV), positive selection rate (PSR), false discovery rate (FDR),  the final number of selected representative samples (FN), and computation time (Time).


# Case 1.1 Implicit normal mixture model (2 variables)

![plot2](https://user-images.githubusercontent.com/64850893/86276183-73326500-bba2-11ea-8fd8-c126cde82109.jpg)

In this case, more and more red points were selected into representative samples set as iteration grew. On the other side, the blue noise samples were ruled out gradually by the AM algorithm, which meant AM could effectively select the representative samples.

# Case 1.2 Implicit normal mixture model (multiple variables)

![plot3](https://user-images.githubusercontent.com/64850893/86276657-3fa40a80-bba3-11ea-85c2-fddd09cad604.jpg)

From the above table, it's noticeable that more complex covariance structure in this experiment did not lead to obvious reduction of accuracy for AM and FAM. By comparing AM with FAM, AM was computationally more efficient than FAM, which was due to a close form of true parameter in this model.


# Case 2.1 Implicit linear mixture model (1 variable)

![plot5](https://user-images.githubusercontent.com/64850893/86277218-54cd6900-bba4-11ea-867d-503cb6a8a3e6.jpg)

For this case, the AM algorithm successfully selected almost all the red representative samples and few blue noise samples after 7 iterations.

# Case 2.2 Implicit linear mixture model (multiple variables)

![plot6](https://user-images.githubusercontent.com/64850893/86277385-9827d780-bba4-11ea-86e2-7e4b22c3286f.jpg)

In this experiment, the performance of all 5 methods was similar. AM  outperformed other methods in terms of DEV and selection indexes. Due to the similar reasons, MLE, EM and K-means failed to obtain accurate estimate and representative samples set. As for FAM, although it remained good selection performance, its estimation accuracy was lower than that in previous experiments.


# Case 3.1 Implicit logistic mixture model (1 variable)

![plot7](https://user-images.githubusercontent.com/64850893/86277737-39af2900-bba5-11ea-80f6-e650a45a2c47.jpg)

In this scenario, the AM algorithm successfully selected a majority of the representative samples with only few noise samples after 4 iterations.


# Case 3.2 Implicit logistic mixture model (multiples variables)

![plot8](https://user-images.githubusercontent.com/64850893/86277895-7ed35b00-bba5-11ea-9c4f-0f28972eb7d1.jpg)

EM, AM and FAM performed well in selecting representative samples in terms of relatively high PSR and low FDR. From the perspective of DEV, EM performed more accurately than AM and FAM. The result of EM was different from that in previous models, which was because the implicit model assumption of EM matched the true model setup. 


# Case 4.1 Implicit Poisson mixture model (1 variable)

![plot9](https://user-images.githubusercontent.com/64850893/86278141-e1c4f200-bba5-11ea-895e-17d92ee26e12.jpg)

In this case, most of the representative samples with few noise samples were selected by AM after 5 iterations.

# Case 4.2 Implicit Poisson mixture model (multiple variables)
![plot10](https://user-images.githubusercontent.com/64850893/86278315-310b2280-bba6-11ea-9a2d-8ebce871ea62.jpg)

In this experiment, AM and FAM performed better than other methods when estimation and selection were together considered. It's noteworthy that, by the above data generation methods, some noise samples were close to the true model, which meant a small amount of noise samples would be chosen as representative samples.

# Conclusion

According to the experiments in 4 types of mixed datasets, where the explicit models were normal distribution, linear regression, logistic regression and Poisson regression, respectively, as well as 3 different correlation structures among variables, the AM algorithm was robust, outperforming other methods in terms of estimation accuracy and selection consistency in most cases.

# Future work

In our future investigation, I will attempt to study the specific influence of initial parameter input in the AM algorithm. Besides, I will design a general and
data-driven threshold rule. In addition, the in-depth theoretical guarantees of AM is another important point in my future work.



