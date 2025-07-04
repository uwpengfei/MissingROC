The following files include the code for data generation, the implementation of our proposed method, and four competitive methods

(1) gendata.R: Contains R functions to generate data for Scenarios 1-6.

(2) supp.R: Contains R functions for calculating ROC curves and AUC.

(3) Our.R: Contains R functions to implement the proposed method. 

(4) ipw.cpp and IPW.R: Contains functions for implementing the IPW method. Note that the code is adapted from the implementation available at https://github.com/wbaopaul/AUC-IV.

(5) Others.R: Contains R functions to implement the IG, VER, and Full methods.


(6) Example.R: Contains sample code for using the proposed method, as well as the IPW, IG, VER, and Full methods. 

(7) seed.txt: Contains 1000 random seeds, with each seed used to generate data for one replication.


For more detials on each above method, please refer to Hu, D*, Yu, T and Li, P (2025) Receiver operating characteristic curve analysis with non-ignorable missing disease status. Canadian Journal of Statistics,Accepted.    
