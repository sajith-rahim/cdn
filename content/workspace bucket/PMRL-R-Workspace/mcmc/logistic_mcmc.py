# -*- coding: utf-8 -*-
"""
Created on Wed Sep  9 01:25:07 2020

@author: sajithar99
Credit : John Clements
"""


class mcmc_logistic_reg:
    
    import numpy as np
    
    def __init__self(self):
        self.raw_beta_distr = np.empty(1)
        self.beta_distr = np.empty(1)
        self.beta_hat = np.empty(1)
        self.cred_ints = np.empty(1)
    
    def inv_logit(self, beta, X):
        ###
        # A function to undo a logit transformation. Translates log-odds to 
        # probabilities.
        ###
        return (np.exp(np.matmul(X, beta.reshape((-1, 1)))) / 
                (1 + np.exp(np.matmul(X, beta.reshape((-1, 1))))))
    
    def normal_log_prior(self, beta, prior_means, prior_stds):
        ###
        # A function to calculate the log prior using a normal prior. The
        # log prior is used to avoid underflow.
        ###
        from scipy.stats import norm
        import numpy as np
        
        return np.sum(norm.logpdf(beta, loc=prior_means.reshape((-1, 1)), 
                                  scale=prior_stds.reshape((-1, 1)))) 
    
    def log_likelihood(self, y, X, beta):
        ###
        # Defines a function to calculate the log likelihood of the betas given 
        # the data. Uses the log likelihood instead of the normal likelihood to 
        # avoid underflow.
        ###        
        return np.sum(y * np.log(self.inv_logit(beta.reshape((-1, 1)), X)) + 
                      (1-y)*np.log((1-self.inv_logit(beta.reshape((-1,1)),X))))
    
    def log_posterior(self, y, X, beta, prior_means, prior_sds):
        ###
        # Defines a function to calculate the log posterior of the betas given 
        # the log likelihood and the log prior assuming it is a normal prior.
        ###      
        return (self.normal_log_prior(beta, prior_means, prior_stds) + 
                self.log_likelihood(y, X, beta))
    
    def mcmc_solver(self, 
                               y, 
                               X, 
                               beta_priors, 
                               prior_stds,
                               jumper_stds, 
                               num_iter,
                               add_intercept,
                               random_seed):
        
        from scipy.stats import norm
        from tqdm import tqdm
        
        # set a random seed
        np.random.seed(random_seed)
        
        # add an intercept if desired
        if add_intercept:
            X_mod = np.append(np.ones(shape=(X.shape[0], 1)),X, 1)
        else:
            X_mod = X
            
        # creates a list of beta indexes to be looped through each iteration
        beta_indexes = [k for k in range(len(beta_priors))]
        
        ###
        # Initialize beta_hat with the priors. It will be have a number of rows
        # equal to the number of beta coefficients and a number of columns
        # equal to the number of iterations + 1 for the prior. Each row will 
        # hold values of a single coefficient. Each column is an iteration
        # of the algorithm.
        ###
        beta_hat = np.array(np.repeat(beta_priors, num_iter+1))
        beta_hat = beta_hat.reshape((beta_priors.shape[0], num_iter+1))
        
        # perform num_iter iterations
        for i in tqdm(range(1, num_iter + 1)):
            
            # shuffle the beta indexes so the order of the coefficients taking 
            # the Metropolis step is random
            np.random.shuffle(beta_indexes)
            
            # perform the sampling for each beta hat sequentially
            for j in beta_indexes:
                
                # generate a proposal beta using a normal distribution and the 
                # beta_j
                proposal_beta_j = beta_hat[j, i-1] + norm.rvs(loc=0, 
                                                              scale=jumper_stds[j], 
                                                              size=1)
                
                # get a single vector for all the most recent betas
                beta_now = beta_hat[:, i-1].reshape((-1, 1))
                
                # copy the current beta vector and insert the proposal beta_j 
                # at the jth index
                beta_prop = np.copy(beta_now)
                beta_prop[j, 0] = proposal_beta_j
                
                # calculate the posterior probability of the proposed beta
                log_p_proposal = self.log_posterior(y, X_mod, beta_prop, 
                                                    beta_now, prior_stds)
                # calculate the posterior probability of the current beta
                log_p_previous = self.log_posterior(y, X_mod, beta_now, 
                                                    beta_now, prior_stds)
                
                # calculate the log of the r-ratio
                log_r = log_p_proposal - log_p_previous
                
                # if r is greater than a random number from a Uniform(0, 1) 
                # distribution add the proposed beta_j to the list of beta_js
                if np.log(np.random.random()) < log_r:
                    beta_hat[j, i] = proposal_beta_j
                # otherwise, just add the old value
                else:
                    beta_hat[j, i] = beta_hat[j, i-1]
        
        # sets the attribute raw_beta_distr to matrix of beta_hat
        self.raw_beta_distr = beta_hat
    
    def trim(self, burn_in):
        ###
        # This function that trims the distribution of beta hats based on the 
        # burn-in rate
        ###
        start = 1 + round(burn_in * self.raw_beta_distr.shape[1] - 1)
        
        # sets the attribute beta_distr to the raw_beta_distr minus the burn-in
        self.beta_distr = self.raw_beta_distr[:, start:-1]
        
    def credible_int(self, alpha=0.05):
        ###
        # Returns the 100*(1-alpha)% credible interval for each coefficient in
        # a 2 by number of coefficients numpy array
        ###
        from numpy import transpose, quantile
        self.cred_ints =  transpose(quantile(self.beta_distr,
                                             q=(alpha/2, 1-alpha/2),
                                             axis=1))
    
    def fit(self, method = 'median'):
        ###
        # Uses the distribution of betas_hat without the burn-in to either give
        # the median, mean, or mode as an estimate of the beta vector
        ###
        
        from numpy import median, mean
        from scipy.stats import mode
        
        if method == 'median':
            beta_hat = median(self.beta_distr, axis=1).reshape((-1,1))
        elif method == 'mean':
            beta_hat = mean(self.beta_distr, axis=1).reshape((-1,1))
        else:
            beta_hat = mode(self.beta_distr, axis=1)[0]
        
        # sets the beta_hat attribute to either the median, mean, or mode
        self.beta_hat = beta_hat
        
    def predict(self, X_new, add_intercept=True, prob=True):
        ###
        # Gives predictions, either in log-odds or probabilities (default)
        ###
        from numpy import matmul
        
        # add an intercept column if desired
        if add_intercept:
            X_mod = np.append(np.ones(shape=(X_new.shape[0], 1)),X_new, 1)
        else:
            X_mod = X_new
        
        # outputs predicted probabilities if prob == True
        if prob:
            predictions = self.inv_logit(self.beta_hat, X_mod)
        # outputs predicted log-odds otherwise
        else:
            predictions = matmul(X_mod, self.beta_hat)
            
        # returns predictions
        return predictions
    
    def predict_class(self, X_new, add_intercept=True, boundary=0.5):
        ###
        # predicts the class of the new observations based on a decision 
        # boundary for probability. If predicted probability > boundary, it
        # belongs to class 1
        ###
        from numpy import where
        # predict the probabilities
        preds = self.predict(X_new, add_intercept, prob=True)
        # set predictions to 1 or 0 based on boundary
        pred_classes = where(preds > boundary, 1, 0)
        
        return pred_classes

#---------------------CANDY POWER DATASET-------------    
    
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# load in field goal data
all_data = pd.read_csv('candy_power.csv')

# list of independent variables in the model
vars_of_interest = ['pricepercent']
# name of dependent variable
target = 'chocolate'

my_data = all_data[[target] + vars_of_interest]
my_data = all_data.dropna()

# make a column vector numpy array for the dependent variable
y = np.array(my_data[target]).reshape((-1,1))
# make a column vector numpy array for the independent variable
X = np.array(my_data[vars_of_interest]).reshape((-1,len(vars_of_interest)))

# make a numpy array for the beta means priors
beta_priors = np.repeat(0.0, len(vars_of_interest)+1) 
# make a numpy array for the beta standard deviation priors
prior_stds = np.repeat(1, len(vars_of_interest))
# make a row vector numpy array for the standard deviation of the jumper distribution
jumper_stds = np.repeat(0.1, len(vars_of_interest)+1) 

# set the number of iterations to run
num_iter = 50000
# set what early portion to burn out
burn_in = 0.1
# set an alpha for the credible interval
cred_int_alpha = 0.05
# set to true to add an intercept to the data
add_intercept = True
# set the random seed
random_seed = 1

# set the decision boundary for classes
boundary = 0.5

# create an instance of the mcmc_logistic_reg
mcmc_log_mod = mcmc_logistic_reg()
# create the untrimmed distribution of beta hats
mcmc_log_mod.mcmc_solver(y, 
                         X,
                         beta_priors, 
                         prior_stds,
                         jumper_stds, 
                         num_iter,
                         add_intercept,
                         random_seed)
# create a distribution of beta_hats without the burn-in
mcmc_log_mod.trim(burn_in)
# use the median of the beta hats distributions as the coefficients for prediction
mcmc_log_mod.fit(method='median')
# create credible intervals for the 
mcmc_log_mod.credible_int(cred_int_alpha)


#---------------------PLOT-------------


plt.figure(figsize=(12, 5), dpi= 80, facecolor='w', edgecolor='k')

plt.subplot(1, 2, 1)
plt.plot(mcmc_log_mod.raw_beta_distr[0], mcmc_log_mod.raw_beta_distr[1])
plt.title('Simulated Raw Joint Distribution of the Coefficients', fontsize=12)
plt.xlabel('Intercept', fontsize=10)
plt.ylabel('Coefficient of Price Percentile', fontsize=10)

plt.subplot(1, 2, 2)
plt.plot(mcmc_log_mod.beta_distr[0], mcmc_log_mod.beta_distr[1])
plt.title('Simulated Joint Distribution of the Coefficients without Burn-in', fontsize=12)
plt.xlabel('Intercept', fontsize=10)
plt.ylabel('Coefficient of Price Percentile', fontsize=10)
plt.show();


plt.figure(figsize=(12, 5), dpi= 80, facecolor='w', edgecolor='k')

plt.subplot(1, 2, 1)
plt.hist(mcmc_log_mod.beta_distr[0], density=True)
plt.title(f'Simulated Distr. of Intercept with {100 * (1-cred_int_alpha)}% Cred. Int.', 
          fontsize=12)
plt.xlabel('Intercept', fontsize=10)
plt.ylabel('Density', fontsize=10)
plt.axvline(x=mcmc_log_mod.cred_ints[0,0], color='r', linestyle='dashed', linewidth=2)
plt.axvline(x=mcmc_log_mod.cred_ints[0,1], color='r', linestyle='dashed', linewidth=2)

plt.subplot(1, 2, 2)
plt.hist(mcmc_log_mod.beta_distr[1], density=True)
plt.title(f'Simulated Distr. of Price Percentile Coef. with {int(100 * (1-cred_int_alpha))}% Cred. Int.', 
          fontsize=12)
plt.xlabel('Price Percentile Coefficient', fontsize=10)
plt.ylabel('Density', fontsize=10)
plt.axvline(x=mcmc_log_mod.cred_ints[1,0], color='r', linestyle='dashed', linewidth=2)
plt.axvline(x=mcmc_log_mod.cred_ints[1,1], color='r', linestyle='dashed', linewidth=2)
plt.show();

# to make the regression line
x_ = np.arange(my_data[vars_of_interest[0]].min(), 
              my_data[vars_of_interest[0]].max(), 
              step=0.01)
x_ = x_.reshape((-1, 1))
y_ = mcmc_log_mod.predict(x_)

plt.figure(figsize=(12, 5), dpi= 80, facecolor='w', edgecolor='k')

plt.plot(x_, y_, '-', color='r')
plt.scatter(X, y)
plt.hlines(xmin=x_[0], xmax=x_[-1], y=boundary, color='g', linestyle='dashed', linewidth=2)
plt.title('Price Percentile vs. Probability of Being Chocolate', fontsize=15)
plt.xlabel('Price Percentile', fontsize=15)
plt.ylabel('Probability of Being Chocolate', fontsize=15)
plt.show;