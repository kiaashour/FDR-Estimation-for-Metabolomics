#Creating Knockoffs from Gaussian Mixture Model

import numpy as np

#Functions required to create knockoffs

#Checking positive definiteness of matrix
#epsilon argument allows for foloating point errors

def is_pos_semi_def(A, epsilon = 1e-10):  
    
    eigs = np.linalg.eig(A)[0]
    min_eig = min(eigs)
    if min_eig >= -epsilon:
        
        return True
    
    return False


#Finding Dk matrix in sampling process

def find_Dk(covariance_matrix, embedding_dimension):
    
    eigs = np.linalg.eig(covariance_matrix)[0]
    min_eig = min(eigs)
    s = min(2*min_eig, 1)
    Dk = np.diag([s]*embedding_dimension)
    return Dk
                  
#Get posterior mixture assignment sample

def posterior_mixture(components, probs):
    
    sample = np.random.choice(components, p=probs)
    return sample



#Create Knockoffs

#Functin to create knockoffs given list of arrays (vectors) and Gaussian mixture model (GMM)

def create_knockoffs(GMM, vectors):
    
    no_comp = len(GMM.weights_)
    embedding_dimension = len(vectors[0])
    covariances = GMM.covariances_
    means = GMM.means_

#Obtaining Dks for knockoff distributiuons

    Dks = []
    for i in range(no_comp):

        cov = covariances[i]
        Dk = find_Dk(cov, embedding_dimension)
        Dks.append(Dk)
        if not is_pos_semi_def(2*cov-Dk):

            print("Joint covariance matrix for knockoff and actual model is not positive semi-definite.")
            print("\n Knockoff generation terminated.")
            return
            
#Getting components of conditional distribution of knockoff

    knock_means_comps_1 = []
    knock_means_comps_2 = []
    knock_covs = []
    Id = np.diag([1]*embedding_dimension)
    for i in range(no_comp):

        Di = Dks[i]
        cov_inv_i = np.linalg.inv(covariances[i]) 
        mean_i = means[i]
        knock_cov = 2*Di - Di@(cov_inv_i@Di)
        knock_mean_comp_1 = Di@(cov_inv_i@mean_i)
        knock_mean_comp_2 = Id - Di@cov_inv_i
        knock_means_comps_1.append(knock_mean_comp_1)
        knock_means_comps_2.append(knock_mean_comp_2)
        knock_covs.append(knock_cov)
        
#Creating knockoffs for each document spectrum

    knockoffs = []
    components = np.arange(no_comp)
    probs = GMM.predict_proba(vectors)
    for i, x in enumerate(vectors):
        
        x_probs = probs[i]
        k_posterior = posterior_mixture(components, x_probs)
        knock_mean = knock_means_comps_1[k_posterior] + knock_means_comps_2[k_posterior]@x
        knock_cov = knock_covs[k_posterior]
        knockoff_sample = np.random.multivariate_normal(knock_mean, knock_cov)
        knockoffs.append(knockoff_sample)
        
    return knockoffs
