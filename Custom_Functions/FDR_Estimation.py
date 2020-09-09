#Custom function for FDR estimation

import pandas as pd
import numpy as np

#FDR Estimation function for decoys

#Produces an array of scores from reference and query spectra, based on a given score function

def scores_array_decoy(query_spectra, reference_spectra, score_function, mass_tolerance = 1.0):
    
    scores = []
    for query in query_spectra:
        
        for reference in reference_spectra:
            
            score = score_function(query, reference)[0] if parent_mass_match_decoy(query, reference, mass_tolerance) else 0
            scores.append(score)
    return np.array(scores)

#Checks whether two spectra are within a certain parent mass tolerance

def parent_mass_match_decoy(spectrum_1, spectrum_2, tolerance):
    
    parent_mass_1 = spectrum_1.get("parent_mass") 
    parent_mass_2 = spectrum_2.get("parent_mass")
    if not parent_mass_1 or not parent_mass_2:
        return False
    
    elif abs(parent_mass_1 - parent_mass_2) <= tolerance:
        return True
    
    return False

#Getting hits for each query spectrum in reference spectra
#Type_query here is whether it's a target or decoy hit. True is for target, False for decoy hit.

def get_hits_decoy(query_spectra, reference_spectra, score_function, type_hit):
    
    no_query = len(query_spectra)
    hits = pd.DataFrame({'Cosine Score':pd.Series([0]*no_query, dtype='float'), 'Reference InChl':pd.Series([0]*no_query, dtype='str'), 'Reference Id':pd.Series([0]*no_query, dtype='str'), 'Query InChl':pd.Series([0]*no_query, dtype='str'), 'Query Id':pd.Series([0]*no_query, dtype='str'), 'Target/Decoy Hit':pd.Series([0]*no_query, dtype='bool'), 'True/False Hit':pd.Series([0]*no_query, dtype='bool')})
    for i in range(no_query):
        
        query_scores = scores_array_decoy([query_spectra[i]], reference_spectra, score_function)
        max_reference = np.argpartition(query_scores, -1)[-1]
        reference_inchikey = reference_spectra[max_reference].get("inchikey")
        reference_inchikey = reference_inchikey if reference_inchikey else np.nan
        reference_id = reference_spectra[max_reference].get("spectrumid")
        reference_id = reference_id if reference_id else np.nan
        query_inchikey = query_spectra[i].get("inchikey") 
        query_inchikey = query_inchikey if query_inchikey else np.nan
        query_id = query_spectra[i].get("spectrumid")
        query_id = query_id if query_id else np.nan
        hit_identity = (str(query_inchikey)[:14] == str(reference_inchikey)[:14]) if type_hit else np.nan
        
        hit = [query_scores[max_reference], reference_inchikey, reference_id, query_inchikey, query_id, type_hit, hit_identity]
        hits.iloc[i, :] = hit
 
    return hits



#FDR estimation functions for knockoffs

#Produces an array of scores from reference and query spectra, based on a given score function

def scores_array_knockoff(query_embeddings, query_docs, reference_embeddings, 
                 reference_docs, score_function, mass_tolerance = 1.0):
    
    scores = []
    for i, query in enumerate(query_embeddings):
        
        for j, reference in enumerate(reference_embeddings):
            
            score = score_function(query, reference) if parent_mass_match_knockoff(query_docs[i], reference_docs[j], mass_tolerance) else 0
            scores.append(score)
    return np.array(scores)

#Checks wherther two spectra are within a certain parent mass tolerance

def parent_mass_match_knockoff(spectrum_doc_1, spectrum_doc_2, tolerance):
    
    parent_mass_1 = spectrum_doc_1._obj.get("parent_mass") 
    parent_mass_2 = spectrum_doc_2._obj.get("parent_mass")
    if not parent_mass_1 or not parent_mass_2:
        return False
    
    elif abs(parent_mass_1 - parent_mass_2) <= tolerance:
        return True
    
    return False

#Getting hits for each query spectrum in reference spectra
#Type_query here is whether it's a target or decoy hit. True is for target hit, False for decoy hit.

def get_hits_knockoff(query_embeddings, query_docs, reference_embeddings, 
             reference_docs, score_function, type_hit):
    
    no_query = len(query_docs)
    hits = pd.DataFrame({'Cosine Score':pd.Series([0]*no_query, dtype='float'), 'Reference InChl':pd.Series([0]*no_query, dtype='str'), 'Reference Id':pd.Series([0]*no_query, dtype='str'), 'Query InChl':pd.Series([0]*no_query, dtype='str'), 'Query Id':pd.Series([0]*no_query, dtype='str'), 'Target/Decoy Hit':pd.Series([0]*no_query, dtype='bool'), 'True/False Hit':pd.Series([0]*no_query, dtype='bool')})
    for i in range(no_query):
        
        query_scores = scores_array_knockoff([query_embeddings[i]], [query_docs[i]], 
                                    reference_embeddings, reference_docs, score_function)
        max_reference = np.argpartition(query_scores, -1)[-1]
        query_inchikey = query_docs[i]._obj.get("inchikey") 
        query_inchikey = query_inchikey if query_inchikey else np.nan
        query_id = query_docs[i]._obj.get("spectrumid")
        query_id = query_id if query_id else np.nan
        reference_inchikey = reference_docs[max_reference]._obj.get("inchikey")
        reference_inchikey = reference_inchikey if (reference_inchikey and type_hit) else np.nan
        reference_id = reference_docs[max_reference]._obj.get("spectrumid")
        reference_id = reference_id if (reference_id and type_hit) else np.nan
        hit_identity = (str(query_inchikey)[:14] == str(reference_inchikey)[:14]) if type_hit else np.nan
        
        hit = [query_scores[max_reference], reference_inchikey, reference_id, query_inchikey, query_id, type_hit, hit_identity]
        hits.iloc[i, :] = hit
 
    return hits

#Note that different functions are used for getting the decoy and knockoff hits.
#This is because the knockoffs are stored as vectors. 
#Future changes may resolve this by storing knockoffs and decoys as spectrum documents.



#Finding true and estimated FDRs

#Estimate FDR from dataframe of sorted hits
#FDR is estimated as no. decoy hits/no. target hits
#Returns np.nan if no. targets is 0

def estimate_FDR(hits):
    
    targets = sum(hits["Target/Decoy Hit"] == 1)
    decoys = sum(hits["Target/Decoy Hit"] == 0)
    if targets != 0:
        
        FDR_estimate = decoys/targets
        return FDR_estimate
    else:
        
        print("No targets found. No FDR estimate available.")
        return np.nan
    
#Calculate FDR from dataframe of sorted hits
#False hits are hits which are not the true identity for associated query spectra
#FDR is estimated as no. false hits/no. target hits
#Returns np.nan if no. correct targets is 0

def calculate_FDR(hits):
    
    no_false_hits = sum((hits["Target/Decoy Hit"] == 1) & (hits["True/False Hit"] == 0))
    no_correct_targets = sum((hits["Target/Decoy Hit"] == 1) & (hits["True/False Hit"] == 1))
    if no_correct_targets != 0:
        
        FDR = no_false_hits/no_correct_targets
        return FDR
    else:
        
        print("No targets found. No FDR estimate available.")
        return np.nan
