%Example script for running EChemv1_1. 

%Import 1000 random molecules from ChEMBL
decoy = importdata('datMat_chembl_Morgan3'); 

%parameters 
thres = 0.95; % see EChem.m for what this means 
percent_var_set = 0.2; %take 20% as verification 


A = importdata('CHEMBL210'); 
numel = size(A,1); 
    
p = randperm(numel); 
numel_verification = round(percent_var_set*numel);
verification = A(1:numel_verification,:);  
training = A(numel_verification+1:end,:);
[false_neg, false_pos, eval, evect] = EChem(training,verification,decoy,thres); 
