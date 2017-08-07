% The EChem algorithm v1.1
% Alpha Lee 15/04/2016
%
% Changes: outputs e-vals and e-vects 
%
%%%%%%%%%%%%%%%%%%%%%%%%
% Input 
% 
% training: training set data, n x p matrix where n is the number of
% compounds and p is the number of bits in the fingerprint
% 
% verification: verification set data, m x p martix where m is the number of
% compounds in the verification set 
%
% decoy: r x p matrix, the set of r molecules that are known to NOT bind 
% to the target. 
% 
% thres: percentile threshold that we use to determine whether a ligand is
% in the binding set. 
%
% e.g. thres = 0.75 means that the ligand is predicted to bind if the 
% distance between its fingerprint and the low-dimensional reconstruction 
% is less than the 75% percentile of the distance between the fingerprints
% in the training set and their low-dimensional reconstructions. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%
% Output 
% 
% false_neg: fraction of verification set that is predicted NOT to bind 
% (false negative) 
% 
% false_pos: fraction of decoy set that is predicted to bind (false positive)
%
%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

function [false_neg, false_pos, eval, evect] = EChem(training,verification,decoy,thres) 
  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%processing the training set 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%get rid of columns with the same entry 
tset_cleaned = training(:,any(diff(training,1))); 

%compute z score 
[tset_cleaned_z, mu, sigma] = zscore(tset_cleaned);  


%get covarience matrix and eigenvalues 
covar = cov(tset_cleaned_z);

[v, d] =eig(covar);
l = diag(d);

%Use the MP bound to get the number of significant eigenvalues  
p = size(tset_cleaned,2); 
n = size(tset_cleaned,1); 

num_eig = length(l(find(l>(1+sqrt(p/n))^2))); 


%get the significnt eigenvectors 
vv = v(:,end-num_eig+1:end); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%processing the verification set 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%mean center and scale the verification set; remove columns that have zero
%variance 

veriset_mu = (verification(:,any(diff(training,1)))-repmat(mu,size(verification,1),1))./repmat(sigma,size(verification,1),1);   

%compute coefficients  
coeff = veriset_mu*vv; 

%project back into the ligand space 
proj_vect = (vv*coeff')';

%the background distribution 
coeff_training = tset_cleaned_z*vv; 
proj_vect_training = (vv*coeff_training')';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%processing the decoy set 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

decoy_mu = (decoy(:,any(diff(training,1)))-repmat(mu,size(decoy,1),1))./repmat(sigma,size(decoy,1),1);   
coeff_decoy = decoy_mu*vv; 
proj_vect_decoy = (vv*coeff_decoy')';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%compute summary statistics 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

norm_test = sqrt(sum((proj_vect-veriset_mu).^2,2));
norm_training = sqrt(sum(abs(proj_vect_training-tset_cleaned_z).^2,2)); 
norm_decoy = sqrt(sum(abs(proj_vect_decoy-decoy_mu).^2,2));

threshold = quantile(norm_training,thres);
 

% compute false negative and false positve 
false_neg = length(norm_test(find(norm_test>threshold)))/length(norm_test);  
false_pos = length(norm_decoy(find(norm_decoy<threshold)))/length(norm_decoy); 
eval = l; 
evect = vv; 

end 