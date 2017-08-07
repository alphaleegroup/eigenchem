# Predicting Protein-Ligand Affinity with a Random Matrix Framework
This repository provides the Matlab script that accompanies our paper entitled "Predicting Protein-Ligand Affinity with a Random Matrix Framework" (A. A. Lee, M. P. Brenner, L. J. Colwell, Proceedings of the National Academy of Sciences, 113, 13564 (2016)). 

## Abstract
Developing computational methods to screen ligands against protein targets is a major challenge for drug discovery. We present a robust mathematical framework, inspired by random matrix theory, which predicts ligand binding to a target given the known ligand set of that target. Our method considers binding prediction as a denoising problem, recognizing that only some of the chemically important features associated with each ligand contribute to binding to a particular receptor. We use correlations among chemical features in the known ligand set, combined with random matrix theory, to eliminate statistically insignificant correlations. Our method outperforms existing algorithms in the literature. We show that our algorithm has the physical interpretation of estimating the ligand–target binding energy.

## Installation
Matlab is needed to run this script. This script has been tested on Matlab R2015b. 

## Usage  
EChem.m is the main function. testGPCR.m is provided as an example of how to use EChem. 
