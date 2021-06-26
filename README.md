# README for MEG Retinotopy code

(Last update by EK: June 25, 2021)

## DESCRIPTION
This project is an MEG retinotopy project, as a collaboration between 
NYU (Eline Kupers, Noah Benson, and Jon Winawer) and the Amsterdam Spinoza Center 
(Akhil Edadan, Maartje de Jong, Wietske Zuiderbaan, Serge Dumoulin).

Historically, this project was initiated by Serge Dumoulin and Jon Winawer, and started by Barrie Klein.
The code from Barrie Klein is archived in this repository under archive/archive_master).


## PREPRINT
Our analyses and results are described in this preprint on BiorXiv. Please cite this paper when using code.

TITLE:		A Population Receptive Field Model of the Magnetoencephalography Response.

AUTHORS:	Eline R Kupers\*, Akhil Edadan\*,  Noah C Benson, Wietske Zuiderbaan, Maartje C de Jong, Serge O Dumoulin\*, Jonathan Winawer\*

(\* indicates shared first or senior authorship)

YEAR:		2020

DOI:		https://doi.org/10.1101/2020.08.28.272534

URL: 		https://www.biorxiv.org/content/10.1101/2020.08.28.272534v1


## GOAL
The aim of this project is to predict visually-triggered MEG responses responses to retinotopic bar stimuli. To do this, we use the population receptive field (pRF) models estimated with functional MRI and combine these with a forward model of current propagation from the cortex to the MEG sensors. 


## MODEL OVERVIEW
1. Preprocess fMRI and MEG data.
2. Derive pRF models on the cortex using fMRI data.
3. Use pRF models to predict cortical responses for MEG stimulus (sweeping bars & blanks).
4. Compute predicted MEG response by multiplying predicted cortical responses with the gain matrix of the forward/head model.
5. Using split-half cross-validation to fit predictions of the predicted MEG responses (concatenated phase-referenced SSVEF responses from 1-s epochs of a given run containg 5 bar sweeps + blanks)
6. Average cross-validated data and model predictions and compute goodness of fit (R-squared) 


## FOLDER STRUCTURE
This code repository contains 6 subfolders specific for this project:
1. analysis --- Folder with code for multiple stages of the analysis:
	- preprocessing: code for MEG, MRI and eye tracking analyses
	- StimRefFwdMdl_analysis: code for stimulus referred MEG forward model, after data are preprocessed.
	- synthetic: code for building a synthetic dataset (still under construction)		

2. data --- Empty folder where we make symbolic links to 3 existing folders:
	- "Retinotopy", a folder with individual subject MEG data
	- "Freesurfer_subjects", a folder with FreeSurfer autorecon segmentation for each subject
	- "brainstorm_db", a folder linking to Brainstorm's Database, containing anat and data folder with subject MEG-MRI alignment and headmodel.

3. figurescripts --- Folder with functions to reproduce main and supplementary manuscript figures, several helper functions and some loose visualization figures.

4. Stimulus --- Code to create and run stimulus for MEG and MRI

## MAIN FUNCTIONS 
* mprf_addPaths.m
	+ Function to add data folders and relevant code (used by ToolboxToolbox)
* mprf_main.m
	+ Function to run analysis from preprocessing MEG & MRI data to forward modeling to fitting
* mprf_makeManuscriptFigures.m
	+ Function to that makes all manuscript figures by calling subfunctions
* mprf_rootPath.m
	+ Function to locate root of this toolbox (must be at this folder level)
* mprf_runAllAnalyses.m
	+ Function to call mprf_main.m with different input variables/settings to run the three analyses in the manuscript: using original pRF params, scaled pRF sizes, or rotated pRF positions


## DATA 
Currently, raw and derivatives of MRI and MEG data are on the Winawerlab server.
Soon, we will deidentify data and permantently store it on a publicly available OSF URL.

Symbolic links point to the following places:
* MEG/MRI Data: server > Projects > MEG > Retinotopy
* Freesurfer:   server > Freesurfer_subjects
* Brainstorm: 	server > Projects > MEG > brainstorm_db


## DEPENDENCIES
This code is Matlab based (written in R2016b), and has the following depencies:
* FreeSurfer (v5.3)
* Brainstorm (https://github.com/brainstorm-tools/brainstorm3)
    last version commit 49a4f9b6 (June 24, 2021)
* FieldTrip (github.com/fieldtrip/fieldtrip)
    last version commit 5cf4cf4ce (March 9, 2020)
* VistaSoft (github.com/vistalab/vistasoft)
    last run version commit d6d0207 (March 17, 2021)
* meg_utils (github.com/WinawerLab/meg_utils)
    last run version commit f0c35f (December 28, 2020)


## HOW TO GET STARTED
(Here I assume MRI/MEG data are preprocessed)

1. Pull code directory from GitHub with the ToolboxToolbox (https://github.com/toolboxhub/toolboxtoolbox) and add to paths
`tbUse('retMeg')`
2. Set up symbolic links in data folder to 'Retinotopy' for MEG/MRI data, 'Freesurfer_subjects' for Freesurfer subject segmentation folders, and 'brainstorm_db' for the Brainstorm database.
3. Define a subject ID:
`subjectID = wlsubj004;`
4. Run run the forward model for subject S1, using the original pRF parameters estimated with fMRI
`results = mprf_main(subjID)`
5. Visualize results:
 `dirPth = loadPaths(subjectID); opts = getOpts; makeFigure4(dirPth, opts);`

## HOW TO REPRODUCE ANALYSES AND FIGURES
1. To run all analyses from the paper:
`mprf_runAllAnalyses.m`
2. Remake all figures from the paper:
`mprf_makeManuscriptFigures.m`


## HOW TO PREPROCESS MRI/MEG DATA
Current data on the OSF webpage (https://osf.io/c3hxj/) are already preprocessed. We tried to share as many data as possible, but because OSF has a storage limit we aren't able to provide raw files (MEG SQD, MRI EPI or T1 niftis, brainstorm database data and anatomy files). We might upload those via a different medium in the future. But in the mean time, if you want to access please reach out to the first or last author of this MEG Retinotopy paper.

* If you want to preprocess MEG data from raw sqd to epoched data, set `opt.doMEGPreproc = true`.
* If you want to preprocess MRI data from mrVista pRF solutions to pRF solutions on FreeSurfer surface, set `opt.doMRIPreproc = true`.
* If you want to go from distortion corrected nifti to pRF solutions in mrVista, see `s_preprocessMRIRetinotopyData.m`
* If you learn how to go from uncorrected EPIs to corrected niftis in BIDS format, see our lab's preprocessing toolbox:  https://github.com/WinawerLab/MRI_tools/blob/master/preprocessing/prisma_preproc.py
* If you want to learn more about MRI-MEG alignment in Brainstorm, then I have to apologize as there is no script. But please see the Brainstorm website for tutorials. https://neuroimage.usc.edu/brainstorm/Tutorials/ChannelFile

## HOW TO RUN FORWARD MODEL
1. Fit model with initially estimated fMRI pRF parameters
`subjID = 'wlsubj004'; mprf_main(subjID);`

2. Fit model for each altered pRF size
`opt = getOpts('perturbOrigPRFs','size'); subjID = 'wlsubj004'; mprf_main(subjID, opt);`

3. Fit model for each altered pRF size
`opt = getOpts('perturbOrigPRFs','position'); subjID = 'wlsubj004'; mprf_main(subjID, opt);`

