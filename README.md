# README for MEG Retinotopy code

(Last update by EK: Oct 7, 2020)

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
This code repository contains 5 subfolders specific for this project:
1. analysis --- Folder with code for multiple stages of the analysis:
	- preprocessing: code for MEG, MRI and eye tracking analyses
	- StimRefFwdMdl_analysis: code for stimulus referred MEG forward model, after data are preprocessed.
	- synthetic: code for building a synthetic dataset (still under construction)		

2. data --- Empty folder where we make symbolic links to 3 existing folders:
	- "Retinotopy", a folder with individual subject MEG data
	- "Freesurfer_subjects", a folder with FreeSurfer autorecon segmentation for each subject
	- "brainstorm_db", a folder linking to Brainstorm's Database, containing anat and data folder with subject MEG-MRI alignment and headmodel.

3. figurescripts --- Folder with functions to reproduce main and supplementary manuscript figures, several helper functions and some loose visualization figures.

4. obsolete --- Folder with old functions from Barrie's legacy, these are not used anymore and should eventually be removed.

5. Stimulus --- Code to create and run stimulus for MEG and MRI

6. Synthetic --- Folder to save synthetic data produced by the separated analysis function

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
We are in the process of deidentifying data and permantently store it on a publicly available OSF URL.

Symbolic links point to the following places:
* MEG/MRI Data: server > Projects > MEG > Retinotopy
* Freesurfer:   server > Freesurfer_subjects
* Brainstorm: 	server > Projects > MEG > brainstorm_db


## PIPELINE
This code is Matlab based, and has the following depencies:
* FreeSurfer (v5.3)
* Vistasoft/mrVista (v??)
* ToolboxToolbox (v??)
* Fieldtrip (v??)
* Brainstorm (v??)
* meg_utils (v??)


## HOW TO GET STARTED
(Here I assume MRI/MEG data are preprocessed)

1. Pull code directory from GitHub with the ToolboxToolbox and add to paths
`tbUse('retMeg')`
2. Set up symbolic links in data folder to MEG/MRI data, Freesurfer subject segmentation folders, Brainstorm database.
3. Define a subject ID:
`subjectID = wlsubj004;`
4. Run run the forward model for subject S1, using the original pRF parameters estimated with fMRI
`results = mprf_main(subjID)`
5. Visualize results:
 `dirPth = loadPaths(subjectID); opts = getOpts; makeFigure4(dirPth, opts);`



## HOW TO PREPROCESS MRI/MEG DATA
** coming soon / under construction **
