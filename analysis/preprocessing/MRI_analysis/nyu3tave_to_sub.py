# Script to project freesurfer average of NYU 3T Retinotopy data 
# onto individual subject's freesurfer mesh.
# Requires neuropythy (https://github.com/noahbenson/neuropythy) 

import neuropythy as ny
import numpy as np
import os

subName = 'wlsubj081_wVitaminE' # 'wlsubj004', 'wlsubj039', 'wlsubj040', 'wlsubj058','wlsubj068', 'wlsubj070', 'wlsubj081_wVitaminE', 'wlsubj106', 'wlsubj109', 'wlsubj111'
dataPth = "/Volumes/server-1/Projects/MEG/Retinotopy/Data/fMRI"
subFolders = "NYU3T_prfresults/sub-wlsubj_avg_param/ses-nyu3t01"
freeSurferPth = "/Volumes/server-1/Freesurfer_subjects/";
sub = ny.freesurfer_subject(os.path.join(freeSurferPth,subName))

# rh_avg_ang = ny.load(os.path.join(dataPth,subFolders,"rh.angle_adj.mgz"))
# lh_avg_ang = ny.load(os.path.join(dataPth,subFolders,"lh.angle_adj.mgz"))
# rh_avg_ecc = ny.load(os.path.join(dataPth,subFolders,"rh.eccen.mgz"))
# lh_avg_ecc = ny.load(os.path.join(dataPth,subFolders,"lh.eccen.mgz"))

# convert to x/y (so we avoid circular interpolation issues)
# (lh_x, lh_y) = (lh_avg_ecc*np.cos(lh_avg_ang),lh_avg_ecc*np.sin(lh_avg_ang))
# (rh_x, rh_y) = (rh_avg_ecc*np.cos(rh_avg_ang),rh_avg_ecc*np.sin(rh_avg_ang))

rh_avg_x = np.average(ny.load(os.path.join(dataPth,subFolders,"rh.x.mgz")),axis=1)
lh_avg_x = np.average(ny.load(os.path.join(dataPth,subFolders,"lh.x.mgz")),axis=1)
rh_avg_y = np.average(ny.load(os.path.join(dataPth,subFolders,"rh.y.mgz")),axis=1)
lh_avg_y = np.average(ny.load(os.path.join(dataPth,subFolders,"lh.y.mgz")),axis=1)

rh_avg_sigma = np.average(ny.load(os.path.join(dataPth,subFolders,"rh.sigma.mgz")),axis=1)
lh_avg_sigma = np.average(ny.load(os.path.join(dataPth,subFolders,"lh.sigma.mgz")),axis=1)

rh_avg_ve = np.average(ny.load(os.path.join(dataPth,subFolders,"rh.vexpl.mgz")),axis=1)
lh_avg_ve = np.average(ny.load(os.path.join(dataPth,subFolders,"lh.vexpl.mgz")),axis=1)

rh_avg_beta = ny.load(os.path.join(dataPth,subFolders,"rh.beta.mgz"))
lh_avg_beta = ny.load(os.path.join(dataPth,subFolders,"lh.beta.mgz"))


# Now, load the fsaverage
fsa = ny.freesurfer_subject(os.path.join(dataPth,"NYU3T_prfresults","fsaverage")) # (can use full path also)

# Now, we just interpolate over:
sub_lh_x = fsa.lh.interpolate(sub.lh, lh_avg_x)
sub_lh_y = fsa.lh.interpolate(sub.lh, lh_avg_y)
sub_rh_x = fsa.rh.interpolate(sub.rh, rh_avg_x)
sub_rh_y = fsa.rh.interpolate(sub.rh, rh_avg_y)
sub_lh_sigma = fsa.lh.interpolate(sub.lh, lh_avg_sigma)
sub_rh_sigma = fsa.rh.interpolate(sub.rh, rh_avg_sigma)
sub_lh_ve = fsa.rh.interpolate(sub.lh, lh_avg_ve)
sub_rh_ve = fsa.rh.interpolate(sub.rh, rh_avg_ve)
sub_lh_beta = fsa.rh.interpolate(sub.lh, lh_avg_beta)
sub_rh_beta = fsa.rh.interpolate(sub.rh, rh_avg_beta)

# Convert back to polar angle and eccen, start from right horz as 0, moving CCW
(lh_ang, lh_ecc) = ny.as_retinotopy({'x':sub_lh_x, 'y':sub_lh_y},'visual')
(rh_ang, rh_ecc) = ny.as_retinotopy({'x':sub_rh_x, 'y':sub_rh_y},'visual')
# (lh_x, lh_y) = ny.as_retinotopy({'x':sub_lh_x, 'y':sub_lh_y},'geographical')
# (rh_x, rh_y) = ny.as_retinotopy({'x':sub_rh_x, 'y':sub_rh_y},'geographical')

# And save out files:
subName = 'wlsubj081'
if not os.path.isdir(os.path.join(dataPth, subName, "NYU3Tave_interp")):
	os.mkdir(os.path.join(dataPth, subName, "NYU3Tave_interp"))
ny.save(os.path.join(dataPth, subName, "NYU3Tave_interp/lh.polar_angle_output-file.mgz"), lh_ang)
ny.save(os.path.join(dataPth, subName, "NYU3Tave_interp/lh.eccentricity_output-file.mgz"), lh_ecc)
ny.save(os.path.join(dataPth, subName, "NYU3Tave_interp/rh.polar_angle_output-file.mgz"), rh_ang)
ny.save(os.path.join(dataPth, subName, "NYU3Tave_interp/rh.eccentricity_output-file.mgz"), rh_ecc)

ny.save(os.path.join(dataPth, subName, "NYU3Tave_interp/lh.x_output-file.mgz"), sub_lh_x)
ny.save(os.path.join(dataPth, subName, "NYU3Tave_interp/rh.x_output-file.mgz"), sub_rh_x)

ny.save(os.path.join(dataPth, subName, "NYU3Tave_interp/lh.y_output-file.mgz"), sub_lh_y)
ny.save(os.path.join(dataPth, subName, "NYU3Tave_interp/rh.y_output-file.mgz"), sub_rh_y)

ny.save(os.path.join(dataPth, subName, "NYU3Tave_interp/lh.sigma_output-file.mgz"), sub_lh_sigma)
ny.save(os.path.join(dataPth, subName, "NYU3Tave_interp/rh.sigma_output-file.mgz"), sub_rh_sigma)

ny.save(os.path.join(dataPth, subName, "NYU3Tave_interp/lh.varexplained_output-file.mgz"), sub_lh_ve)
ny.save(os.path.join(dataPth, subName, "NYU3Tave_interp/rh.varexplained_output-file.mgz"), sub_rh_ve)

ny.save(os.path.join(dataPth, subName, "NYU3Tave_interp/lh.beta_output-file.mgz"), sub_lh_beta)
ny.save(os.path.join(dataPth, subName, "NYU3Tave_interp/rh.beta_output-file.mgz"), sub_rh_beta)
