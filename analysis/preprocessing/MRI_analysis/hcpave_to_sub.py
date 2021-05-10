# Script to project freesurfer average of HCP 7T Retinotopy data 
# onto individual subject's freesurfer mesh.
# Requires neuropythy (https://github.com/noahbenson/neuropythy) 

import neuropythy as ny
import numpy as np
import os

subName = 'wlsubj030' # 'wlsubj004', 'wlsubj039', 'wlsubj040', 'wlsubj058','wlsubj068', 'wlsubj070', 'wlsubj081', 'wlsubj106', 'wlsubj109', 'wlsubj111'
dataPth = "/Volumes/server-1/Projects/MEG/Retinotopy/Data/fMRI"
freeSurferPth = "/Volumes/server-1/Freesurfer_subjects/";
sub = ny.freesurfer_subject(os.path.join(freeSurferPth,subName))

rh_avg_ang = ny.load(os.path.join(dataPth,"HCP_prfresultsmgz/999999/rh.fit1_ang.mgz"))
lh_avg_ang = ny.load(os.path.join(dataPth,"HCP_prfresultsmgz/999999/lh.fit1_ang.mgz"))
lh_avg_ecc = ny.load(os.path.join(dataPth,"HCP_prfresultsmgz/999999/lh.fit1_ecc.mgz"))
rh_avg_ecc = ny.load(os.path.join(dataPth,"HCP_prfresultsmgz/999999/rh.fit1_ecc.mgz"))

rh_avg_sigma = ny.load(os.path.join(dataPth,"HCP_prfresultsmgz/999999/rh.fit1_rfsize.mgz"))
lh_avg_sigma = ny.load(os.path.join(dataPth,"HCP_prfresultsmgz/999999/lh.fit1_rfsize.mgz"))
lh_avg_beta = ny.load(os.path.join(dataPth,"HCP_prfresultsmgz/999999/lh.fit1_gain.mgz"))
rh_avg_beta = ny.load(os.path.join(dataPth,"HCP_prfresultsmgz/999999/rh.fit1_gain.mgz"))

lh_avg_ve = ny.load(os.path.join(dataPth,"HCP_prfresultsmgz/999999/lh.fit1_R2.mgz"))
rh_avg_ve = ny.load(os.path.join(dataPth,"HCP_prfresultsmgz/999999/rh.fit1_R2.mgz"))

# Convert deg to radians
lh_avg_ang = np.pi/180 * (lh_avg_ang)
rh_avg_ang = np.pi/180 * (rh_avg_ang)

# convert to x/y (so we avoid circular interpolation issues)
(lh_x, lh_y) = (lh_avg_ecc*np.cos(lh_avg_ang),lh_avg_ecc*np.sin(lh_avg_ang))
(rh_x, rh_y) = (rh_avg_ecc*np.cos(rh_avg_ang),rh_avg_ecc*np.sin(rh_avg_ang))

# Now, load the fsaverage
fsa = ny.freesurfer_subject(os.path.join(freeSurferPth,'fsaverage')) # (can use full path also)
# and whatever subject you plan to use
# sub = ny.freesurfer_subject(subName) # (can use full path also)

# Now, we just interpolate over:
sub_lh_x = fsa.lh.interpolate(sub.lh, lh_x)
sub_lh_y = fsa.lh.interpolate(sub.lh, lh_y)
sub_rh_x = fsa.rh.interpolate(sub.rh, rh_x)
sub_rh_y = fsa.rh.interpolate(sub.rh, rh_y)
sub_lh_sigma = fsa.lh.interpolate(sub.lh, lh_avg_sigma)
sub_lh_beta = fsa.lh.interpolate(sub.lh, lh_avg_beta)

sub_rh_sigma = fsa.rh.interpolate(sub.rh, rh_avg_sigma)
sub_rh_beta = fsa.rh.interpolate(sub.rh, rh_avg_beta)

sub_lh_ve = fsa.rh.interpolate(sub.lh, lh_avg_ve)
sub_rh_ve = fsa.rh.interpolate(sub.rh, rh_avg_ve)

# Convert back to polar angle and eccen, start from right horz as 0, moving CCW
(lh_ang, lh_ecc) = ny.as_retinotopy({'x':sub_lh_x, 'y':sub_lh_y},'visual')
(rh_ang, rh_ecc) = ny.as_retinotopy({'x':sub_rh_x, 'y':sub_rh_y},'visual')
# (lh_x, lh_y) = ny.as_retinotopy({'x':sub_lh_x, 'y':sub_lh_y},'geographical')
# (rh_x, rh_y) = ny.as_retinotopy({'x':sub_rh_x, 'y':sub_rh_y},'geographical')

# And save out files:
ny.save(os.path.join(dataPth, subName, "hcpave_interp/lh.polar_angle_output-file.mgz"), lh_ang)
ny.save(os.path.join(dataPth, subName, "hcpave_interp/lh.eccentricity_output-file.mgz"), lh_ecc)
ny.save(os.path.join(dataPth, subName, "hcpave_interp/rh.polar_angle_output-file.mgz"), rh_ang)
ny.save(os.path.join(dataPth, subName, "hcpave_interp/rh.eccentricity_output-file.mgz"), rh_ecc)

ny.save(os.path.join(dataPth, subName, "hcpave_interp/lh.x_output-file.mgz"), sub_lh_x)
ny.save(os.path.join(dataPth, subName, "hcpave_interp/rh.x_output-file.mgz"), sub_rh_x)

ny.save(os.path.join(dataPth, subName, "hcpave_interp/lh.y_output-file.mgz"), sub_lh_y)
ny.save(os.path.join(dataPth, subName, "hcpave_interp/rh.y_output-file.mgz"), sub_rh_y)

ny.save(os.path.join(dataPth, subName, "hcpave_interp/lh.sigma_output-file.mgz"), sub_lh_sigma)
ny.save(os.path.join(dataPth, subName, "hcpave_interp/rh.sigma_output-file.mgz"), sub_rh_sigma)

ny.save(os.path.join(dataPth, subName, "hcpave_interp/lh.beta_output-file.mgz"), sub_lh_beta)
ny.save(os.path.join(dataPth, subName, "hcpave_interp/rh.beta_output-file.mgz"), sub_rh_beta)

ny.save(os.path.join(dataPth, subName, "hcpave_interp/lh.varexplained_output-file.mgz"), sub_lh_ve)
ny.save(os.path.join(dataPth, subName, "hcpave_interp/rh.varexplained_output-file.mgz"), sub_rh_ve)
