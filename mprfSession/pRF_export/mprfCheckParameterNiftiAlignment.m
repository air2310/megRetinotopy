function mprfCheckParameterNiftiAlignment(cls_nifti, parameter_nifti_file_name)

% Check that all non-zero values of the parameter nifti are in voxels with
% nonzero class values. This will be true if the images are aligned.

par_ni = niftiRead(parameter_nifti_file_name);
idx = par_ni.data ~= 0;

assert(all(cls_nifti.data(idx)))

% Check the headers
assertEqual(par_ni.qto_xyz, cls_nifti.qto_xyz);

% There were instances when nifti alignments were off by a very small
% amount. The above line was changed for dealing with this machine
% precision
tol = 10.^-2;
assertElementsAlmostEqual(par_ni.qto_xyz, cls_nifti.qto_xyz,'absolute',tol)

end