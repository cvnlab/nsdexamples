%% Example 1: Basic exploration of the NSD data files

%% Introduction

% In this script, we are going to do some initial exploration of the
% types of data files that are made available as part of the prepared
% NSD data. For this, we are going to use ITK-SNAP (as a volume viewer)
% and freeview (as a surface viewer).
%
% Skills/concepts: 
% - How to use ITK-SNAP and freeview
% - The nature of volumes and surfaces
% - The different types of FreeSurfer surfaces
% - MNI and fsaverage spaces



%% Use ITK-SNAP to explore volume-format data

% Load the following file into ITK-SNAP:
%   ~/nsd/nsddata/ppdata/subj01/anat/T1_0pt8_masked.nii.gz
%%

% Play around in ITK-SNAP and know how to:
% - Navigate / Zoom / Pan
% - Change the colormap and color range
% - View the NIFTI header information (e.g. resolution, origin, orientation)
%%

% Load as a second image the following file:
%   ~/nsd/nsddata/ppdata/subj01/anat/roi/visualsulc.nii.gz
%%

% Play around and know how to:
% - Flip through different volumes ([ or ])
% - Overlay volumes using a thresholding approach
% - Overlay volumes using transparency
% - Toggle overlays on/off (w)
%%

% Draw an ROI and save the ROI out as a NIFTI file called:
%   testroi.nii.gz
% Then load that file back into ITK-SNAP.
%%

% Start over. Load T1_0pt8_masked.nii.gz as the main image,
% and load roi/visualsulc.nii.gz as the segmentation image.
%%

% Using Label Editor, hide all of the labels and show only
% Label 3 (fusiform gyrus). Play around with the 3D volume rendering.
% Information on the meaning of the values in visualsulc is in:
%   ~/nsd/nsddata/freesurfer/fsaverage/label/visualsulc.mgz.ctab
%%

% Quit ITK-SNAP and try launching it from the command-line:
%   cd ~/nsd/nsddata/ppdata/subj01/anat/
%   itksnap T1_0pt8_masked.nii.gz
%%

% There are many other volume viewers out there.
% If you are curious, you can try, for example, mricro.



%% Inspect NSD anatomical data (anat)

% The prepared NSD data come in 3 anatomical spaces (0.5 mm, 0.8 mm, and 1.0 mm).
% These different spaces are exactly coincident and share exactly the same 
% field-of-view and origin. (The origin is placed at the center of each volume.)

% Load the 0.5-mm and 1.0-mm T1 volumes as a main image and an additional image,
% respectively. Compare the two volumes.
%%

% Load the 0.5-mm and 1.0-mm T1 volumes in two separate ITK-SNAP sessions with
% synchronized views (see Preferences). Play around a little.
%%

% Assess the co-registration quality between the 0.8-mm T1 and 0.8-mm T2.
% Notice that compared to the T2, the T1 has a little bit of dropout in 
% the ventral part of the frontal lobe.
%%

% Inspect the Kastner atlas (roi/Kastner2015.nii.gz) overlaid on the 0.8-mm T1.
% Information on the meaning of the values is in:
%   ~/nsd/nsddata/freesurfer/fsaverage/label/Kastner2015.mgz.ctab
%%

% Explore and visualize some additional volumes:
% - aseg (anatomical segmentation from FreeSurfer)
% - brainmask (liberal brain mask used to de-identify data and reduce file sizes)
% - EPI_to_anat1pt0 (mean EPI volume that has been warped to the anat1pt0 space)
%%

% Thus far, we have been visualizing data from a single subject.
% To facilitate comparison across subjects, we can put subjects in a common space.
% A typical one is MNI space, and the prepared NSD data already include some
% MNI transformations for your convenience.

% To determine the transformation between an individual subject and MNI space,
% the T1 from that subject has been nonlinearly warped to an MNI template.
% Furthermore, a version of each subject's T1 that reflects this warping (and
% resampling) has been saved. Let's use ITK-SNAP to quickly look at the MNI
% template and each individual subject's warped T1 volume. 
%   cd ~/nsd/nsddata/
%   itksnap -g templates/MNI152_T1_1mm.nii.gz -o ppdata/subj*/anat/T1_to_MNI.nii.gz
% How good is the alignment?
%%



%% Inspect NSD functional data (func1pt8mm)

% The prepared NSD data come in 2 functional spaces (1.8 mm, 1.0 mm).
% These spaces are aligned to each other, but have slightly different field-of-views 
% and origins (see 'nsddata description' for details). Note that the functional spaces
% are not the same as the anatomical spaces.

% Load floc_facestval.nii.gz (t-value for the contrast of faces vs. non-faces from
% the floc category localizer experiment) and overlay it on mean.nii.gz.
%%

% Alternatively, overlay the t-values on T1_to_func1pt8mm.nii.gz (this is a
% version of the T1 that has already been warped to the func1pt8mm space).
%%

% Explore and visualize some additional volumes:
% - R2 (a measure of signal quality from the NSD experiment)
% - prf_eccentricity (estimate of pRF eccentricity from the prf experiment)
% - valid (fraction of NSD scan sessions in which valid data were recorded)
% - R2_session??.nii.gz (R2 from individual NSD scan sessions)
%   Hint: itksnap -g mean.nii.gz -o R2_*.nii.gz
%%



%% Use freeview to explore surface-format data

% Compared to volume data, surface data are a bit trickier to
% visualize and manage. We will use FreeSurfer's freeview as
% a simple surface viewer, though there are many alternatives.

% Load subj01's left hemisphere inflated surface:
%   ~/nsd/nsddata/freesurfer/subj01/surf/lh.inflated
%%

% Play around in freeview and know how to:
% - Navigate / Rotate / Zoom / Pan
%%

% Load as an overlay (generic):
%   ~/nsd/nsddata/freesurfer/subj01/label/lh.flocfacestval.mgz
%%

% Play around and know how to:
% - Change the colormap and color range
%%

% Set the Render to 'Surface & mesh' to see the faces and vertices
% that comprise a surface.
%%

% Set the Render back to Surface and load in:
%   ~/nsd/nsddata/freesurfer/subj01/label/lh.visualsulc.mgz
% Develop a good visualization of this atlas.
% Save a screenshot.
%%

% Now visualize the same data (visualsulc) on some 
% other versions of the cortical surface:
%   lh.white, lh.pial, lh.sphere
% Use the same color settings as in the lh.inflated case and
% save some screenshots to facilitate comparison.
%%

% Finally, to illustrate the folding-based registration
% that underlies FreeSurfer's fsaverage surface, take
% screenshots of the following (be careful to keep
% the camera view fixed):
%   (1) subj01's lh.sphere     - curvature
%   (2) subj01's lh.sphere     - Kastner2015
%   (3) subj01's lh.sphere.reg - curvature
%   (4) subj01's lh.sphere.reg - Kastner2015
%   (5) fsaverage's lh.sphere  - curvature
%   (6) fsaverage's lh.sphere  - Kastner2015
% Compare these screenshots and notice how sphere.reg is a 
% surface such that visualizing subj01's curvature on that surface
% yields a map that looks well matched to fsaverage's curvature. 
%%

% Like MNI, fsaverage is a common space that can be used to 
% compare results across subjects. MNI is a volume-based space
% that is semi-accurate for cortex and is most accurate for 
% subcortical structures. In contrast, fsaverage is a surface-based
% space that is applicable only to cortex.


