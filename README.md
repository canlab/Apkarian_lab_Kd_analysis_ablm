# ablm
Wager lab fork of the Apkarian group's ABLM toolbox.

INSTALLING
-------------------------------------------
To build, run `make ablm`.

If you get this error, below:
xcrun: error: invalid active developer path (/Library/Developer/CommandLineTools), missing xcrun at: /Library/Developer/CommandLineTools/usr/bin/xcrun

...then you need to install XCode (on OS X mac)
In terminal:
xcode-select --install

To make accessible from anywhere:
1. Add the ablm directory to your path, 
2. Specify full path when executing, 
3. Copy ablm binary (file called ablm) to your /usr/local/bin folder. it is stand-alone.


TO RUN IT ON IMAGE FILES
-------------------------------------------

It runs only on a single 4-D NIFTI .nii image file.

Get help: 
./ablm

Help contents:

Usage: ablm  -i <input_image_name>   -m <mask_image_name>  -x <roi_image_name> -o <output_image_name>   -t <correlation threshold>  -s <correlation sign>  -r <radius>   -g <link region> -a -d
     -i <input_image_name>: input_image path and name; the name is in format of .nii instead of .nii.gz
     -m <mask_image_name>: calculate correlations of every voxel inside the mask with every other voxel only inside this mask; default mask image is the average input image over time; if want to use defaut mask image, ignore -m.
     -x <roi_image_name>: calculate correlations of every voxel inside the roi with every other voxel only outside this roi but inside the mask.
     -o <output_image_name>: output_image path and name; the name is in format of .nii instead of .nii.gz
     -t <correlation threshold>: the voxel is counted only if the correlation is beyond correlation threshold; default value is 0.25 and its range is between 0 and 1.
     -s <correlation sign>: 1: only count the voxel with greater than positive correlation threshold; -1: only count the voxel with less than negative correlation threshold; 0: count the voxel with above both conditions; default value is 0.
     -r <sphere radius>: sphere radius and its units is voxel; default value is 3.
     -g <link region>: 1: only search the voxels outside of the sphere in mask image; -1: only search the voxels inside of the sphere in mask image ; 0: search all masked voxels; default value is 0.
     -a : generate an average correlation map with correlation threshold > 0.
     -d : generate an average parameter map with correlation threshold > thres .
     Ex. 1. default mask image, R = 0.25, absolute R, whole brain: ablm -i filt_rsn_01.nii -o link_rsn_01.nii
     Ex. 2.  Hippocampus mask image, R = 0.3, positive R, whole brain: ablm -i filt_rsn_01.nii -m Hippo_mask.nii -o link_rsn_01.nii -t 0.3 -s 1
     Ex. 3.  Hippocampus mask image, R = 0.3, positive R, inside of sphere r=3 voxels link: ablm -i filt_rsn_01.nii -m Hippo_mask.nii -o link_rsn.nii -t 0.3 -s 1 -r 3 -g -1
     Ex. 4.  Gray_matter mask image, Hippocampus roi image, R = 0.3, positive R: ablm -i filt_rsn_01.nii -m GM_mask.nii -x Hippo_roi.nii -o link_rsn.nii -t 0.3 -s 1
     Ex. 5.  default mask image to generate an average correlation map: ablm -i filt_rsn_01.nii -o ave_map.nii -a 
     Ex. 6.  default mask image to generate an average distance parameter map with correlation > 0.35: ablm -i filt_rsn_01.nii -o D_map.nii -d -t 0.35
