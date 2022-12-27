% Calculates the scalp to cortex distance based on FSL BET generated
% surfaces 
% STEPS
% 1) Load scalp nifti and find non-zero surface mask voxels
% 2) Load brain envelop nifti and find non-zero surface mask voxels
% 3) Find voxel from brain envelop closest to target 
% 4) Find shortest distance from voxel in 3) to scalp voxel

% adpated by Xin to work on the standard brain.
% refer to this post for MNI to voxel coordinate transform:
% https://www.alivelearn.net/?p=1434
% Dec-23-2022.


%% Set path and load common variables
clear;
addpath('NifTI_20140122/');

stimulate_target_name = ["L_Fp", "R_Fp", "L_aMFG", "R_aMFG", "L_pMFG", "R_pMFG", ...
                         "R_IFJ", "R_FEF", "R_M1", "R_preSMA", "R_IPL"];
                     
stimulate_target_mni = [-24, 60, -2;
                        24, 60, -2
                        -32, 42, 34; 
                        30, 50, 26;
                        -38, 22, 38;
                        46, 26, 38;
                        46, 8, 28;
                        34, 6, 62;
                        40, -18, 64;
                        6, 2, 68;
                        48, -54, 46]; 
                    
% V = spm_vol('MNI152_T1_2mm.nii.gz');  
% V.mat is the transformation matrix, but there seems to be a small error. 
% We correct it by define the transform mat manually.

T = [-2,    0,    0,   90;
      0,    2,    0, -126;
      0,    0,    2,  -72;
      0,    0,    0,    1];

stimulate_target_voxel = mni2cor(stimulate_target_mni, T);
                    
%% find FSL BET processed nifti files
brain_env_filename = 'standard_bet_overlay.nii.gz';
scalp = [];

brain_env = load_untouch_nii(brain_env_filename);
% binarize the img to create envelope sphere:
max_val = max(max(brain_env.img(:,:,50)));
brain_env.img(brain_env.img < max_val) = 0;
brain_env.img(brain_env.img >= max_val) = 1;

brain_env_x = [];
brain_env_y = [];
brain_env_z = [];
for kz = 1:size(brain_env.img,3)
    [X, Y] = find(brain_env.img(:,:,kz)>0);
    Z = ones(length(X),1)*kz;
    brain_env_x = [brain_env_x; X];
    brain_env_y = [brain_env_y; Y];
    brain_env_z = [brain_env_z; Z];
end

brain_env_voxel = stimulate_target_voxel;
% voxel size in each dimension (get with fslinfo):
dim_x = 2;
dim_y = 2;
dim_z = 2;

for i = 1: length(stimulate_target_name)
    min_dist = 999; %an impossible large distance
    for k = 1:length(brain_env_x)
        currdist = sqrt(((brain_env_x(k) - stimulate_target_voxel(i, 1)) * dim_x).^2 + ...
            ((brain_env_y(k) - stimulate_target_voxel(i, 2)) * dim_y).^2 + ...
            ((brain_env_z(k) - stimulate_target_voxel(i, 3)) * dim_z).^2);
        if currdist < min_dist
            min_dist = currdist;
            brain_env_voxel(i,:) = [brain_env_x(k), brain_env_y(k), brain_env_z(k)];
        end
    end

    if min_dist == 999
        error('no voxel found');
    end
end

brain_env_mni = cor2mni(brain_env_voxel, T);

writetable(table(stimulate_target_name', brain_env_mni), 'brain_env_mni.csv')
writetable(table(stimulate_target_name', brain_env_voxel), 'brain_env_voxel.csv')


