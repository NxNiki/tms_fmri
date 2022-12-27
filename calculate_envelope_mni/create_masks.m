% create masks based on mni coordinates on the brain envelop:

% note we got problems on mask with radius >= 11 mm.
% related discussion: https://www.jiscmail.ac.uk/cgi-bin/webadmin?A2=FSL;3431dce9.1409
% try to use nltools in python.

% fslmaths /usr/local/fsl/data/standard/MNI152_T1_2mm -mul 0 -add 1 -roi 127 1 168 1 110 1 0 1 mask_F3_point -odt float
% fslmaths mask_F3_point -kernel sphere 20 -fmean mask_F3_sphere -odt float

clear
clc

setenv('PATH', [getenv('PATH') ':/usr/local/fsl/bin']);

opts = detectImportOptions('brain_env_voxel.csv');
brain_env_voxel = readtable('brain_env_voxel.csv',opts);

% sites = {'L_Fp'; 'R_Fp'; 'L_aMFG'; 'R_aMFG'; 'L_pMFG'; 'R_pMFG'; 'R_IFJ'; 'R_FEF'; 'R_M1'; 'R_preSMA'; 'R_IPL';};

output_dir='brain_envelope_mask/';
if ~exist(output_dir, 'dir')
    mkdir(output_dir)
end

for st=1:height(brain_env_voxel)
    
%     x_cor = brain_env_voxel{st, 2};
%     y_cor = brain_env_voxel{st, 3};
%     z_cor = brain_env_voxel{st, 4};

%     tmpoint = [output_dir, brain_env_voxel{st,1}{1,1}, '_tmpoint'];
%     cmd0 = ['fslmaths MNI152_T1_2mm_brain_mask.nii -mul 0 -add 1 -roi ' ...
%         num2str(x_cor) ' 1 ' num2str(y_cor) ' 1 ' num2str(z_cor) ' 1 0 1 ' tmpoint ' -odt float'];
%     system(cmd0);

%     % 6mm:
%     out1 = [output_dir, brain_env_voxel{st,1}{1,1}, '_6mm_tmp'];
%     cmd1 = ['fslmaths ' tmpoint '.nii.gz -kernel sphere 6 -fmean -bin ' out1 ' -odt float'];
%     system(cmd1);
% 
    out11 = [output_dir, brain_env_voxel{st,1}{1,1}, '_6mm'];
%     cmd11 = ['fslmaths ', out1, ' -mul MNI152_T1_2mm_brain_mask.nii ' out11];
%     system(cmd11);
%     delete([out1 '.nii.gz'])
%     
%     % 10mm:
%     out2 = [output_dir, brain_env_voxel{st,1}{1,1}, '_10mm_tmp'];
%     cmd2 = ['fslmaths ' tmpoint '.nii.gz -kernel sphere 10 -fmean -bin ' out2 ' -odt float'];
%     system(cmd2);
% 
    out22 = [output_dir, brain_env_voxel{st,1}{1,1}, '_10mm'];
%     cmd22 = ['fslmaths ' out2 ' -mul MNI152_T1_2mm_brain_mask.nii ' out22];
%     system(cmd22);
%     delete([out2 '.nii.gz'])
%     
%     % 14mm:
%     out3 = [output_dir, brain_env_voxel{st,1}{1,1}, '_14mm_tmp'];
%     cmd3 = ['fslmaths ' tmpoint '.nii.gz -kernel sphere  -fmean -bin ' out3 ' -odt float'];
%     system(cmd3);
% 
    out33 = [output_dir, brain_env_voxel{st,1}{1,1}, '_14mm'];
%     cmd33 = ['fslmaths ' out3 ' -mul MNI152_T1_2mm_brain_mask.nii ' out33];
%     system(cmd33);
%     %delete([out3 '.nii.gz'])
    
    % 10mm - 6mm:
    out21 = [output_dir, brain_env_voxel{st,1}{1,1}, '_10-6mm'];
    cmd21 = ['fslmaths ' out22 ' -sub ' out11 ' ' out21 ' -odt float'];
    system(cmd21)
    
    % 14mm - 10mm:
    out32 = [output_dir, brain_env_voxel{st,1}{1,1}, '_14-10mm'];
    cmd32 = ['fslmaths ' out33 ' -sub ' out22 ' ' out32 ' -odt float'];
    system(cmd32)

    %delete([tmpoint '.nii.gz'])
    
end
