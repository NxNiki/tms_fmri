% Calculates the scalp to cortex distance based on FSL BET generated
% surfaces

% NO affine transform necessary when using .MRI file from visor

% STEPS
% 1) Load scalp nifti and find non-zero surface mask voxels
% 2) Load brain envelop nifti and find non-zero surface mask voxels
% 3) Find voxel from brain envelop closest to target
% 4) Find shortest distance from voxel in 3) to scalp voxel

% loop through all subjects
% loop through all files
% loop through all targets

%% Set path and load common variables
username=getenv('USERNAME');

% load files with visor ouptut
mainDir = ['C:\Users\', username, '\Box\Data_Archive\cc_targets2dicom'];
visor_folders = dir([mainDir '\CAUSCON_*_VISOR']);

scriptDir = 'G:\My Drive\POSTDOC Research\Collaborations\CausCon\';

% Load the voxel dimensions associated with each file in visor folder order
load([scriptDir, 'voxdim_162visorfiles.mat']); %voxdim

% Load subject IDs
load([scriptDir, 'scan2_subID.mat']); %subID

% Load targets for each subject in MRI coordinates
% load([scriptDir, 'mri_targets_all_updated.mat']);
% load([scriptDir, 'mri_targets_all_wMNI.mat']);
load([scriptDir, 'mri_targets_all_wMNI_updated.mat']); % mri_targets_all is loaded here.

% Load FSL processed masks
anatDir = 'F:\Linux_Shared\CausCon\VisorT1_iso2nii\Selected\FSL_BET\';
anatFiles = dir([anatDir '*_VISOR_brain_outskin_mesh.nii.gz']);

missingMNIsubs = {'CAUSCON_1035_RC3679_VISOR';...
    'CAUSCON_2003_RC8393_VISOR';...
    'CAUSCON_2073_RC8466_VISOR';...
    'CAUSCON_4116_RC6080_VISOR';...
    'CAUSCON_4206_RC9999_VISOR'};

%% Load .asc matching table
match_asc = readtable('cc_target2dicom_asc_select.csv');
for k = 1:size(visor_folders,1)
    if ~strcmp(visor_folders(k).name,match_asc.foldername{k})
        error(['mismatched entry: ', num2str(k)]);
    elseif k == size(visor_folders,1)
        disp('all entries matched!');
    end
end

%% Loop through all available T1
scalp2cortex_BET = [];

%%
for ksub = 1:size(visor_folders,1)
    
    subjname = visor_folders(ksub).name;
    disp(['ksub ', num2str(ksub), ', ', subjname]);
    subjfolder = [visor_folders(ksub).folder, '\', subjname];
    
    strParts = regexp(subjname,'_','split');
    scalp2cortex_BET(ksub).subID = strParts{2};
    
    %%
    scalp2cortex_BET(ksub).targetnames = [];
    scalp2cortex_BET(ksub).loc = [];
    scalp2cortex_BET(ksub).dist = [];
    
    %% find FSL BET processed nifti files
    anat_found = 0;
    scalpFilename = [];
    brainEnvFilename = [];
    scalp = [];
    brainEnv = [];
    
    for k = 1:size(anatFiles,1)
        
        if strfind(anatFiles(k).name,['_', scalp2cortex_BET(ksub).subID, '_'])
            anat_found = anat_found + 1;
            
            if anat_found > 1
                anatFiles(k).name
                error('More than one set of anatomy files found');
            end
            
            scalpFilename = anatFiles(k).name;
            brainEnvFilename = strrep(anatFiles(k).name, 'outskin', 'overlay');
            
        end
        
    end
    
    if ~anat_found
        disp('No anatomy files found');
    end
    
    %% Find target files
    tfile_found = 0;
    kasc_row_select = [];
    kasc_col_select = [];
    MNICoord = [];
    
    % find subject target files based on first asc file
    for kasc_row = 1:size(mri_targets_all,1)
        
        if strfind(mri_targets_all{kasc_row,1}.name, ['_', scalp2cortex_BET(ksub).subID, '_'])
            tfile_found = tfile_found + 1;
            
            if tfile_found > 1
                disp(mri_targets_all{kasc_row,1}.name);
                error('More than one set of target files found');
            end
            
            kasc_row_select = kasc_row;
            
        end
        
        % find selected asc file based on matching table
        if ~isempty(kasc_row_select)
            if strcmp(visor_folders(ksub).name, 'CAUSCON_1027_RC2545_VISOR')
                %exception for concatenating sessions
                for kasc_col = 1:size(mri_targets_all,2)
                    if ~isempty(mri_targets_all{kasc_row_select,kasc_col})
                        if strcmp(mri_targets_all{kasc_row_select,kasc_col}.name, match_asc.ascfile1{ksub})
                            kasc1 = kasc_col;
                        elseif strcmp(mri_targets_all{kasc_row_select,kasc_col}.name, match_asc.ascfile2{ksub})
                            kasc2 = kasc_col;
                        elseif strcmp(mri_targets_all{kasc_row_select,kasc_col}.name, match_asc.ascfile3{ksub})
                            kasc3 = kasc_col;
                        end
                    end
                end
                
                targetnames = [mri_targets_all{kasc_row_select,kasc1}.targetnames;...
                    mri_targets_all{kasc_row_select,kasc2}.targetnames;...
                    mri_targets_all{kasc_row_select,kasc3}.targetnames];
                visorCoord = [mri_targets_all{kasc_row_select,kasc1}.pointNeuroVX;...
                    mri_targets_all{kasc_row_select,kasc2}.pointNeuroVX;...
                    mri_targets_all{kasc_row_select,kasc3}.pointNeuroVX];
                MNICoord = [mri_targets_all{kasc_row_select,kasc1}.coordMNI;...
                    mri_targets_all{kasc_row_select,kasc2}.coordMNI;...
                    mri_targets_all{kasc_row_select,kasc3}.coordMNI];
                
            else %one ascfile used
                kasc_col_select = 1; %if no specific match, take the only file
                if ~isempty(match_asc.ascfile1{ksub})
                    for kasc_col = 1:size(mri_targets_all,2)
                        if ~isempty(mri_targets_all{kasc_row_select,kasc_col})
                            if strcmp(mri_targets_all{kasc_row_select,kasc_col}.name, match_asc.ascfile1{ksub})
                                kasc_col_select = kasc_col;
                            end
                        end
                    end
                end
                
                targetnames = mri_targets_all{kasc_row_select,kasc_col_select}.targetnames;
                visorCoord = mri_targets_all{kasc_row_select,kasc_col_select}.pointNeuroVX;
                
                % skip if missing MNI transform
                if ~ismember(visor_folders(ksub).name, missingMNIsubs)
                    MNICoord = mri_targets_all{kasc_row_select,kasc_col_select}.coordMNI;
                end
            end
            
        end
        
    end
    
    if ~tfile_found
        disp('No associated visor target file found');
    end
    
    %% exception files
    exception = 0;
    %     skip_files = {'CAUSCON_2044_RC4167_VISOR','CAUSCON_2108_RC10494_VISOR',...
    %         'CAUSCON_4199_RC10455_VISOR'};
    skip_files = {'CAUSCON_1023_RC2324_VISOR', 'CAUSCON_3090_RC6292_VISOR'};
    
    for kskip = 1:size(skip_files,2)
        if strcmp(subjname, skip_files{kskip}) %exceptions
            exception = 1;
        end
    end
    
    if exception
        disp('Other problems found');
    end
    
    %% Process files with all information
    if ~(anat_found & tfile_found & ~exception)
        disp(['Skipping ', scalp2cortex_BET(ksub).subID]);
    else
        disp(['scalpFilename: ', scalpFilename]);
        disp(['brainEnvFilename: ', brainEnvFilename]);
        %         disp(['targetFilename: ', mri_targets_all{kfile_select,1}.name]);
        
        % 'load_untouch_nii' used instead of 'load_nii' due to rotation and
        % shearing in affine matrix of header
        scalp = load_untouch_nii([anatDir, scalpFilename]);
        brainEnv = load_untouch_nii([anatDir, brainEnvFilename]);
        
        %% Extract only the non-zero voxels of scalp
        scalp_X = [];
        scalp_Y = [];
        scalp_Z = [];
        for kz = 1:size(scalp.img,3)
            [X, Y] = find(scalp.img(:,:,kz)>0);
            Z = ones(length(X),1)*kz;
            scalp_X = [scalp_X; X];
            scalp_Y = [scalp_Y; Y];
            scalp_Z = [scalp_Z; Z];
        end
        
        numvox(ksub,1) = size(scalp.img,1);
        numvox(ksub,2) = size(scalp.img,2);
        numvox(ksub,3) = size(scalp.img,3);
        
        % Define size of voxels to estimate distance in mm
        % order was swapped in image, so has to be done here
        dimX = voxdims(ksub,2);
        dimY = voxdims(ksub,1);
        dimZ = voxdims(ksub,3);
        
        %% Extract only the non-zero voxels of brain envelop
        brainEnv_X = [];
        brainEnv_Y = [];
        brainEnv_Z = [];
        for kz = 1:size(brainEnv.img,3)
            [X, Y] = find(brainEnv.img(:,:,kz)>0);
            Z = ones(length(X),1)*kz;
            brainEnv_X = [brainEnv_X; X];
            brainEnv_Y = [brainEnv_Y; Y];
            brainEnv_Z = [brainEnv_Z; Z];
        end
        
        %%
        scalp2cortex_BET(ksub).targetnames = targetnames;
        targetLoc = cat(2,visorCoord,ones(size(visorCoord,1),1))';
        
        %% Loop through all targets in file
        for ktarget = 1:size(targetLoc,2)
            %     ktarget = 1;
            
            %% Find scalp location closest to the target coordinate
            
            mindist = 999; %an impossible large distance
            mindist_loc = [];
            
            for k = 1:length(scalp_X)
                currdist = sqrt(((scalp_X(k)-targetLoc(1,ktarget))*dimX).^2 + ...
                    ((scalp_Y(k)-targetLoc(2,ktarget))*dimY).^2 + ...
                    ((scalp_Z(k)-targetLoc(3,ktarget))*dimZ).^2);
                if currdist < mindist
                    mindist = currdist;
                    mindist_loc(1,:) = [scalp_X(k) scalp_Y(k) scalp_Z(k)];
                end
            end
            
            if mindist == 999
                error('no voxel found');
            end
            
            %% Find brain envelop location closest to the target coordinate
            mindist2 = 999; %an impossible large distance
            mindist2_loc = [];
            
            for k = 1:length(brainEnv_X)
                currdist = sqrt(((brainEnv_X(k)-targetLoc(1,ktarget))*dimX).^2 + ...
                    ((brainEnv_Y(k)-targetLoc(2,ktarget))*dimY).^2 + ...
                    ((brainEnv_Z(k)-targetLoc(3,ktarget))*dimZ).^2);
                if currdist < mindist2
                    mindist2 = currdist;
                    mindist2_loc(1,:) = [brainEnv_X(k) brainEnv_Y(k) brainEnv_Z(k)];
                end
            end
            
            if mindist2 == 999
                error('no voxel found');
            end
            
            %% Find scalp distance closest to brain envelop location
            
            mindist3 = 999; %an impossible large distance
            mindist3_loc = [];
            
            for k = 1:length(scalp_X)
                currdist = sqrt(((scalp_X(k)-mindist2_loc(1))*dimX).^2 + ...
                    ((scalp_Y(k)-mindist2_loc(2))*dimY).^2 + ...
                    ((scalp_Z(k)-mindist2_loc(3))*dimZ).^2);
                if currdist < mindist3
                    mindist3 = currdist;
                    mindist3_loc(1,:) = [scalp_X(k) scalp_Y(k) scalp_Z(k)];
                end
            end
            
            if mindist3 == 999
                error('no voxel found');
            end
            
            %% Update loop variable
            scalp2cortex_BET(ksub).targets = targetLoc(1:3,:)';
            scalp2cortex_BET(ksub).targetsMNI = MNICoord;
            scalp2cortex_BET(ksub).targetnames = targetnames;
            
            scalp2cortex_BET(ksub).loc(ktarget,:) = mindist3_loc;
            scalp2cortex_BET(ksub).dist(ktarget,1) = mindist3;
            
            scalp2cortex_BET(ksub).loc_scalp2target(ktarget,:) = mindist_loc;
            scalp2cortex_BET(ksub).dist_scalp2target(ktarget,1) = mindist;
            
            % this is the 
            scalp2cortex_BET(ksub).loc_brainEnv2target(ktarget,:) = mindist2_loc;
            scalp2cortex_BET(ksub).dist_brainEnv2target(ktarget,1) = mindist2;
            
            scalp2cortex_BET(ksub).loc_scalp2brainEnv(ktarget,:) = mindist3_loc;
            scalp2cortex_BET(ksub).dist_scalp2brainEnv(ktarget,1) = mindist3;
            
        end
    end
end

%% write selected coordinates to text file

outpath = 'F:\Linux_Shared\CausCon\VisorT1_iso2nii\trgtBrainEnv_MRIvox\';

fieldSelect = 'loc_brainEnv2target';
suffixStr = 'trgtBrainEnv';

for id = 1:size(scalp2cortex_BET,2)
    
    if ~isempty(scalp2cortex_BET(1,id).targetnames)
        curr_file = [outpath, 'CAUSCON_', scalp2cortex_BET(1,id).subID,'_', suffixStr,'.txt'];
        targetnames = scalp2cortex_BET(1,id).targetnames;
        eval(['X = scalp2cortex_BET(1,id).', fieldSelect, '(:,1);']);
        eval(['Y = scalp2cortex_BET(1,id).', fieldSelect, '(:,2);']);
        eval(['Z = scalp2cortex_BET(1,id).', fieldSelect, '(:,3);']);
        curr_table = table(targetnames,X,Y,Z);
        writetable(curr_table, curr_file);
    end
    
end

%% Cacluate corresponding MNI coordinates for each target file
% Done through FSL

%% Import MNI coordinates for each target into 'mri_targets_all' variable
MNI_coordMain = 'F:\Linux_Shared\CausCon\MNI_trgtBrainEnv\';
MNI_coordDirs = dir([MNI_coordMain, 'CAUSCON_*']);

suffixStr = 'trgtBrainEnv';
newfieldname = 'loc_brainEnv2targetMNI';

skip_files = {'CAUSCON_2003_trgtBrainEnv'};

MNIfile_found = nan(1,size(scalp2cortex_BET,2));

for id = 1:size(scalp2cortex_BET,2)

    if ~isempty(scalp2cortex_BET(1,id)) %entry is not empty
        MNIfile_found(1,id) = 0;
        fileID = ['CAUSCON_', scalp2cortex_BET(1,id).subID,'_', suffixStr];
        
        if ~ismember(fileID,skip_files) %entry not a skipped file
            % find corresponding MNI
            for kMNI = 1:size(MNI_coordDirs,1)
                
                if strcmp(fileID,MNI_coordDirs(kMNI).name)
                    MNIfile = [MNI_coordMain, MNI_coordDirs(kMNI).name, ...
                        '\', MNI_coordDirs(kMNI).name, '_MNI.txt']
                    T_coordMNI = readtable(MNIfile);
                    
                    eval(['scalp2cortex_BET(1,id).', newfieldname,...
                        ' = table2array(T_coordMNI(:,2:4));']);
                    
                    MNIfile_found(1,id) = 1;
                    
                end
            end
        end
    end
    
end
