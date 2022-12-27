cfg_file = 'Ind_Target_Node_11sites.mat';
load(cfg_file);

% EC.nod.color_thresh = 5;
% EC.nod.color_map_high = 10;
% EC.nod.color_map_type = 2;
EC.nod.CM(1, :) = [0, 0, 204]/255;
EC.nod.CM(2, :) = [255, 51, 51]/255;
% EC.nod.CM(21, :) = [255, 51, 51]/255;

save('tmp_cfg.mat', 'EC');

files = dir('node_mni*.node');
for i = 1: length(files)
    f_name = files(i).name;
    fig_name = [f_name(1:end-5), '.png'];
    
    BrainNet_MapCfg('BrainMesh_ICBM152.nv', f_name, ...
                    'tmp_cfg.mat', fig_name);
end