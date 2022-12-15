% extract scalp to distance information and export to csv file:
clear all;

input_file = 'data_subject_info/targetInfo_11sites_rmol.mat';
load(input_file);

result = {'site', 'subject', 'scalp_dist'};
for i = 1:length(targetInfo_rmol)
    curr = cell(targetInfo_rmol(i).num_sub, 3);
    curr(:,1) = {targetInfo_rmol(i).name};
    curr(:,2) = num2cell(targetInfo_rmol(i).sublist);
    curr(:,3) = num2cell(targetInfo_rmol(i).dist_scalp2brainEnv);
    
    result = [result; curr];
end

% Convert cell to a table and use first row as variable names
T = cell2table(result(2:end,:), 'VariableNames', result(1,:));
% Write the table to a CSV file
writetable(T, 'data_subject_info/dist_to_scalp.csv')


