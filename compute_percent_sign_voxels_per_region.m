
function compute_percent_sign_voxels_per_region(niak_path, outputs_path, output_name, atlas_name, atlas_filepath, atlas_info_filepath, p_filepath, p_value_sign, t_filepath)
%% Isabelle Lajoie 2023-09-26
% This function compute the number/percentage of significant voxels per region of the provided atlas,  
% based on the p-maps, and compute the significant t-value averaged per region.
%
% Inputs: 
%  niak_path - path toward niak
%  outputs_path - path for outputs. Will create it if does not exist 
%  output_name - name to add to the generated files' name
%  atlas_name - i) 'allen' : atlas_info_filepath must include "id_141" and "name" columns 
%               ii)'cerebra' : atlas_info_filepath must include "LH Label", "RH Label" and "Label Name" columns
%               iii) other atlas : atlas_info_filepath must include "ID" and "Name" columns
%  atlas_filepath - atlas mnc filepath
%  atlas_info_filepath - atlas csv filepath containing regions' ID and name
%  p_filepath - mnc file of p-value
%  p_value_sign - p value below which we consider significance 
%  t_filepath - mnc file of t-stats
%
% The function generates the following files (names starting with "percent_sign_vox"):
%  A csv file containing the computed stats
%  A mnc file with each atlas' region taking a value between 0 and 100,
%  representing the obtained percentage of significant voxels
%%
    addpath(genpath(niak_path));
    %addpath(genpath('/home/cic/lajisa/vandalab/niak'));

    function [ROI_nb_voxels, nb_sign_voxels, perc_sign_vox, mean_tval] = compute_ROI_sign_voxel(p_mask, roi_mask, sign_tmaps)
    % This function compute percent of significant voxels for a specific region
    % Inputs:
    %   p_mask - binary p mask (FDR corrected or not) representing significance
    %   roi_mask - binary ROI mask
    %   sign_tmaps - tmaps masked with p_mask
    % Outputs:
    %   ROI_nb_voxels - number of voxels in ROI
    %   nb_sign_voxels - number of significant voxels in the ROI
    %   perc_sign_vox - percentage of significant voxels in the ROI
    %   mean_tval - tval averaged over the significant voxels in the ROI
        ROI_nb_voxels = sum(roi_mask,'all');
        nb_sign_voxels = sum(p_mask(roi_mask),'all');
        perc_sign_vox = round((nb_sign_voxels / ROI_nb_voxels * 100),1);
        mean_tval = round((sum(sign_tmaps(roi_mask),'all')/nb_sign_voxels),2);
    end

    if ~exist(outputs_path, 'dir')
        mkdir(folderName);
    end

    [hdr,atlas]=niak_read_minc(atlas_filepath);
    [~,p]=niak_read_minc(p_filepath);
    [~,tmaps]=niak_read_minc(t_filepath);
    
    atlas_info=readtable(atlas_info_filepath,'Delimiter',',');
    
    % Output brain maps
    out_maps = atlas .* 0;
    out_maps_filename = fullfile(outputs_path, strcat('percent_sign_vox_', output_name, '_', atlas_name,'.mnc'));
    
    % Output csv file
    output_header = {'Atlas', 'ID', 'Name', 'Voxel_count', 'Nb_Sign_Vox', 'Perc_Sign_Vox', 'Mean_t_values'};
    output_filepath = fullfile(outputs_path, strcat('percent_sign_vox_', output_name, '_', atlas_name,'.csv'));
    
    % p mask based on the input p_filepath and p_value_sign
    p_mask = (p<p_value_sign) & (p>0);
    % Significant tmaps 
    sign_tmaps = p_mask.*tmaps;
   
    atlas_name = lower(atlas_name);
    if contains(atlas_name, 'allen');
        %% Allen atlas
        disp('Compute percent sign voxels for Allen atlas: will use "id_141" (left ID) and "name" columns. Right ID = "id_141" value + 1000.')
        output_array = cell(numel(atlas_info.id_141)*2, numel(output_header));
        for i = 1:numel(atlas_info.id_141)
            % Left and Right hemisphere
            L_roi_ID = atlas_info.id_141(i);
            R_roi_ID = atlas_info.id_141(i) + 1000;
            roi_ID = {L_roi_ID, R_roi_ID};
            roi_masks = {atlas==L_roi_ID, atlas==R_roi_ID};
            hemisphere_label = {'L','R'};

            for j = 1:numel(roi_masks)
                roi_mask = roi_masks{j};
                roi_name = strcat(hemisphere_label(j),'_',strrep(atlas_info.name(i),' ',''));
                
                [ROI_nb_voxels, nb_sign_voxel, perc_sign_vox, mean_tval] = compute_ROI_sign_voxel(p_mask, roi_mask, sign_tmaps);
                output_array((i-1)*2+j,:) = {atlas_name, roi_ID{j}, roi_name, ROI_nb_voxels, nb_sign_voxel, perc_sign_vox, mean_tval};
               
                out_maps(atlas==roi_ID{j}) = perc_sign_vox;
            end
        end
    elseif contains(atlas_name, 'cerebra');
        %% CerebrA atlas
        disp('Compute percent sign voxels for CerebrA atlas: will use "LH Label", "RH Label" and "Label Name" columns.')
        output_array = cell(numel(atlas_info.RHLabel)*2, numel(output_header));
        for i = 1:numel(atlas_info.RHLabel)
            % Left and Right hemisphere
            L_roi_ID = atlas_info.LHLabel(i);
            R_roi_ID = atlas_info.RHLabel(i);
            roi_ID = {L_roi_ID, R_roi_ID};
            roi_masks = {atlas==L_roi_ID, atlas==R_roi_ID};
            hemisphere_label = {'L','R'};
            for j = 1:numel(roi_masks)
                roi_mask = roi_masks{j};
                roi_name = strcat(hemisphere_label(j),'_',strrep(atlas_info.LabelName(i),' ',''));

                [ROI_nb_voxels, nb_sign_voxel, perc_sign_vox, mean_tval] = compute_ROI_sign_voxel(p_mask, roi_mask, sign_tmaps); 
                output_array((i-1)*2+j,:) = {atlas_name, roi_ID{j}, roi_name, ROI_nb_voxels, nb_sign_voxel, perc_sign_vox, mean_tval};
                
                out_maps(atlas==roi_ID{j}) = perc_sign_vox;
            end
        end
    else
        %% Others
        disp('Compute percent sign voxels for a generic atlas: will use "ID" and "Name" columns.')
        output_array = cell(numel(atlas_info.ID), numel(output_header));
        for i = 1:numel(atlas_info.ID)
            roi_ID = atlas_info.ID(i);
            roi_name = atlas_info.Name(i);
            roi_mask = atlas==roi_ID;

            [ROI_nb_voxels, nb_sign_voxel, perc_sign_vox, mean_tval] = compute_ROI_sign_voxel(p_mask, roi_mask, sign_tmaps); 
            output_array(i,:) = {atlas_name, roi_ID, roi_name, ROI_nb_voxels, nb_sign_voxel, perc_sign_vox, mean_tval};
            
            out_maps(atlas==roi_ID) = perc_sign_vox;
        end    
    end

    %% Sort by percentage of significant voxel (highest first) and save results
    output_table = array2table(output_array, 'VariableNames', output_header);
    sortedTable = sortrows(output_table, 'Perc_Sign_Vox', 'descend'); 
    writetable(sortedTable,	output_filepath)
    
    %% Save brain maps containing the percent of significant voxels per ROI
    hdr.file_name = out_maps_filename;
    niak_write_minc(hdr, out_maps);
end
