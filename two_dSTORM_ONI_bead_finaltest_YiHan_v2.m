% v2: if the field has <2 NPs pass the test, delete the mat file. (03/01/2024)

addpath(genpath('M:\Yi-Han\ONI_Storm\Yi-Han_edited_code\MATLAB_library'))

% Call the main function to start the process
FWHMtest= start();

% clear; clc;

temp = 1;
for i = 1:length(FWHMtest)
    if FWHMtest(i).pass_FWHMthr_ratio_corr~=1
        % disp(FWHMtest(i).SML_file)
        % disp(strcat('    pass_FWHMthr_ratio_corr =', num2str(FWHMtest(i).pass_FWHMthr_ratio_corr)))
        if length(str2num(FWHMtest(i).pass_FWHMthr_par_num))>=2
            test_result = process1(FWHMtest(i).folder, FWHMtest(i).NP1_fif_file, FWHMtest(i).SML_file, str2num(FWHMtest(i).pass_FWHMthr_par_num));
            FWHMtest2(temp)= test_result;
            temp=temp+1;
            clear test_result
        elseif length(str2num(FWHMtest(i).pass_FWHMthr_par_num))<2
            temp = split(FWHMtest(i).SML_file, '.csv');
            driftcorr_matrix = strcat(FWHMtest(i).folder, '\output\', char(temp(1)), '_driftcorr.mat');
            if exist(driftcorr_matrix,"file") ==2
                delete(driftcorr_matrix)
            end  
        end
    end
end





function FWHMtest = start()
    % Get the list of files and nanoparticles
    [folder_name, file_list] = processTextAndFiles();
     
    % Loop through each file set and process
    for i = 1:length(file_list)
        file_set = file_list{i};
        csv_file = file_set{1};
        tif_file = file_set{2};
        % Convert cell array of strings to a numeric array
        select_NP = cellfun(@str2num, file_set{3});
        % Call the process function for each set
        test_result = process1(folder_name, tif_file, csv_file, select_NP);
        FWHMtest(i)= test_result;
        clear test_result
    end
end

% function [folder_name, file_list] = processTextAndFiles()
%     % Use UI to get the text file
%     [txtFileName, txtPath] = uigetfile('*.txt', 'Select the text file');
%     if isequal(txtFileName, 0)
%         error('No file selected');
%     end
%     txtFileFullPath = fullfile(txtPath, txtFileName);
% 
%     % Read and parse the text file
%     [fileID, errmsg] = fopen(txtFileFullPath, 'r');
%     if fileID < 0
%         error(errmsg);
%     end
%     data = textscan(fileID, '%s', 'Delimiter', '\n');
%     fclose(fileID);
%     lines = data{1};
% 
%     % Process each line
%     file_map = containers.Map;
%     for i = 1:length(lines)
%         line = lines{i};
%         % Splitting line into test name and numbers
%         parts = strsplit(line, ':');
%         if numel(parts) < 2
%             continue; % Skip lines that don't have a ':' character
%         end
%         test_name = strtrim(parts{1}); % Extract test name
%         nums = strtrim(parts{2}); % Extract number string
%         nums = strsplit(nums, ' '); % Split by space
%         nums = nums(~cellfun('isempty', nums)); % Remove any empty strings
% 
%         if length(nums) > 1 % Include line if more than one particle number
%             file_map(test_name) = nums;
%         end
%     end
% 
%     % Select folder
%     folder_name = uigetdir;
%     if folder_name == 0
%         error('No folder selected');
%     end
% 
%     % Filter TIFF and CSV files
%     tif_files = dir(fullfile(folder_name, '*_NP_f1.tif'));
%     csv_files = dir(fullfile(folder_name, '*_SMLs.csv'));
% 
%     % Create output list
%     file_list = {};
%     keys = file_map.keys;
%     for i = 1:length(keys)
%         key = keys{i};
%         tif_name = strcat(key, '_NP_f1.tif');
%         csv_name = strcat(key, '_SMLs.csv');
% 
%         % Check if both TIFF and CSV files exist
%         if any(strcmp({tif_files.name}, tif_name)) && any(strcmp({csv_files.name}, csv_name))
%             tif_full_path = fullfile(folder_name, tif_name);
%             csv_full_path = fullfile(folder_name, csv_name);
%             file_list{end+1} = {csv_full_path, tif_full_path, file_map(key)};
%         end
%     end
% 
%     % Display the file list in a structured format
%     disp('Selected files and PNGs:');
%     for i = 1:length(file_list)
%         disp(['Test Name: ', key]);
%         disp(['CSV File: ', file_list{i}{1}]);
%         disp(['TIFF File: ', file_list{i}{2}]);
%         selectedParticles = strjoin(file_list{i}{3}, ', '); % Remove leading comma
%         disp(['Selected Particles: ', selectedParticles]);
%         disp('-----------------------------------');
%     end
% end


% function process1(raw_tif_filename, SMLs_filename, select_NP )
%      if ~isnumeric(select_NP)
%         error('select_NP must be a numeric array.');
%      end
% 
%     % Add necessary paths (adjust as needed)
%     addpath(genpath('M:\Yen-Cheng\YCC_Matlab'));
%     disp(['Raw Tif File: ', raw_tif_filename]);
%     disp(['SML Filename: ', SMLs_filename]);
%     disp('Selected Nano-Particles:');
%     disp(select_NP);
% 
%     % Initialization
% 
%     offset = 0;
%     first_frame_NP = 11;
%     first_frame_NP_after_prebleach = 1011;
%     first_frame_GFP = 1;
%     start_frame = 1010; % GFP 0-9, pre-bleach 10-1009 frames
%     end_frame = 21009; % total 40k frames
%     start_frame = start_frame + offset;
%     end_frame = end_frame + offset;
%     minSMLs_driftcorr = (end_frame - start_frame + 1) * 0.7;
%     nm_per_pixel = 116.999998688698;
%     precision_xy_thr = 20;
%     precision_z_thr = inf;
%     GFP_loc_thr = 200;
%     GFP_SMLs_stdev_thr = 30;
%     GFP_SMLs_stdev_thr_z = 60;
%     FWHMthr = 20;
%     search_radius_nm = 200; % particle size = 200 nm
%     runDBSCAN = false;
%     DBSCAN_eps_1 = 15; % 20 nm cluster size
%     if_2D = true;
%     if_outputimage = true;
%     clus_max_axis_thr = inf;
%     clu1_counts = 0;
%     doparse = true;
% 
%      % Adjusting TopDir and Mypath
%     [TopDir, name, ext] = fileparts(raw_tif_filename); % Extract the directory and filename separately
%     raw_tif_file = [name, ext]; % Reconstruct filename with extension
% 
% 
%     matlaboutbase = '\output\';
%     disp(['TopDir: ', TopDir]);
%     disp(['matlaboutbase: ', matlaboutbase]);
%     mkdir(strcat(TopDir, matlaboutbase, 'image_output'));
%     addpath(genpath(TopDir)); % add the specified folder and subfolders to the path.
% 
%     % Remove .csv extension from SMLs_filename if it exists
%     [~, name, ~] = fileparts(SMLs_filename); % Extract the name without extension
%     SMLs_filename = name;
% 
%     parse_particle_filename = strcat(SMLs_filename, '_NP_Particles');
%     delete_NP=[ ...
%            ];
% 
%     % Load particle info (adjust as needed)
%     particleFilePath = strcat(TopDir, matlaboutbase, parse_particle_filename, '.mat');
%     if exist(particleFilePath, 'file')
%         load(particleFilePath);
%     else
%         error('Particle file not found: %s', particleFilePath);
%     end
% 
%     % Load NP low res image
%     Im_NP = imread(fullfile(TopDir, strcat('*',raw_tif_file))); % Corrected path
%     % imwrite(Im_NP,strcat(TopDir,raw_tif_filename,'NP_f1.tiff'))
%     % tfilebase_all={'20190917 B-';...
%     %                '20190919 B-'};
% 
% 
% 
%     % output x,y,frame
%     output_ini_f=start_frame;
%     output_fin_f=end_frame;
%     output_total_f=end_frame-start_frame+1;
%     frame=(output_ini_f:1:output_fin_f)';
%     logic_frame=zeros((end_frame-start_frame+1),1);
% 
%     output_pixx=[];
%     output_pixy=[];
%     output_pixx_map=[];
%     output_pixy_map=[];
%     output_xlimspix=[];
%     output_ylimspix=[];
%     output_par_num=[];
% 
%     before_prebleach_pix=[];
% 
%     win_size=20; %x frames average
% 
%     for par=1:size(AllParticles,1)
% 
%     %     if sum(delete_NP==AllParticles(par).partnum)>0
%     %         continue
%     %     end
% 
%         if sum(select_NP==AllParticles(par).partnum)==0
%             continue
%         end
% 
% 
%         toutput_x=nan((end_frame-start_frame+1),1);
%         toutput_y=nan((end_frame-start_frame+1),1);
%         toutput_x_map=nan((end_frame-start_frame+1),1);
%         toutput_y_map=nan((end_frame-start_frame+1),1);
% 
%         tframe=AllParticles(par).SMLLabel.Frame;
% 
%         t_tiff_pix=AllParticles(par).partlocpix;
% 
%         txRaw=AllParticles(par).SMLLabel.XRaw_pix_;
%         tyRaw=AllParticles(par).SMLLabel.YRaw_pix_;
% 
%         tx=AllParticles(par).SMLLabel.X_pix_;
%         ty=AllParticles(par).SMLLabel.Y_pix_;
% 
%         if size(tx,1)<minSMLs_driftcorr
%             continue
%         end
% 
%         for i=1:size(tframe,1)
%             if (tframe(i))<(start_frame) ||... %remove SMLs from prebleaching time
%                   (tframe(i))>end_frame      %remove NP SMLs from GFP illumination time
%                 continue
%             end
%             logic_frame(tframe(i)+1-start_frame)=1;
%             toutput_x(tframe(i)+1-start_frame)=tx(i);
%             toutput_y(tframe(i)+1-start_frame)=ty(i);
%             toutput_x_map(tframe(i)+1-start_frame)=tx(i)-t_tiff_pix(1);
%             toutput_y_map(tframe(i)+1-start_frame)=ty(i)-t_tiff_pix(2);
%         end
% 
%         output_pixx=cat(2,output_pixx,toutput_x);
%         output_pixy=cat(2,output_pixy,toutput_y);
% 
%         output_pixx_map=cat(2,output_pixx_map,nanmean(toutput_x_map));
%         output_pixy_map=cat(2,output_pixy_map,nanmean(toutput_y_map));
% 
%         output_xlimspix=cat(1,output_xlimspix,AllParticles(par).xlimspix);
%         output_ylimspix=cat(1,output_ylimspix,AllParticles(par).ylimspix);
% 
%         output_par_num=cat(1,output_par_num,AllParticles(par).partnum);
% 
%         before_prebleach_pix=cat(1,before_prebleach_pix,t_tiff_pix);
%     end
%     center_pixx=nanmean(output_pixx,1);
%     center_pixy=nanmean(output_pixy,1);
%     local_output_pixx=output_pixx-ones(size(output_pixx,1),1)*center_pixx;
%     local_output_pixy=output_pixy-ones(size(output_pixy,1),1)*center_pixy;
%     mean_pixx=nanmean(local_output_pixx,2);
%     mean_pixy=nanmean(local_output_pixy,2);
% 
%     mean_pixx_map=output_pixx_map;
%     mean_pixy_map=output_pixy_map;
% 
%     plot_merge_drift=figure;
%     meanx_nm=mean_pixx*nm_per_pixel;
%     meany_nm=mean_pixy*nm_per_pixel;
%     pixel_thispart_drift=[meanx_nm, meany_nm];
%     c = linspace(1,10,length(pixel_thispart_drift(:,1)));
%     ax2_drift = axes;
%     plot_clus_drift=scatter(ax2_drift,pixel_thispart_drift(:,1), pixel_thispart_drift(:,2), 1, c, 'filled');
%     xlabel('nm')
%     ylabel('nm')
%     hold on
%     colormap(ax2_drift,'jet')
%     OutDir_partnum_drift = strcat(TopDir, matlaboutbase, 'image_output\', SMLs_filename);
%     saveas(plot_merge_drift, strcat(OutDir_partnum_drift, '_NP_drift.png'));
%     close(plot_merge_drift);
% 
%     fixed_win_x=zeros(floor(output_total_f/win_size),1);
%     fixed_win_y=zeros(floor(output_total_f/win_size),1);
% 
%     %calculate the drift from the initial xy
%     for frame_win=1:floor(output_total_f/win_size)
%         fixed_win_x(frame_win)=nanmean(mean_pixx((frame_win-1)*win_size+1:(frame_win)*win_size));
%         fixed_win_y(frame_win)=nanmean(mean_pixy((frame_win-1)*win_size+1:(frame_win)*win_size));
%     end
%     ini_x=fixed_win_x(1);
%     ini_y=fixed_win_y(1);
% 
%     %drift correction
%     driftcorr_matrix_x=nan(size(output_pixx,1),1);
%     driftcorr_matrix_y=nan(size(output_pixy,1),1);
%     driftcorr_output_pixx=nan(size(output_pixx));
%     driftcorr_output_pixy=nan(size(output_pixy));
%     for frame_win=1:floor(output_total_f/win_size)
%         tdriftx=fixed_win_x(frame_win)-ini_x;
%         tdrifty=fixed_win_y(frame_win)-ini_y;
% 
%         driftcorr_output_pixx((frame_win-1)*win_size+1:(frame_win)*win_size,:)=...
%             output_pixx((frame_win-1)*win_size+1:(frame_win)*win_size,:)-tdriftx*ones(win_size,size(output_pixx,2));
%         driftcorr_output_pixy((frame_win-1)*win_size+1:(frame_win)*win_size,:)=...
%             output_pixy((frame_win-1)*win_size+1:(frame_win)*win_size,:)-tdrifty*ones(win_size,size(output_pixy,2));
% 
%         driftcorr_matrix_x((frame_win-1)*win_size+1:(frame_win)*win_size,1)= -tdriftx*ones(win_size,1);
%         driftcorr_matrix_y((frame_win-1)*win_size+1:(frame_win)*win_size,1)= -tdrifty*ones(win_size,1);
% 
%     end
%     %remainder from the win_size
%     if rem(output_total_f,win_size)~=0
%       driftcorr_output_pixx((frame_win)*win_size:end,:)=...
%             output_pixx((frame_win)*win_size:end,:)-tdriftx*ones(rem(output_total_f,win_size),size(output_pixx,2));
%       driftcorr_output_pixy((frame_win)*win_size:end,:)=...
%             output_pixy((frame_win)*win_size:end,:)-tdrifty*ones(rem(output_total_f,win_size),size(output_pixy,2));
% 
%       driftcorr_matrix_x((frame_win)*win_size:end,:)= -tdriftx*ones(rem(output_total_f,win_size),1);
%       driftcorr_matrix_y((frame_win)*win_size:end,:)= -tdrifty*ones(rem(output_total_f,win_size),1);
%     end
% 
%     driftcorr_matrix=[driftcorr_matrix_x driftcorr_matrix_y];
%     save(strcat(TopDir, matlaboutbase, SMLs_filename, '_driftcorr.mat'), 'driftcorr_matrix');
% 
%     total_NP_num=0;
%     pass_FWHMthr_num_ori=0;
%     pass_FWHMthr_num_corr=0;
% 
%     Im_low=min(min(Im_NP));
%     Im_high=max(max(Im_NP))*0.4;
% 
%     %%NP loc from Image after prebleach
%     % Parse SMLs with NP presenting; Peak finder to locate centroids in GFP image
%     prebleach_tif_file = strrep(raw_tif_file, '_NP_f1.tif', '_NP_f1001.tif');
%     Im_NP_after_prebleach = imread(strcat(TopDir, '\', prebleach_tif_file));
%     peaksizepix = 5;
%     peakintensitythresh =300;
%     [centr_lowres_NP_after_prebleach] = one_channel_partID_includeSide(Im_NP_after_prebleach, peaksizepix, peakintensitythresh, []);
% 
%     % NP peaks
%     % centr_lowres is
%     % [x,y,meanback,sumsig,corrsumsig,SNR,resolvpow,perc_backpix]
%     centr_lowres_NP_filt_after_prebleach = centr_lowres_NP_after_prebleach;
%     % filter resolvpow>0.8
%     centr_lowres_NP_filt_after_prebleach = centr_lowres_NP_filt_after_prebleach(centr_lowres_NP_filt_after_prebleach(:,7)>0.8,:);
%     % filter perc_backpix>50%
%     centr_lowres_NP_filt_after_prebleach = centr_lowres_NP_filt_after_prebleach(centr_lowres_NP_filt_after_prebleach(:,8)>50,:);
%     % filter SNR>1.5
%     centr_lowres_NP_filt_after_prebleach = centr_lowres_NP_filt_after_prebleach(centr_lowres_NP_filt_after_prebleach(:,6)>5,:);
% 
%     NPpix_after_prebleach = centr_lowres_NP_filt_after_prebleach(:,[1,2]);
% 
%     prebleach_driftcorr_x=[];
%     prebleach_driftcorr_y=[];
%         for par=1:size(output_pixx,2)
% 
% 
% 
%             %% correction of drift in prebleach
%             tsearch_pixx=before_prebleach_pix(par,1)+[-1,1]*search_radius_nm/nm_per_pixel;
%             tsearch_pixy=before_prebleach_pix(par,2)+[-1,1]*search_radius_nm/nm_per_pixel;
% 
%             tpos=NPpix_after_prebleach((NPpix_after_prebleach(:,1) > tsearch_pixx(1)...
%                                       & NPpix_after_prebleach(:,1) < tsearch_pixx(2)...
%                                       & NPpix_after_prebleach(:,2) > tsearch_pixy(1)...
%                                       & NPpix_after_prebleach(:,2) < tsearch_pixy(2)),:);
% 
%             if size(tpos,1)==0
%                 disp('no NP after prebleach step'); continue
%             elseif size(tpos,1)>1
%                 disp('too dense of NP in search field'); continue
%             end
% 
%             prebleach_driftcorr_x=cat(2,prebleach_driftcorr_x,tpos(1)-before_prebleach_pix(par,1));
%             prebleach_driftcorr_y=cat(2,prebleach_driftcorr_y,tpos(2)-before_prebleach_pix(par,2));
%         end
%             prebleach_driftcorr_x
%             prebleach_driftcorr_y
% 
%            %% image output
% 
%         for par=1:size(output_pixx,2)
%             total_NP_num=total_NP_num+1;
% 
%             if if_outputimage
% 
%             %scatter + low res image
%             plot_merge_ori=figure;
%             pixel_thispart_ori=[output_pixx(:,par)-mean_pixx_map(par), output_pixy(:,par)-mean_pixy_map(par)];
%     %          pixel_thispart_ori=[output_pixx(:,par), output_pixy(:,par)];
% 
%             ax1 = axes;
%             plot_SMLs=scatter(ax1,pixel_thispart_ori(:,1), pixel_thispart_ori(:,2), 20, [1,0,0], '.');
%             hold on
%             colormap(ax1,'jet')
%             scale_bar_coord_ori=[output_xlimspix(par,2)-120/nm_per_pixel output_ylimspix(par,2)-20/nm_per_pixel;...
%                 output_xlimspix(par,2)-20/nm_per_pixel output_ylimspix(par,2)-20/nm_per_pixel];
%             plot_scalebar=plot(scale_bar_coord_ori(1:2,1),scale_bar_coord_ori(1:2,2),'-c', 'LineWidth', 3);hold on
%     %             pixel_clus = cell(size(clus1,1),1);
%     %             for k = 1: size(clus1,1)
%     %                 pixel_clus = clus1{k,1}.Points*transform_matrix+shift_matrix;
%     %                 plot_clus=plot(alphaShape(pixel_clus(:,1),pixel_clus(:,2)),'FaceAlpha',0.1,'FaceColor',[1,0,0],'EdgeAlpha',0);
%     %                 hold on
%     %             end
% 
%             plot_partnum = imshow(Im_NP, [Im_low, Im_high*0.7]);
%             axisgfp = gca;
%             axisgfp.XLim = output_xlimspix(par,:);
%             axisgfp.YLim = output_ylimspix(par,:);
%             axisgfp.Visible = 'off';
%             hold off
%             alpha(0.5)
% 
%             OutDir_partnum_ori = strcat(TopDir, matlaboutbase, 'image_output\', SMLs_filename, ...
%                 '-part-', num2str(output_par_num(par)));
%             saveas(plot_merge_ori, strcat(OutDir_partnum_ori, '_NP_merge_ori.png'));
%             close(plot_merge_ori);
%             end
% 
% 
% 
%             t_pixx=(output_pixx(:,par)-center_pixx(par))*nm_per_pixel;
%             t_pixy=(output_pixy(:,par)-center_pixy(par))*nm_per_pixel;
%             pdx=fitdist(t_pixx,'Normal');
%             pdy=fitdist(t_pixy,'Normal');
% 
%             pdx_fwhm= pdx.sigma*2*sqrt(2*log(2));
%     %             pdx_ci = paramci(pdx);
%     %             pdx_fwhmCI=round(pdx_ci(:,2)'*2*sqrt(2*log(2)),1);
% 
%             pdy_fwhm= pdy.sigma*2*sqrt(2*log(2));
%     %             pdy_ci = paramci(pdy);
%     %             pdy_fwhmCI=round(pdy_ci(:,2)'*2*sqrt(2*log(2)),1);
%             if pdx_fwhm<FWHMthr && pdy_fwhm<FWHMthr
%                pass_FWHMthr_num_ori=pass_FWHMthr_num_ori+1; 
%             end
% 
%             if if_outputimage
%             % histogram of fiducial markers
%             plot_hist_ori = figure;
%             plot_hist_ori.Units = 'inches';
%             plot_hist_ori.Position = [1 1 4 4];
%             plot_hist_ori.PaperPosition = [0 0 4 4];
% 
%             ax3_ori_x = axes('Position', [0.15 0.58 0.8 0.4]);
%             p_histx=histogram(t_pixx,40,'Normalization','pdf','FaceColor',[.5 .5 .5]);hold on
%     %             xlabel('nm');
%             ylabel('x pdf');
%             str={'FWHM',strcat(num2str(round(pdx_fwhm,1)),'nm')};
%     %                 strcat('CI: ',num2str(pdx_fwhmCI(1)),'-',num2str(pdx_fwhmCI(2)))};
%             peak_ctr_x=max(p_histx.BinCounts)/sum(p_histx.BinCounts)/p_histx.BinWidth;
%             text(pdx.mu,peak_ctr_x-0.01,str,'Color','red','FontSize',12);
%             xgrid = linspace(p_histx.BinEdges(1),p_histx.BinEdges(end),100)';
%             pdfx_Est = pdf(pdx,xgrid);
%             line(xgrid,pdfx_Est)
% 
%             ax3_ori_y = axes('Position', [0.15 0.12 0.8 0.4]);
%             p_histy=histogram(t_pixy,40,'Normalization','pdf','FaceColor',[.5 .5 .5]);hold on
%             xlabel('nm');
%             ylabel('y pdf');
%             str={'FWHM',strcat(num2str(round(pdy_fwhm,1)),'nm')};
%     %                strcat('CI: ',num2str(pdy_fwhmCI(1)),'-',num2str(pdy_fwhmCI(2)))};
%             peak_ctr_y=max(p_histy.BinCounts)/sum(p_histy.BinCounts)/p_histy.BinWidth;
%             text(pdy.mu,peak_ctr_y-0.01,str,'Color','red','FontSize',12);
%             ygrid = linspace(p_histy.BinEdges(1),p_histy.BinEdges(end),100)';
%             pdfy_Est = pdf(pdy,ygrid);
%             line(ygrid,pdfy_Est)
% 
%             OutDir_partnum_ori = strcat(TopDir, matlaboutbase, 'image_output\', SMLs_filename, ...
%                 '-part-', num2str(output_par_num(par)));
%             saveas(plot_hist_ori, strcat(OutDir_partnum_ori, '_NP_hist_ori.png'));
%             close(plot_hist_ori);
%             end
% 
% 
%             if if_outputimage
%             %scatter + low res image
%             plot_merge_corr=figure;
%     %         pixel_thispart_corr=[driftcorr_output_pixx(:,par), driftcorr_output_pixy(:,par)];
%              pixel_thispart_corr=[driftcorr_output_pixx(:,par)-mean_pixx_map(par), driftcorr_output_pixy(:,par)-mean_pixy_map(par)];
% 
%             ax1 = axes;
%             plot_SMLs=scatter(ax1,pixel_thispart_corr(:,1), pixel_thispart_corr(:,2), 20, [1,0,0], '.');
%             hold on
%             colormap(ax1,'jet')
%             scale_bar_coord_corr=[output_xlimspix(par,2)-120/nm_per_pixel output_ylimspix(par,2)-20/nm_per_pixel;...
%                 output_xlimspix(par,2)-20/nm_per_pixel output_ylimspix(par,2)-20/nm_per_pixel];
%             plot_scalebar=plot(scale_bar_coord_corr(1:2,1),scale_bar_coord_corr(1:2,2),'-c', 'LineWidth', 3);hold on
%     %             pixel_clus = cell(size(clus1,1),1);
%     %             for k = 1: size(clus1,1)
%     %                 pixel_clus = clus1{k,1}.Points*transform_matrix+shift_matrix;
%     %                 plot_clus=plot(alphaShape(pixel_clus(:,1),pixel_clus(:,2)),'FaceAlpha',0.1,'FaceColor',[1,0,0],'EdgeAlpha',0);
%     %                 hold on
%     %             end
% 
%             plot_partnum = imshow(Im_NP, [Im_low, Im_high*0.7]);
%             axisgfp = gca;
%             axisgfp.XLim = output_xlimspix(par,:);
%             axisgfp.YLim = output_ylimspix(par,:);
%             axisgfp.Visible = 'off';
%             hold off
%             alpha(0.5)
% 
%             OutDir_partnum_corr = strcat(TopDir, matlaboutbase, 'image_output\', SMLs_filename, ...
%                 '-part-', num2str(output_par_num(par)));
%             saveas(plot_merge_corr, strcat(OutDir_partnum_corr, '_NP_merge_corr.png'));
%             close(plot_merge_corr);
%             end
% 
% 
% 
% 
%             t_pixx=(driftcorr_output_pixx(:,par)-center_pixx(par))*nm_per_pixel;
%             t_pixy=(driftcorr_output_pixy(:,par)-center_pixy(par))*nm_per_pixel;
%             pdx=fitdist(t_pixx,'Normal');
%             pdy=fitdist(t_pixy,'Normal');
% 
%             pdx_fwhm= pdx.sigma*2*sqrt(2*log(2));
%     %             pdx_ci = paramci(pdx);
%     %             pdx_fwhmCI=round(pdx_ci(:,2)'*2*sqrt(2*log(2)),1);
% 
%             pdy_fwhm= pdy.sigma*2*sqrt(2*log(2));
%     %             pdy_ci = paramci(pdy);
%     %             pdy_fwhmCI=round(pdy_ci(:,2)'*2*sqrt(2*log(2)),1);
%             if pdx_fwhm<FWHMthr && pdy_fwhm<FWHMthr
%                pass_FWHMthr_num_corr=pass_FWHMthr_num_corr+1; 
%             end
% 
%             if if_outputimage
%             % histogram of fiducial markers
%             plot_hist_corr = figure;
%             plot_hist_corr.Units = 'inches';
%             plot_hist_corr.Position = [1 1 4 4];
%             plot_hist_corr.PaperPosition = [0 0 4 4];
% 
%             ax3_corr_x = axes('Position', [0.15 0.58 0.8 0.4]);
%             p_histx=histogram(t_pixx,40,'Normalization','pdf','FaceColor',[.5 .5 .5]);hold on
%     %             xlabel('nm');
%             ylabel('x pdf');
%             str={'FWHM',strcat(num2str(round(pdx_fwhm,1)),'nm')};
%     %                 strcat('CI: ',num2str(pdx_fwhmCI(1)),'-',num2str(pdx_fwhmCI(2)))};
%             peak_ctr_x=max(p_histx.BinCounts)/sum(p_histx.BinCounts)/p_histx.BinWidth;
%             text(pdx.mu,peak_ctr_x-0.01,str,'Color','red','FontSize',12);
%             xgrid = linspace(p_histx.BinEdges(1),p_histx.BinEdges(end),100)';
%             pdfx_Est = pdf(pdx,xgrid);
%             line(xgrid,pdfx_Est)
% 
%             ax3_corr_y = axes('Position', [0.15 0.12 0.8 0.4]);
%             p_histy=histogram(t_pixy,40,'Normalization','pdf','FaceColor',[.5 .5 .5]);hold on
%             xlabel('nm');
%             ylabel('y pdf');
%             str={'FWHM',strcat(num2str(round(pdy_fwhm,1)),'nm')};
%     %                strcat('CI: ',num2str(pdy_fwhmCI(1)),'-',num2str(pdy_fwhmCI(2)))};
%             peak_ctr_y=max(p_histy.BinCounts)/sum(p_histy.BinCounts)/p_histy.BinWidth;
%             text(pdy.mu,peak_ctr_y-0.01,str,'Color','red','FontSize',12);
%             ygrid = linspace(p_histy.BinEdges(1),p_histy.BinEdges(end),100)';
%             pdfy_Est = pdf(pdy,ygrid);
%             line(ygrid,pdfy_Est)
% 
%             OutDir_partnum_corr = strcat(TopDir, matlaboutbase, 'image_output\', SMLs_filename, ...
%                 '-part-', num2str(output_par_num(par)));
%             saveas(plot_hist_corr, strcat(OutDir_partnum_corr, '_NP_hist_corr.png'));
%             close(plot_hist_corr);
%             end
% 
% 
%         end 
% 
% 
%           disp(strcat('prebleach_driftcorr_x = ', num2str(prebleach_driftcorr_x)))
%           disp(strcat('prebleach_driftcorr_y = ', num2str(prebleach_driftcorr_y)))
%           disp(strcat('pass_FWHMthr_num_ori/total_NP_num = ', num2str(pass_FWHMthr_num_ori/total_NP_num)))
%           disp(strcat('pass_FWHMthr_num_corr/total_NP_num = ', num2str(pass_FWHMthr_num_corr/total_NP_num)))
%           disp(strcat('pass_FWHMthr_num_corr = ', num2str(pass_FWHMthr_num_corr)))
% 
% 
%     end

