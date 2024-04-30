addpath(genpath('M:\Yen-Cheng\YCC_Matlab'));

colors = [ [0,0,0]; [1,0,0]; [0,1,0]; [0,0,1]; [1,0,1]; [0,1,1];
           [0,0,0]; [1,0,0]; [0,1,0]; [0,0,1]; [1,0,1]; [0,1,1] ];
% overlay factors
clc;clear;

nm_per_pixel=116.999998688698;
%total_pixel=196;
%x_off=100;
%y_off=0; 
precision_xy_thr=20;
precision_z_thr=inf;
GFP_loc_thr=200;
GFP_SMLs_stdev_thr=30;
GFP_SMLs_stdev_thr_z=60;

search_radius_nm = 200; %particle size = 200 nm

%minSMLs_2 = 5;
runDBSCAN = false;
DBSCAN_eps_1 = 15; %20 nm cluster size
%DBSCAN_eps_2 = 20; %30 nm cluster size
if_2D = true;
if_outputimage = true;
clus_max_axis_thr=inf;
clu1_counts=0;

doparse = true;
first_frame_NP=11;
start_frame_NP=1011;
end_frame_NP=21010;
offset=0;
start_frame_NP=start_frame_NP+offset;
end_frame_NP=end_frame_NP+offset;
mid_frame_NP=(end_frame_NP-start_frame_NP+1)/2+start_frame_NP;
minimal_SMLs_NP=(end_frame_NP-start_frame_NP+1)*0.9;

[raw_tif_filename, Mypath] = uigetfile({'*.tif; *.tiff';'*.xlsx';'*.mat';'*.*'}, 'Select NP_f1.tif files', 'MultiSelect','on');

TopDir = Mypath;
matlaboutbase = 'output\';
mkdir(strcat(Mypath, matlaboutbase, 'image_output'));
addpath(genpath(Mypath)); % add the specified folder and subfolders to the path.
tfilebase_all={'bead121720_beadSML';...  %% why do we need this?
%               '20190905 A-';...
%                '20191001 A-';...
%                '20191010 A-';...
               };                       

% tfilebase_all={'20190917 B-';...
%                '20190919 B-'};


% Files to input
tStart = tic;   %Stopwatch start

for name=1:size(tfilebase_all,1)
    
    tfilebase=tfilebase_all{name,:};
    for acq=1:length(raw_tif_filename) 
        %read low res images
        % first frame: NP, second frame: GFP
        
        temp = split(raw_tif_filename(acq), '_NP_f1.tif'); 
        SMLs_filename = strcat(char(temp(1)), '_SMLs');

        Im_NP = imread(strcat(TopDir,char(raw_tif_filename(acq))));
%         imwrite(Im_NP,strcat(TopDir,raw_tif_filename,'NP_f1.tiff'))
        %Im_NP = imread(raw_tif_filename,2);
        
        if exist(strcat(TopDir,char(raw_tif_filename(acq))), 'file') == 0
          % File does not exist
          % Skip to bottom of loop and continue with the loop
          continue;
        end
        [maxY,maxX]=size(Im_NP);

        %% Parse SMLs with NP presenting; Peak finder to locate centroids in NP image
        peaksizepix = 3;
        peakintensitythresh =500;
%        [centr_lowres] = one_channel_partID(Im_GFP, peaksizepix, peakintensitythresh, []);
        [centr_lowres_NP] = one_channel_partID_includeSide(Im_NP, peaksizepix, peakintensitythresh, []);
            % one_channel_partID_includeSide[Image, spot_size, spot_threshold, spot_coordinate].
            % spot_coordinate could be empty. 
            % output would be [ x, y, mean background, Sum of particle
            % intensity, particle int - background, SNR, resolvpow,
            % perc_backpix]

        %% NP peaks
        % centr_lowres is
        % [x,y,meanback,sumsig,corrsumsig,SNR,resolvpow,perc_backpix]
        centr_lowres_NP_filt = centr_lowres_NP;
        % filter resolvpow>0.8
        centr_lowres_NP_filt = centr_lowres_NP_filt(centr_lowres_NP_filt(:,7)>0.8,:);
        % filter perc_backpix>50%
        centr_lowres_NP_filt = centr_lowres_NP_filt(centr_lowres_NP_filt(:,8)>50,:);
        % filter SNR>1.5
        centr_lowres_NP_filt = centr_lowres_NP_filt(centr_lowres_NP_filt(:,6)>3,:);

        %% Filter GFP peaks
%         % centr_lowres is
%         % [x,y,meanback,sumsig,corrsumsig,SNR,resolvpow,perc_backpix]
%         centr_lowres_filt = centr_lowres;
%         % filter resolvpow>0.8
%         centr_lowres_filt = centr_lowres_filt(centr_lowres_filt(:,7)>0.8,:);
%         % filter perc_backpix>50%
%         centr_lowres_filt = centr_lowres_filt(centr_lowres_filt(:,8)>50,:);
%         % filter SNR>1.5
%         centr_lowres_filt = centr_lowres_filt(centr_lowres_filt(:,6)>1.5,:);
%         % remove overlap NP and GFP (distance<5)
%         k=1;
%         GFP_NP_logic=true(size(centr_lowres_filt(:,1),1),1);
%         for i=1:size(centr_lowres_filt(:,1))
%             GFP_NP_logic(i) = true;
%             for m=1:size(centr_lowres_NP_filt(:,1))
%                 d_NP_GFP=norm([centr_lowres_filt(i,1),centr_lowres_filt(i,2)]-[centr_lowres_NP_filt(m,1),centr_lowres_NP_filt(m,2)]);
%                if d_NP_GFP <5 && d_NP_GFP >0.5
%                   GFP_NP_logic(i)=false;
%                end
%             end
%         end
%         centr_lowres_filt2 = centr_lowres_filt(GFP_NP_logic,:);
%         clear GFP_NP_logic


        %% Create peak detection and filtering plots
        Im_low=min(min(Im_NP));
        Im_high=max(max(Im_NP))*0.75;
        pfil = plot_spot_filtering(centr_lowres_NP, centr_lowres_NP_filt);
        pdet = plot_spot_detection(Im_NP, centr_lowres_NP, centr_lowres_NP_filt, Im_low, Im_high);


        %% Check flag to parse ASCII file
        % True - save/close peak filtering/detection images and parse ASCII
        % False - end script
        OutDir= strcat(TopDir, matlaboutbase, SMLs_filename);
        if doparse
            saveas(pfil, strcat(OutDir, 'NP_part-filtering.png'));
            close(pfil);
            saveas(pdet, strcat(OutDir, 'NP_part-detection.png'));
            close(pdet);
%             saveas(pdet_no_overlap, strcat(OutDir, '_part-detection-noNPoverlap.png'));
%             close(pdet_no_overlap);
        else
            continue;
        end

        %%
        %read SMLs from csv file
%         ssds = tabularTextDatastore(SMLs_filename,'SelectedVariableNames',{'XRaw_pix_','YRaw_pix_'},'ReadSize','file');
%         SMLs_info_rawXY = read(ssds);
        Fi_SML = strcat(TopDir, SMLs_filename,'.csv');
        ds_SML = datastore(Fi_SML);
        
        AllParticles = [];       
        for partnum = 1:size(centr_lowres_NP_filt,1)
            XRaw_Right=centr_lowres_NP_filt(partnum,1)+ search_radius_nm/nm_per_pixel;
            XRaw_Left=centr_lowres_NP_filt(partnum,1)- search_radius_nm/nm_per_pixel;
            YRaw_Down=centr_lowres_NP_filt(partnum,2)+ search_radius_nm/nm_per_pixel;
            YRaw_Up=centr_lowres_NP_filt(partnum,2)- search_radius_nm/nm_per_pixel;
            
            if XRaw_Right > maxX || XRaw_Left < 0 || YRaw_Down > maxY || YRaw_Up < 0
                continue
            end
            xlimspix=[XRaw_Left XRaw_Right];
            ylimspix=[YRaw_Up YRaw_Down];
            thispart = findSMLs_mapred_ONI(ds_SML, xlimspix, ylimspix, Mypath);
            
            if isempty(thispart)
                continue
            end
            % precision selection for both x and drift correction
            thispart_mat=table2array(thispart);
%             select_thispart_mat=thispart_mat(thispart.XPrecision_nm_ < precision_xy_thr &...
%                                          thispart.YPrecision_nm_ < precision_xy_thr,:);
                                     
            %search again due to some images with bigger drifting
            select_thispart_mid=thispart_mat(thispart.Frame > (mid_frame_NP-100) &...
                                         thispart.Frame < (mid_frame_NP+100),:);
            
            if isempty(select_thispart_mid)
                continue
            end
            % XRaw, YRaw at midpoint
            new_part_loc=[mean(select_thispart_mid(:,13)),mean(select_thispart_mid(:,14))];
            new_xlimspix=new_part_loc(1)+[-1,1]*search_radius_nm/nm_per_pixel;
            new_ylimspix=new_part_loc(2)+[-1,1]*search_radius_nm/nm_per_pixel;
            thispart = findSMLs_mapred_ONI(ds_SML, new_xlimspix, new_ylimspix, Mypath);
            
            clear thispart_mat select_thispart_mat
            % precision selection for both x and drift correction
            thispart_mat=table2array(thispart);
            select_thispart_mat=thispart_mat(thispart.XPrecision_nm_ < precision_xy_thr &...
                                         thispart.YPrecision_nm_ < precision_xy_thr,:);
            
            select_thispart = array2table(select_thispart_mat,...
                'VariableNames',ds_SML.VariableNames);
            
            clear thispart_mat 
            thispartloc=[mean(select_thispart.X_nm_),mean(select_thispart.Y_nm_)];
            
            
            shiftx=new_part_loc(1)-centr_lowres_NP_filt(partnum,1);
            shifty=new_part_loc(2)-centr_lowres_NP_filt(partnum,2);
            


            
            %Save Global data into Particle
            Particle.Fi_SML = Fi_SML;
            Particle.nmperpix = nm_per_pixel;
            Particle.SMLsearchpix = search_radius_nm/nm_per_pixel;
            %Save Local data into Particle
            Particle.partnum = partnum;
            Particle.partlocpix = centr_lowres_NP_filt(partnum,[1,2]);
            Particle.xlimspix = xlimspix;
            Particle.ylimspix = ylimspix;
            Particle.partlocnm = thispartloc;
%             Particle.twocolor = twocol;
            Particle.NumSMLLabel1 = size(select_thispart,1);
            Particle.SMLLabel = select_thispart;

            AllParticles = [AllParticles; Particle];
            
            if size(select_thispart,1)<minimal_SMLs_NP
                continue
            end
            
             if if_outputimage
                             
                %scatter + low res image
                plot_merge=figure;
                pixel_thispart=[select_thispart.XRaw_pix_, select_thispart.YRaw_pix_];
        %       pixel_thispart_corr=[driftcorr_output_pixx(:,par)-mean_pixx_map(par), driftcorr_output_pixy(:,par)-mean_pixy_map(par)];

                ax1 = axes;
                plot_SMLs=scatter(ax1,pixel_thispart(:,1)-shiftx, pixel_thispart(:,2)-shifty, 20, [1,0,0], '.');
                hold on
                colormap(ax1,'jet')
                scale_bar_coord=[xlimspix(1,2)-120/nm_per_pixel ylimspix(1,2)-20/nm_per_pixel;...
                    xlimspix(1,2)-20/nm_per_pixel ylimspix(1,2)-20/nm_per_pixel];
                plot_scalebar=plot(scale_bar_coord(1:2,1),scale_bar_coord(1:2,2),'-c', 'LineWidth', 3);hold on
        %             pixel_clus = cell(size(clus1,1),1);
        %             for k = 1: size(clus1,1)
        %                 pixel_clus = clus1{k,1}.Points*transform_matrix+shift_matrix;
        %                 plot_clus=plot(alphaShape(pixel_clus(:,1),pixel_clus(:,2)),'FaceAlpha',0.1,'FaceColor',[1,0,0],'EdgeAlpha',0);
        %                 hold on
        %             end

                plot_partnum = imshow(Im_NP, [Im_low, Im_high*0.7]);
                axisgfp = gca;
                axisgfp.XLim = xlimspix;
                axisgfp.YLim = ylimspix;
                axisgfp.Visible = 'off';
                hold off
                alpha(0.5)

                OutDir_partnum = strcat(TopDir, matlaboutbase, 'image_output\', SMLs_filename, ...
                    '-part-', num2str(partnum));
                saveas(plot_merge, strcat(OutDir_partnum, '_NP_SML.png'));
                close(plot_merge);
                
              end
           
        end
        
        save(strcat(OutDir, '_NP_Particles.mat'), 'AllParticles');
    end    
end
tEnd = toc(tStart); % Stopwatch stop
disp(strcat('total time spent: ', num2str(floor(tEnd/3600)), ' hours ', ...
num2str(floor(rem(tEnd,3600)/60)), ' minutes and ', num2str(rem(tEnd,60)), ' seconds'))  % Show how long the processing takes



function [tblparse] = findSMLs_mapred_ONI(ds, xlims, ylims, Outputdir)
% extract entries from datastore that lie within x- and ylims
%   Detailed explanation goes here

% parse on 'X Raw pix' and 'Y Raw pix' cols 13 and 14
inbetween = @(data) data{:,13} > xlims(1,1) & data{:,13} < xlims(1,2) &...
                    data{:,14} > ylims(1,1) & data{:,14} < ylims(1,2);

% see
% 'https://www.mathworks.com/help/matlab/import_export/simple-data-subsetting-using-mapreduce.html'
configuredMapper = ...
    @(data, info, intermKVStore) subsettingMapperGeneric(data, info, ...
    intermKVStore, inbetween);

result = mapreduce(ds, configuredMapper, @subsettingReducer,...
    'OutputFolder', strcat(Outputdir,'MapReduceFiles'));

a = readall(result);
tblparse = a.Value{1};

%delete MapReduceFiles\*.mat
end