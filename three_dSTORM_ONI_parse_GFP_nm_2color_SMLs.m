addpath(genpath('M:\Yen-Cheng\YCC_Matlab'));
addpath(genpath('M:\Yi-Han\ONI_Storm\Yi-Han_edited_code\MATLAB_library'))

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
% here since I need to do drift correction to find the SMLs, I did
% 2*search_radius first and then drift correct and then search again within
% 1*search_radius
minimal_SMLs_NP=15000;
%minSMLs_2 = 5;
runDBSCAN = false;
DBSCAN_eps_1 = 15; %20 nm cluster size
%DBSCAN_eps_2 = 20; %30 nm cluster size
if_2D = true;
if_outputimage = true;
clus_max_axis_thr=inf;
clu1_counts=0;

doparse = true;

AF647_channel= 1;
CF568_channel= 0;
GFP_channel= 2;

dye0 = 'CF568';
dye1 = 'AF647';

offset=0;
first_frame_NP=11;
first_frame_NP_after_prebleach=1011;
first_frame_NP_after_prebleach=first_frame_NP_after_prebleach+offset;
first_frame_GFP=1;
last_frame_GFP=10;
total_frame=20000;

[raw_tif_filename, Mypath] = uigetfile({'*.tif; *.tiff';'*.xlsx';'*.mat';'*.*'}, 'Select "_GFP_f1.tif" files', 'MultiSelect','on');
% raw_tif_filename = '510_A_HXB2_2G12-4';
% SMLs_filename = '510_A_HXB2_2G12_SML-5';
% TopDir = '510_A_HXB2\';
TopDir = Mypath;
matlaboutbase = 'output\';
if ischar(raw_tif_filename)==1 % change char array into cell array. This is required when only one file exists in "raw_tif_filename".
    raw_tif_filename = {raw_tif_filename};
end
% driftcorr_matrix_name = '510_A_HXB2_2G12_SML-5_driftcorr';

tfilebase_all={'bead121720_beadSML';...
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
    %for acq=1:10 %:20
        %read low res images
        % first frame: NP, second frame: GFP

%         Im_NP = imread(strcat(TopDir,raw_tif_filename),first_frame_NP);
        temp = split(raw_tif_filename(acq), '_GFP_f1.tif'); 
        SMLs_filename = strcat(char(temp(1)), '_SMLs');
        
        driftcorr_matrix_name = strcat(SMLs_filename, '_driftcorr');

        % load GFP-found 1color SMLs
        load(strcat(TopDir, matlaboutbase, driftcorr_matrix_name,'.mat'));


        Im_GFP = imread(strcat(TopDir,char(raw_tif_filename(acq))));
        %Im_NP = imread(raw_tif_filename,2);
        
        if exist(strcat(TopDir,char(raw_tif_filename(acq))), 'file') == 0
          % File does not exist
          % Skip to bottom of loop and continue with the loop
          continue;
        end
        [maxY,maxX]=size(Im_GFP);

        %% Parse SMLs with GFP presenting; Peak finder to locate centroids in GFP image
        peaksizepix = 5;
        peakintensitythresh =200;
%         [centr_lowres] = one_channel_partID(Im_GFP, peaksizepix, peakintensitythresh, []);
        [centr_lowres_GFP] = one_channel_partID_includeSide(Im_GFP, peaksizepix, peakintensitythresh, []);

        %% GFP peaks
        % centr_lowres is
        % [x,y,meanback,sumsig,corrsumsig,SNR,resolvpow,perc_backpix]
        centr_lowres_GFP_filt = centr_lowres_GFP;
        % filter resolvpow>0.8
        centr_lowres_GFP_filt = centr_lowres_GFP_filt(centr_lowres_GFP_filt(:,7)>0.8,:);
        % filter perc_backpix>50%
        centr_lowres_GFP_filt = centr_lowres_GFP_filt(centr_lowres_GFP_filt(:,8)>50,:);
        % filter SNR>1.5
        centr_lowres_GFP_filt = centr_lowres_GFP_filt(centr_lowres_GFP_filt(:,6)>1.5,:);

%         %% Filter GFP peaks
%         % centr_lowres is
%         % [x,y,meanback,sumsig,corrsumsig,SNR,resolvpow,perc_backpix]
%         centr_lowres_filt = centr_lowres;
%         % filter resolvpow>0.8
%         centr_lowres_filt = centr_lowres_filt(centr_lowres_filt(:,7)>0.8,:);
%         % filter perc_backpix>50%
%         centr_lowres_filt = centr_lowres_filt(centr_lowres_filt(:,8)>50,:);
%         % filter SNR>1.5
%         centr_lowres_filt = centr_lowres_filt(centr_lowres_filt(:,6)>5,:);
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
        Im_low=min(min(Im_GFP));
        Im_high=max(max(Im_GFP))*0.75;
        pfil = plot_spot_filtering(centr_lowres_GFP, centr_lowres_GFP_filt);
        pdet = plot_spot_detection(Im_GFP, centr_lowres_GFP, centr_lowres_GFP_filt, Im_low, Im_high);
        


        %% Check flag to parse ASCII file
        % True - save/close peak filtering/detection images and parse ASCII
        % False - end script
%         if ~exist(strcat(matlaboutbase, SMLs_filename), 'dir');
%             mkdir(strcat(matlaboutbase, SMLs_filename));
%             addpath(fullfile(pwd, strcat(matlaboutbase, SMLs_filename)));
%         end
        %OutDir= strcat(TopDir, matlaboutbase, SMLs_filename,'_');
        OutDir= strcat(TopDir, matlaboutbase, SMLs_filename);
        if doparse
            saveas(pfil, strcat(OutDir, '_GFP_part-filtering.png'));
            close(pfil);
            saveas(pdet, strcat(OutDir, '_GFP_part-detection.png'));
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
        
        AllGFPParticles = [];
        AllParticles_0 = [];       
        AllParticles_1 = [];       
        for partnum = 1:size(centr_lowres_GFP_filt,1)
            size(centr_lowres_GFP_filt,1);
            XRaw_Right=centr_lowres_GFP_filt(partnum,1)+ 2*search_radius_nm/nm_per_pixel;
            XRaw_Left=centr_lowres_GFP_filt(partnum,1)- 2*search_radius_nm/nm_per_pixel;
            YRaw_Down=centr_lowres_GFP_filt(partnum,2)+ 2*search_radius_nm/nm_per_pixel;
            YRaw_Up=centr_lowres_GFP_filt(partnum,2)- 2*search_radius_nm/nm_per_pixel;
            
            % if XRaw_Right > maxX || XRaw_Left < 0 || YRaw_Down > maxY || YRaw_Up < 0
            % Select particles located in the middle region only
            if XRaw_Right > maxX || XRaw_Left < 0 || YRaw_Down > maxY || YRaw_Up < 0

                continue
            end
            xlimspix=[XRaw_Left XRaw_Right];
            ylimspix=[YRaw_Up YRaw_Down];
            %low res GFP to GFP SMLs
            thispartRaw = findSMLs_mapred_ONI(ds_SML, xlimspix, ylimspix, Mypath);
            
            if isempty(thispartRaw)
                continue
            end
            
            thispartRaw_mat=table2array(thispartRaw);
            select_thispartRaw_mat=thispartRaw_mat(thispartRaw.XPrecision_nm_ < precision_xy_thr &...
                                         thispartRaw.YPrecision_nm_ < precision_xy_thr &...
                                         thispartRaw.Channel == GFP_channel,:);
            if isempty(select_thispartRaw_mat)
                continue
            end
            thispartRaw = array2table(select_thispartRaw_mat,...
                                 'VariableNames',ds_SML.VariableNames);
            clear select_thispartRaw_mat
            
            %Save Global GFP SMLs data into Particle
            GFPParticles.Fi_SML = Fi_SML;
            GFPParticles.nmperpix = nm_per_pixel;
            GFPParticles.SMLsearchpix = search_radius_nm/nm_per_pixel;
            %Save Local GFP SMLs data into Particle
            GFPParticles.partnum = partnum;
            GFPParticles.partlocpix = centr_lowres_GFP_filt(partnum,[1,2]);
            GFPParticles.xlimspix = xlimspix;
            GFPParticles.ylimspix = ylimspix;
%             Particle.partlocnm = thispartloc;
%             Particle.twocolor = twocol;
            GFPParticles.NumGFPSML = size(thispartRaw,1);
            GFPParticles.GFPSML = thispartRaw;

            AllGFPParticles = [AllGFPParticles; GFPParticles];

            GFP_loc_x_nm=thispartRaw.X_nm_;
            GFP_loc_y_nm=thispartRaw.Y_nm_;
            
            GFP_loc_x_pix=thispartRaw.X_nm_/nm_per_pixel;
            GFP_loc_y_pix=thispartRaw.Y_nm_/nm_per_pixel;
            
            shiftx=GFP_loc_x_nm(first_frame_GFP)/nm_per_pixel-centr_lowres_GFP_filt(partnum,1);
            shifty=GFP_loc_y_nm(first_frame_GFP)/nm_per_pixel-centr_lowres_GFP_filt(partnum,2);
            
            % search before drift correction, 2*search_radius_nm
            GFP_loc_xlimsnm_beforedrift = GFP_loc_x_nm(first_frame_GFP)+ [-2,2]*search_radius_nm;
            GFP_loc_ylimsnm_beforedrift = GFP_loc_y_nm(first_frame_GFP)+ [-2,2]*search_radius_nm;
            
            if GFP_loc_xlimsnm_beforedrift(2) > maxX*search_radius_nm || GFP_loc_xlimsnm_beforedrift(1) < 0 ...
                || GFP_loc_ylimsnm_beforedrift(2) > maxY*search_radius_nm || GFP_loc_ylimsnm_beforedrift(1) < 0
                continue
            end
            % SMLs from Channel 0 
            thispart = findSMLs_mapred_ONI_nm_colormapping(ds_SML, GFP_loc_xlimsnm_beforedrift, GFP_loc_ylimsnm_beforedrift, Mypath);
            if isempty(thispart)
                continue
            end
            % precision selection for both x and drift correction
            thispart_0_mat=table2array(thispart);
            select_thispart_0_mat=thispart_0_mat(thispart.XPrecision_nm_ < precision_xy_thr &...
                                         thispart.YPrecision_nm_ < precision_xy_thr &...
                                         thispart.Channel == CF568_channel &...
                                         thispart.Frame >= (first_frame_NP_after_prebleach-1) &...
                                         thispart.Frame < (first_frame_NP_after_prebleach+total_frame-1),:);
            if isempty(select_thispart_0_mat)
                continue
            end
            %drift correction
            for i=1:size(select_thispart_0_mat,1)
                % exclude signals from other channels

                tframe=select_thispart_0_mat(i,2); %frame format from ONI starts from 0
                tframe=tframe+1-first_frame_NP_after_prebleach+1;
                %drift correction for X_nm_(3), Y_nm_(4), X_pix_(8), Y_pix_(9), XRaw_pix_(13),
                %YRaw_pix_(14)
                select_thispart_0_mat(i,3)=select_thispart_0_mat(i,3)+ driftcorr_matrix(tframe,1)*nm_per_pixel;
                select_thispart_0_mat(i,4)=select_thispart_0_mat(i,4)+ driftcorr_matrix(tframe,2)*nm_per_pixel;
                select_thispart_0_mat(i,8)=select_thispart_0_mat(i,8)+ driftcorr_matrix(tframe,1);
                select_thispart_0_mat(i,9)=select_thispart_0_mat(i,9)+ driftcorr_matrix(tframe,2);
                select_thispart_0_mat(i,13)=select_thispart_0_mat(i,13)+ driftcorr_matrix(tframe,1);
                select_thispart_0_mat(i,14)=select_thispart_0_mat(i,14)+ driftcorr_matrix(tframe,2);
            end
            select_thispart_0 = array2table(select_thispart_0_mat,...
                'VariableNames',ds_SML.VariableNames);
            clear thispart_mat
            
            select_thispart_0_mat=table2array(select_thispart_0);
            
            % search after drift correction, 1*search_radius_nm
            GFP_loc_xlimsnm_afterdrift=GFP_loc_x_nm(first_frame_GFP)+ [-1,1]*search_radius_nm;
            GFP_loc_ylimsnm_afterdrift=GFP_loc_y_nm(first_frame_GFP)+ [-1,1]*search_radius_nm;            
            
            GFP_loc_xlimspix_afterdrift=GFP_loc_x_pix(first_frame_GFP)+ [-1,1]*search_radius_nm/nm_per_pixel;
            GFP_loc_ylimspix_afterdrift=GFP_loc_y_pix(first_frame_GFP)+ [-1,1]*search_radius_nm/nm_per_pixel;
            
            select_thispart_0_mat=select_thispart_0_mat(select_thispart_0.X_nm_ > GFP_loc_xlimsnm_afterdrift(1) &...
                                         select_thispart_0.X_nm_ < GFP_loc_xlimsnm_afterdrift(2) &...
                                         select_thispart_0.Y_nm_ > GFP_loc_ylimsnm_afterdrift(1) &...
                                         select_thispart_0.Y_nm_ < GFP_loc_ylimsnm_afterdrift(2),:);
            if isempty(select_thispart_0_mat)
                continue
            end
            % clear select_thispart
            % select_thispart = array2table(select_thispart_0_mat,...
            %     'VariableNames',ds_SML.VariableNames);
            
            clear select_thispart_0_mat 
            thispartloc=[mean(select_thispart_0.X_nm_),mean(select_thispart_0.Y_nm_)];
            
            %drift correction
            


            
            %Save Global AF647 SMLs data into Particle
            Particle_0.Fi_SML = Fi_SML;
            Particle_0.nmperpix = nm_per_pixel;
            Particle_0.SMLsearchpix = search_radius_nm/nm_per_pixel;
            %Save Local data into Particle
            Particle_0.partnum = partnum;
            Particle_0.partlocpix = [GFP_loc_x_pix(1),GFP_loc_y_pix(1)];
            Particle_0.part_GFPintensity = centr_lowres_GFP_filt(partnum,5);
            Particle_0.xlimspix = GFP_loc_xlimspix_afterdrift;
            Particle_0.ylimspix = GFP_loc_ylimspix_afterdrift;
            Particle_0.xlimsnm = GFP_loc_xlimsnm_afterdrift;
            Particle_0.ylimsnm = GFP_loc_ylimsnm_afterdrift;
            Particle_0.partlocnm = [thispartRaw.X_nm_(1),thispartRaw.Y_nm_(1)];
%             Particle.twocolor = twocol;
            Particle_0.NumSMLLabel1 = size(select_thispart_0,1);
            Particle_0.SMLLabel = select_thispart_0;

            AllParticles_0 = [AllParticles_0; Particle_0];
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            thispart_1_mat=table2array(thispart);
            select_thispart_1_mat=thispart_1_mat(thispart.XPrecision_nm_ < precision_xy_thr &...
                                         thispart.YPrecision_nm_ < precision_xy_thr &...
                                         thispart.Channel == AF647_channel &...
                                         thispart.Frame >= (first_frame_NP_after_prebleach-1) &...
                                         thispart.Frame < (first_frame_NP_after_prebleach+total_frame-1),:);
            if isempty(select_thispart_1_mat)
                continue
            end
            %drift correction
            for i=1:size(select_thispart_1_mat,1)
                % exclude signals from other channels

                tframe=select_thispart_1_mat(i,2); %frame format from ONI starts from 0
                tframe=tframe+1-first_frame_NP_after_prebleach+1;
                %drift correction for X_nm_(3), Y_nm_(4), X_pix_(8), Y_pix_(9), XRaw_pix_(13),
                %YRaw_pix_(14)
                select_thispart_1_mat(i,3)=select_thispart_1_mat(i,3)+ driftcorr_matrix(tframe,1)*nm_per_pixel;
                select_thispart_1_mat(i,4)=select_thispart_1_mat(i,4)+ driftcorr_matrix(tframe,2)*nm_per_pixel;
                select_thispart_1_mat(i,8)=select_thispart_1_mat(i,8)+ driftcorr_matrix(tframe,1);
                select_thispart_1_mat(i,9)=select_thispart_1_mat(i,9)+ driftcorr_matrix(tframe,2);
                select_thispart_1_mat(i,13)=select_thispart_1_mat(i,13)+ driftcorr_matrix(tframe,1);
                select_thispart_1_mat(i,14)=select_thispart_1_mat(i,14)+ driftcorr_matrix(tframe,2);
            end
            select_thispart_1 = array2table(select_thispart_1_mat,...
                'VariableNames',ds_SML.VariableNames);
            clear thispart_mat
            
            select_thispart_1_mat=table2array(select_thispart_1);
            
            % search after drift correction, 1*search_radius_nm
            GFP_loc_xlimsnm_afterdrift=GFP_loc_x_nm(first_frame_GFP)+ [-1,1]*search_radius_nm;
            GFP_loc_ylimsnm_afterdrift=GFP_loc_y_nm(first_frame_GFP)+ [-1,1]*search_radius_nm;            
            
            GFP_loc_xlimspix_afterdrift=GFP_loc_x_pix(first_frame_GFP)+ [-1,1]*search_radius_nm/nm_per_pixel;
            GFP_loc_ylimspix_afterdrift=GFP_loc_y_pix(first_frame_GFP)+ [-1,1]*search_radius_nm/nm_per_pixel;
            
            select_thispart_1_mat=select_thispart_1_mat(select_thispart_1.X_nm_ > GFP_loc_xlimsnm_afterdrift(1) &...
                                         select_thispart_1.X_nm_ < GFP_loc_xlimsnm_afterdrift(2) &...
                                         select_thispart_1.Y_nm_ > GFP_loc_ylimsnm_afterdrift(1) &...
                                         select_thispart_1.Y_nm_ < GFP_loc_ylimsnm_afterdrift(2),:);
            if isempty(select_thispart_1_mat)
                continue
            end
            % clear select_thispart
            % select_thispart = array2table(select_thispart_1_mat,...
            %     'VariableNames',ds_SML.VariableNames);
            
            clear select_thispart_1_mat 
            thispartloc=[mean(select_thispart_1.X_nm_),mean(select_thispart_1.Y_nm_)];
            
            %drift correction
            


            
            %Save Global AF647 SMLs data into Particle
            Particle_1.Fi_SML = Fi_SML;
            Particle_1.nmperpix = nm_per_pixel;
            Particle_1.SMLsearchpix = search_radius_nm/nm_per_pixel;
            %Save Local data into Particle
            Particle_1.partnum = partnum;
            Particle_1.partlocpix = [GFP_loc_x_pix(1),GFP_loc_y_pix(1)];
            Particle_1.part_GFPintensity = centr_lowres_GFP_filt(partnum,5);
            Particle_1.xlimspix = GFP_loc_xlimspix_afterdrift;
            Particle_1.ylimspix = GFP_loc_ylimspix_afterdrift;
            Particle_1.xlimsnm = GFP_loc_xlimsnm_afterdrift;
            Particle_1.ylimsnm = GFP_loc_ylimsnm_afterdrift;
            Particle_1.partlocnm = [thispartRaw.X_nm_(1),thispartRaw.Y_nm_(1)];
%             Particle.twocolor = twocol;
            Particle_1.NumSMLLabel1 = size(select_thispart_1,1);
            Particle_1.SMLLabel = select_thispart_1;

            AllParticles_1 = [AllParticles_1; Particle_1];
            
%             if size(select_thispart,1)<minimal_SMLs_GFP
%                 continue
%             end
            
             if if_outputimage 
                
                xlimspix(1)=GFP_loc_xlimspix_afterdrift(1)-shiftx-0.5;
                xlimspix(2)=GFP_loc_xlimspix_afterdrift(2)-shiftx+0.5;
                ylimspix(1)=GFP_loc_ylimspix_afterdrift(1)-shifty-0.5;
                ylimspix(2)=GFP_loc_ylimspix_afterdrift(2)-shifty+0.5;
                %scatter + low res image
                plot_merge=figure;
                thispart_nm_0=[select_thispart_0.X_nm_, select_thispart_0.Y_nm_];
                pixel_thispart_0=thispart_nm_0/nm_per_pixel;
        %       pixel_thispart_corr=[driftcorr_output_pixx(:,par)-mean_pixx_map(par), driftcorr_output_pixy(:,par)-mean_pixy_map(par)];

                thispart_nm_1=[select_thispart_1.X_nm_, select_thispart_1.Y_nm_];
                pixel_thispart_1=thispart_nm_1/nm_per_pixel;

                ax1 = axes;

                
                colormap(ax1,'turbo')
                scale_bar_coord=[xlimspix(1,2)-120/nm_per_pixel ylimspix(1,2)-20/nm_per_pixel;...
                    xlimspix(1,2)-20/nm_per_pixel ylimspix(1,2)-20/nm_per_pixel];
                plot_scalebar=plot(scale_bar_coord(1:2,1),scale_bar_coord(1:2,2),'-c', 'LineWidth', 3);hold on
        %             pixel_clus = cell(size(clus1,1),1);
        %             for k = 1: size(clus1,1)
        %                 pixel_clus = clus1{k,1}.Points*transform_matrix+shift_matrix;
        %                 plot_clus=plot(alphaShape(pixel_clus(:,1),pixel_clus(:,2)),'FaceAlpha',0.1,'FaceColor',[1,0,0],'EdgeAlpha',0);
        %                 hold on
        %             end
                % legend(dye0,dye1);
                plot_partnum = imshow(Im_GFP, [Im_low, Im_high*0.7]);
                axisgfp = gca;
                axisgfp.XLim = xlimspix(1,:);
                axisgfp.YLim = ylimspix(1,:);
                axisgfp.Visible = 'on'; % show coordinate
                alpha(0.4)
                hold on
                plot_SMLs0=scatter(pixel_thispart_0(:,1)-shiftx, pixel_thispart_0(:,2)-shifty , 20, [1,0.61,0], 'filled', 'MarkerFaceAlpha',0.7);
                hold on
                scatter(pixel_thispart_1(:,1)-shiftx, pixel_thispart_1(:,2)-shifty, 20, [1,0,0], 'filled', 'MarkerFaceAlpha',0.7);
                
                hold off
                

                OutDir_partnum = strcat(TopDir, matlaboutbase, 'image_output\GFP_2color\');
                if exist(OutDir_partnum, 'dir') == 0
                    mkdir(OutDir_partnum);
                end
                addpath(genpath(OutDir_partnum))
                
                saveas(plot_merge, strcat(OutDir_partnum, SMLs_filename, '-part-', num2str(partnum), '_GFP_2colorSMLs.png'));
                close(plot_merge);
              end
           
        end
        save(strcat(OutDir, '_GFP_Particles_GFPSMLs.mat'), 'AllGFPParticles'); % Information of GFP particles only.
        save(strcat(OutDir, '_GFP_Particles_0.mat'), 'AllParticles_0'); % Information of SMLs on each GFP particle.
        save(strcat(OutDir, '_GFP_Particles_1.mat'), 'AllParticles_1'); % Information of SMLs on each GFP particle.

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

function [tblparse] = findSMLs_mapred_ONI_nm_colormapping(ds, xlims, ylims, Outputdir)
% extract entries from datastore that lie within x- and ylims
%   Detailed explanation goes here

% parse on 'X nm' and 'Y nm' cols 3 and 4
% 'X nm' and 'Y nm' are nm loc from 2-color mapping of ONI
inbetween = @(data) data{:,3} > xlims(1,1) & data{:,3} < xlims(1,2) &...
                    data{:,4} > ylims(1,1) & data{:,4} < ylims(1,2);

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
