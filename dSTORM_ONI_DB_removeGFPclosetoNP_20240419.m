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
photon_thr=1000;
continuous_ratio_threshold=1;
thr_too_close_pix=5;

search_radius_nm = 200; %particle size = 200 nm
maxSMLs_1=3000;
minSMLs_1=[20];% 90 120 180];
% minSMLs_2 = 5;
runDBSCAN = false;
DBSCAN_eps_1 = 15; %15 nm cluster size
%DBSCAN_eps_2 = 20; %20 nm cluster size
if_2D = true;
if_pdist = true;
if_outputimage = false; %plot and save the images
clus_max_axis_thr=inf;
clu1_counts=0;

doparse = true;
first_frame_NP=11;
first_frame_GFP=1;
total_frame=20000;

first_frame_NP_after_prebleach=1011;


% Automatically plot the clustering ratio for all the input data. 
% Yes = 1; No = 0.
plot_clus_r = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Yi-Han %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Sample = {'A', 'B', 'C', 'D', 'E', 'F'}; % Manually type in the sample name.

% Change the value of following two variables based on previous analysis.
% If there's no previous analysis, leave these two value as 1
stand_SML_n = [1 2 3 4 5 6];
Curr_SML_n = [1 2 3 4 5 6];

SMLs_result(2:length(Sample)+1,1) = Sample;
for s = 1:length(Sample)
    if s == 1
        [MyGFPfiles,Mypath] = uigetfile({'*.tif; *.tiff';'*.xlsx';'*.mat';'*.*'}, ...
        strcat('Select "GFP_f1.tif" files for sample', " ", Sample(s)), 'MultiSelect','on');
    else
        [MyGFPfiles,Mypath] = uigetfile({'*.tif; *.tiff';'*.xlsx';'*.mat';'*.*'}, ...
        strcat('Select "GFP_f1.tif" files for sample', " ", Sample(s)), 'MultiSelect','on', Mypath);
    end

    AllPath.(strcat('Mypath_', Sample{s})) = Mypath;
    AllGFPfiles.(strcat('MyGFPfiles_', Sample{s})) = MyGFPfiles;

end
% Use seperate loops to avoid user waiting to select the files.
name1 = fieldnames(AllPath);
name2 = fieldnames(AllGFPfiles);
for s = 1:length(Sample)
    
    Mypath = AllPath.(name1{s});
    MyGFPfiles = AllGFPfiles.(name2{s});


    if ischar(MyGFPfiles)==1 % change char array into cell array. This is required when only one file exists in "MyGFPfiles".
        MyGFPfiles = {MyGFPfiles};
    end
    
 
    
    
    MyGFPdata=cell(2,length(MyGFPfiles));
    img_out = strcat(Mypath,'analysis\', 'output\', 'image_output\');
    clusinfo = strcat(Mypath,'analysis\', 'output\', 'clusinfo\');
    mkdir(img_out)
    mkdir(clusinfo)
    addpath(genpath(Mypath)); % add the specified folder and subfolders to the path.


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Yi-Han %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    SMLs_filename_mat = '_GFP_Particles';
    GFPSMLs_filename_mat = '_GFP_Particles_GFPSMLs';
    NPSMLs_filename_mat = '_NP_Particles';

    
    matlaboutbase = strcat(Mypath, 'output042024\');
    
    % Automatically analyze data for all minimum SMLs
    for i = 1:length(minSMLs_1)
        CurrminSML = minSMLs_1(i);
        cumdist_AF647=[];
        all_photons=[];
        compileSMLs_clus=[];
        all_photon_trace=[];
        all_continuous_ratio=[];
        all_dist_btw_f=[];
        thr_dist_btw_f=[];
        all_par_remove=[];

        for acq=1:length(MyGFPfiles) 
            temp = split(MyGFPfiles{acq}, '.');
            MyGFPdata{1,acq} = temp(1);    
            MyGFPdata{2,acq} = imread(strcat(Mypath,MyGFPfiles{acq}));
            temp1 = split(temp(1), '_GFP_f1');  % filename with number
            temp2 = split(char(temp1(1)), '-');
            if isequal(temp2{1} ,'Yi') &  length(temp2) ==2 % function specifically needs for filename starts with Yi-Han 
                sample_name = strjoin(temp2,'-');
            else
                sample_name = strjoin(temp2(1:(end-1)),'-'); % % filename without number
            end
    
            if contains(sample_name,'AB_') % if the filename contain 'AB_'
                sample_name = split(sample_name, 'AB_');
                sample_name = string(sample_name(2));
            end
            base_name = strcat(char(temp1{1}), '_SMLs');
    
            
            par_to_remove=[];
            %read low res images
            % first frame: NP, second frame: GFP
            if exist(MyGFPfiles{1,acq}, 'file') == 0
              % File does not exist
              % Skip to bottom of loop and continue with the loop
              continue;
            end
    %         if exist(strcat(TopDir,num2str(acq)), 'file') == 0
    %           % File does not exist
    %           % Skip to bottom of loop and continue with the loop
    %           continue;
    %         end   
    %       Im_NP = imread(strcat(TopDir,raw_tif_filename),first_frame_NP);
            Im_GFP = imread(MyGFPfiles{1,acq});
    %        imwrite(Im_GFP,strcat(TopDir,raw_tif_filename,'GFP_f1.tiff'))
    %        Im_GFP = imread(strcat(TopDir,raw_tif_filename,'GFP_f1.tiff'),1);
            %Im_NP = imread(raw_tif_filename,2);
    
    
    
            try 
                % load NP SMLs
                load(strcat(matlaboutbase,...
                    base_name,NPSMLs_filename_mat,'.mat'));
                AllNPParticles=AllParticles;
                % load GFP SMLs
                load(strcat(matlaboutbase,...
                    base_name,GFPSMLs_filename_mat,'.mat'));
                [maxY,maxX]=size(Im_GFP);
                % load AF647 SMLs
                load(strcat(matlaboutbase,...
                    base_name,SMLs_filename_mat,'.mat'));
            catch ME
                try
                    fil_num = split(base_name, '_SMLs');
                    fil_num = split(fil_num(1), '-');
                    % load NP SMLs
                    load(strcat(matlaboutbase,...
                        sample_name, '_SMLs-', char(fil_num(2)),NPSMLs_filename_mat,'.mat'));
                    AllNPParticles=AllParticles;
                    % load GFP SMLs
                    load(strcat(matlaboutbase,...
                        sample_name, '_SMLs-', char(fil_num(2)),GFPSMLs_filename_mat,'.mat'));
                    [maxY,maxX]=size(Im_GFP);
                    % load AF647 SMLs
                    load(strcat(matlaboutbase,...
                        sample_name, '_SMLs-', char(fil_num(2)),SMLs_filename_mat,'.mat'));
                catch
                    % load NP SMLs
                    NPSMLs_name = strcat('Select "NP_Particles" files for', " ", base_name);
                    [NPSMLs_file,NPSMLs_path] = uigetfile({'*.mat';'*.*'}, ...
                        NPSMLs_name, 'MultiSelect','off', Mypath);
                    load(fullfile(NPSMLs_path, NPSMLs_file));
                    AllNPParticles=AllParticles;
                    % load GFP SMLs
                    GFPSMLs_name = strcat('Select "GFP_Particles_GFPSMLs" files for', " ", base_name);
                    [GFPSMLs_file,GFPSMLs_path] = uigetfile({'*.mat';'*.*'}, ...
                        GFPSMLs_name, 'MultiSelect','off',NPSMLs_path);
                    load(fullfile(GFPSMLs_path, GFPSMLs_file));
                    [maxY,maxX]=size(Im_GFP);
                    % load AF647 SMLs
                    SMLs_name = strcat('Select "GFP_Particles" files for', " ", base_name);
                    [SMLs_file,SMLs_path] = uigetfile({'*.mat';'*.*'}, ...
                        SMLs_name, 'MultiSelect','off',NPSMLs_path);
                    load(fullfile(SMLs_path, SMLs_file));
                                    
                end
    
    
            end
    
            GFP_pix=[];
    
            for GFPnum=1:size(AllParticles,1)
                GFP_pix = cat(1,GFP_pix,AllParticles(GFPnum).partlocpix );    
            end
            dist_mat=squareform(pdist(GFP_pix));
            n = size(dist_mat,1);
            dist_mat(1:(n+1):end) = nan;
            [r_too_close,c_too_close]=find(dist_mat<thr_too_close_pix);
    
            if ~isempty(r_too_close)
                par_to_remove=cat(1,par_to_remove,r_too_close);
            end
    
            GFP_round_pix = round(GFP_pix);
    
            %remove GFPs when it's too close to NPs, using X_pix_ to match with
            %GFP coord
    
    
            NP_pix=[];
    
            for NPnum=size(AllNPParticles,1)
    
                t_r=find(AllNPParticles(NPnum).SMLLabel.Frame==(first_frame_NP-1));
                if isempty(t_r)
                    continue
                end
                NP_pix = cat(1,NP_pix,[AllNPParticles(NPnum).SMLLabel.X_pix_(t_r) AllNPParticles(NPnum).SMLLabel.Y_pix_(t_r)]);
    
                tNP_round_pix = round([AllNPParticles(NPnum).SMLLabel.X_pix_(t_r) AllNPParticles(NPnum).SMLLabel.Y_pix_(t_r)]);
    
                tsub_GFP_NP = GFP_round_pix-tNP_round_pix;
                tdist_GFP_NP = sqrt(sum(tsub_GFP_NP.*tsub_GFP_NP,2));
    
                %subtraction of GFP if GFP is too close to NP, threshold: thr_too_close_pix
                t_close_GFPnum = find(tdist_GFP_NP<thr_too_close_pix);
                if isempty(t_close_GFPnum)
                    continue
                else
                   for GFPnum=1:size(t_close_GFPnum,1) 
    
                       if AllNPParticles(NPnum).NumSMLLabel1 > maxSMLs_1
                          par_to_remove=cat(1,par_to_remove,t_close_GFPnum);
                       end
                   end
                end
    
            end
            par_to_remove=unique(par_to_remove);
            AllParticles(par_to_remove,:)=[];
    
            all_par_remove=cat(1,all_par_remove,par_to_remove);
    
            %% Create peak detection and filtering plots
            Im_low=min(min(Im_GFP));
            Im_high=max(max(Im_GFP))*0.75;
    
    
            %%
            %read SMLs from .mat file
    
    
            for par = 1:size(AllParticles,1)
    
    
                SMLs_to_lowres_corr=AllParticles(par).partlocpix-AllGFPParticles(par).partlocpix;
                thispart_info=zeros(1,6); 
                %1: num of SMLs, 2: clus num, 3:mean clus density, 4:mean clus
                %area, 5: average SMLs per clus, 6: GFP intensity
                thispart=AllParticles(par).SMLLabel;
                thispart_mat=table2array(thispart);
    
                select_thispart_mat=thispart_mat(AllParticles(par).SMLLabel.Frame<(first_frame_NP_after_prebleach+total_frame-1),:);
                %threshold: photons, column 16
                select_thispart_mat=select_thispart_mat(AllParticles(par).SMLLabel.Photons>photon_thr,:);
                
                % Randomly select the SMLs
                select_thispart_mat = rand_select(select_thispart_mat, stand_SML_n(s), Curr_SML_n(s)); 
    
                AllParticles(par).NumSMLLabel1=size(select_thispart_mat,1);
                AllParticles(par).SMLLabel= array2table(select_thispart_mat,...
                                            'VariableNames',thispart.Properties.VariableNames);
    
                if AllParticles(par).NumSMLLabel1>maxSMLs_1 || ...
                        AllParticles(par).NumSMLLabel1==0
                    continue
                end
    
                % GFP shift 
                % GFP moving
                % GFP intensity
                %%%%%%% SML continuity filter %%%%%%%% 
    
                thispart_frame=AllParticles(par).SMLLabel.Frame; % total frame number of SMLs appear on this virion
                tdiff=diff(thispart_frame);
    
                xy1_pix_x=AllParticles(par).SMLLabel.X_pix_-SMLs_to_lowres_corr(1);
                xy1_pix_y=AllParticles(par).SMLLabel.Y_pix_-SMLs_to_lowres_corr(2);
                xy1_pix=[xy1_pix_x,xy1_pix_y];
                xy1_nm=xy1_pix*nm_per_pixel;
    
                t_continuous_count=0;
                par_dist_when_contin=[];
    
                for tt=1:(size(thispart_frame)-1)
                    if tdiff(tt)==1
                       tdistance=pdist([xy1_nm(tt);xy1_nm(tt+1)]);
                       par_dist_when_contin=cat(1,par_dist_when_contin,tdistance);
                       if tdistance<precision_xy_thr
                           t_continuous_count=t_continuous_count+1;
                       end
                    end
    
                end
                all_dist_btw_f=cat(1,all_dist_btw_f,par_dist_when_contin);
    
                t_continuous_ratio=t_continuous_count/size(thispart_frame,1);
    
                if t_continuous_ratio>continuous_ratio_threshold
                    continue
                end
                all_continuous_ratio=cat(1,all_continuous_ratio,t_continuous_ratio);
                thr_dist_btw_f=cat(1,thr_dist_btw_f,par_dist_when_contin);
                %%%%%%% SML continuity filter %%%%%%%% 
    
    
                %%%%%%% SML photon count filter %%%%%%%%
                thispart_photons=AllParticles(par).SMLLabel.Photons;
                tphoton_trace = zeros(total_frame,1);
                for tframe = 1: total_frame
                   if any(thispart_frame==(first_frame_NP_after_prebleach+tframe-1))
                      ind = find(thispart_frame==(first_frame_NP_after_prebleach+tframe-1));
                      if size(ind,1)>1
                          ind=ind(1);
                      end
                      tphoton_trace(tframe,1) = thispart_photons(ind);
                   end   
                end
                all_photon_trace=cat(2,all_photon_trace,tphoton_trace);
                %%%%%%% SML photon count filter %%%%%%%%
    
                thispart_info(1)=AllParticles(par).NumSMLLabel1;




                pdist_per_par = pdist(xy1_nm)';
                cumdist_AF647 = cat(1, cumdist_AF647, pdist_per_par);
                thispart_info(6)=AllParticles(par).part_GFPintensity;
                tphotons=AllParticles(par).SMLLabel.Photons;
                all_photons=cat(1,all_photons,tphotons);
                if AllParticles(par).NumSMLLabel1<= CurrminSML
                    thispart_info(2)=0;
                    thispart_info(3)=0;
                    thispart_info(4)=0;
                    thispart_info(5)=0;
                end

                [par1, clus1] = clusteranalysis_2D_perclus_den_area_peri(xy1_nm, DBSCAN_eps_1, CurrminSML);

                thispart_info(2)=par1(1,2);
                thispart_info(3)=par1(1,4); %density
                thispart_info(4)=par1(1,5); %area
                thispart_info(5)=sum(cell2mat(clus1(:,4)))/par1(1,2); %ave SMLs

                compileSMLs_clus = cat(1, compileSMLs_clus, thispart_info);

                xlimspix(1)=AllParticles(par).xlimspix(1)-SMLs_to_lowres_corr(1)-0.5;
                xlimspix(2)=AllParticles(par).xlimspix(2)-SMLs_to_lowres_corr(1)+0.5;
                ylimspix(1)=AllParticles(par).ylimspix(1)-SMLs_to_lowres_corr(2)-0.5;
                ylimspix(2)=AllParticles(par).ylimspix(2)-SMLs_to_lowres_corr(2)+0.5;

                if if_outputimage 


                    %scatter + low res image
                    plot_merge=figure;
                    pixel_thispart=xy1_pix;
            %       pixel_thispart_corr=[driftcorr_output_pixx(:,par)-mean_pixx_map(par), driftcorr_output_pixy(:,par)-mean_pixy_map(par)];

                    ax1 = axes;
                    plot_SMLs=scatter(ax1,pixel_thispart(:,1), pixel_thispart(:,2), 20, [1,0,0], '.');
                    hold on
                    colormap(ax1,'jet')
                    scale_bar_coord=[xlimspix(1,2)-120/nm_per_pixel ylimspix(1,2)-20/nm_per_pixel;...
                        xlimspix(1,2)-20/nm_per_pixel ylimspix(1,2)-20/nm_per_pixel];
                    plot_scalebar=plot(scale_bar_coord(1:2,1),scale_bar_coord(1:2,2),'-c', 'LineWidth', 3);hold on
                    pixel_clus = cell(size(clus1,1),1);
                    for k = 1: size(clus1,1)
                        pixel_clus = clus1{k,1}.Points/nm_per_pixel;
                        plot_clus=plot(alphaShape(pixel_clus(:,1),pixel_clus(:,2)),'FaceAlpha',0.1,'FaceColor',[1,0,0],'EdgeAlpha',0);
                        hold on
                    end

                    plot_partnum = imshow(Im_GFP, [Im_low, Im_high*0.7]);
                    axisgfp = gca;
                    axisgfp.XLim = xlimspix(1,:);
                    axisgfp.YLim = ylimspix(1,:);
                    axisgfp.Visible = 'off';
                    hold off
                    alpha(0.5)

                    OutDir_partnum = strcat(img_out, base_name, ...
                        '-part-', num2str(par),' DB_',num2str(DBSCAN_eps_1),'nm_',num2str(CurrminSML),'SMLs');
                    saveas(plot_merge, strcat(OutDir_partnum, '_GFP_1colorSMLs.png'));
                    close(plot_merge);
                end
            end
            
        end


        %OutDirSum = strcat(matlaboutbase, 'clusinfo\', base_name,' DB',num2str(DBSCAN_eps_1),'nm',num2str(minSMLs_1),'SMLs');
        OutDirSum = strcat(clusinfo, sample_name,'_DB',num2str(DBSCAN_eps_1),'nm_',num2str(CurrminSML),'SMLs');
        save(strcat(OutDirSum, '_compile_SMLsandClus.mat'),'compileSMLs_clus');
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Yi-Han %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%% Save excel file with different SMLs numbers in different sheets %%%%%%%% 
        NewcompileSML = GFPfilter(compileSMLs_clus, 100000);    % delete clusters with intensities above the th
        title = {("SMLs") ("clus #") ("mean density") ("mean area") ("SMLs per clus") ("Vpr-GFP")};
        NewcompileSML2 = [title; num2cell(NewcompileSML)]; 
        %%% set the shared parent folder of two samples %%%
        parent_path = find_parent_dir(AllPath);
        
        writecell(NewcompileSML2,strcat(parent_path, sample_name, '_clus_result', '.xlsx'), ...
            'sheet',strcat(sample_name, num2str(CurrminSML),'SMLs'))
        SMLs_result{1, i+1} = minSMLs_1(i);
        SMLs_result{s+1, i+1} = NewcompileSML;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Yi-Han %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        if if_pdist
           save(strcat(OutDirSum, '_pdist_af647_af647.mat'),'cumdist_AF647', '-v7.3');
        end
        save(strcat(OutDirSum, '_allphotons.mat'),'all_photons');
        save(strcat(OutDirSum, '_allphotontrace.mat'),'all_photon_trace');
        save(strcat(OutDirSum, '_all_contin_ratio.mat'),'all_continuous_ratio');
        save(strcat(OutDirSum, '_all_dist_btw_frames.mat'),'all_dist_btw_f');
        save(strcat(OutDirSum, '_thr_dist_btw_frames.mat'),'thr_dist_btw_f');
        save(strcat(OutDirSum, '_row_too_close_to_remove.mat'),'all_par_remove');
        
    end
    %%% display median of SML per Virion
    M = median(NewcompileSML(:, 1),"omitmissing");
    disp(sample_name)
    disp(strcat("Median of SML per virion = ", num2str(M)))

end

if plot_clus_r ==1
    Clustering_ratio(minSMLs_1, sample, SMLs_result, parent_path, sample_name)
end


function comp = GFPfilter(data, MaxInt) % delete all of the data with intensity above "MaxInt"
    comp = data;
    for i = 1:length(data)
        if data(i,6) >= MaxInt
            comp(i,1:5) = NaN;
        end
    end

end