clear all
clc;

addpath(genpath('M:\Yen-Cheng\YCC_Matlab'))

%%%%%%%% Select and categorize files %%%%%%%%
[Myfiles,Mypath] = uigetfile({'*.dv';'*.csv';'*.mat';'*.*'}, 'Select files for all samples', 'MultiSelect','on');

temp = split(Myfiles, '_R3D.dv');
filebase = string(temp(1, :, 1));

temp2 = split(filebase, '-');
temp3 = [];
for i = 1 : length(filebase)
    temp3 = [temp3, strjoin(temp2(:, i, 1:(end-1)),'-')];
    disp(strjoin(temp2(:, i, 1:(end-1)),'-'))
end 
samplebase = unique(temp3(1, :));

%%% Create a cell for filename organization %%%
for j = 1:length(samplebase)
    a = 1;
    for i = 1:length(Myfiles)
        if strfind(temp{:,i,1}, samplebase(j))==1
            Org_files{j,a} = Myfiles{i};
            a=a+1;
        end
    end
end

%%%%%%%% Select and categorize files %%%%%%%%

% tfilebase = 'SVA_1B-'; 

matlaboutbase = strcat(Mypath,'output\');
mkdir(matlaboutbase);
addpath(genpath(Mypath)); % add the specified folder and subfolders to the path.

glob_SE = strel('diamond',1);

icrop=0;
plotfov = 1;
plotpart = 0;
onetime = 0;

corrsumsig_for_selection=300;  %% Change this parameter

for samp = 1:length(samplebase)
    for acq = 1 : length(Org_files(samp,:))

    %     if acq>=10
    %         tfilebase = 'SVA_1B-';
    %     end
    
        if exist(Org_files{samp,acq}, 'file') == 0
          % File does not exist
          % Skip to bottom of loop and continue with the loop
          continue;
        end
        temp = split(Org_files{samp,acq}, '_R3D.dv');
        tfilebase = char(temp{1});
        tsamplebase = split(tfilebase, '-');
        tsamplebase = strjoin(tsamplebase(1:(end-1)),'-');

            

       reader=bfGetReader(strcat(Mypath, Org_files{samp,acq}));
       I_GFP=bfGetPlane(reader,1);
       I_561=bfGetPlane(reader,2);
       I_647=bfGetPlane(reader,3);

    %    if icrop==1
    %       icropX=19;
    %       icropY=901;
    %       icropXhW=18;
    %       icropYhW=50;
    %       
    %       I_GFP=I_GFP(icropY-icropYhW:icropY+icropYhW, icropX-icropXhW:icropX+icropXhW);
    %    end

       %% ID GFP spots
       % part_GFP is
       % [x,y,meanback,sumsig,corrsumsig,SNR,resolvpow,perc_backpix]
       [ part_GFP ] = one_channel_partID( I_GFP, 3, corrsumsig_for_selection, [] );

       %% filter 488 spots ; create plots ; find intensities in 561 and 647 ; clearvars
       part_GFP_filt = part_GFP;
       %filter resolvpow>0.6
       part_GFP_filt = part_GFP_filt(part_GFP_filt(:,7)>0.6,:);
       %filter perc_backpix>40%
       part_GFP_filt = part_GFP_filt(part_GFP_filt(:,8)>40,:);
       %filter SNR>1.5
       part_GFP_filt = part_GFP_filt(part_GFP_filt(:,6)>1.5,:);
       %Discard top 15% brightest 488 spots
       tcorrsumsig = sort(part_GFP(:,5));

       tcut = tcorrsumsig(floor(0.85*size(tcorrsumsig,1)),1);
       part_GFP_filt = part_GFP_filt(part_GFP_filt(:,5)<tcut,:);

       if isempty(part_GFP_filt)
           continue
       end

       if plotfov
        ty = 1:1:size(part_GFP,1); ty=ty'; ty=100*ty/size(part_GFP,1);
        tcorrsumsig = sort(part_GFP(:,5));
        tSNR = sort(part_GFP(:,6));
        tresolvpow = sort(part_GFP(:,7));
        tperc_backpix = sort(part_GFP(:,8));

        tyf = 1:1:size(part_GFP_filt,1); tyf=tyf'; tyf=100*tyf/size(part_GFP_filt,1);
        tcorrsumsigf = sort(part_GFP_filt(:,5));
        tSNRf = sort(part_GFP_filt(:,6));
        tresolvpowf = sort(part_GFP_filt(:,7));
        tperc_backpixf = sort(part_GFP_filt(:,8));   

        h1 = figure;
        scatter(tcorrsumsig,ty); hold on
        scatter(tcorrsumsigf,tyf); hold off
        h2 = figure;
        scatter(tSNR,ty); hold on
        scatter(tSNRf,tyf); hold off
        h3 = figure;
        scatter(tresolvpow,ty); hold on
        scatter(tresolvpowf,tyf); hold off
        h4 = figure;
        scatter(tperc_backpix,ty); hold on
        scatter(tperc_backpixf,tyf); hold off
        h5 = figure;
        h5.Units = 'inches';
        h5.Position = [5 1 8 8];
        h5.PaperPosition = [0 0 8 8];         
        u1 = uipanel('Title','CorrSumSig','position',[0.01,0.51,0.48,0.48]);
        u2 = uipanel('Title','SNR','position',[0.51,0.51,0.48,0.48]);
        u3 = uipanel('Title','ResolvPow','position',[0.01,0.01,0.48,0.48]);
        u4 = uipanel('Title','% Bkgrnd pix','position',[0.51,0.01,0.48,0.48]);   
        set(get(h1,'Children'),'parent',u1);
        set(get(h2,'Children'),'parent',u2);
        set(get(h3,'Children'),'parent',u3);
        set(get(h4,'Children'),'parent',u4);
        close(h1,h2,h3,h4);
        figsavestr = strcat(matlaboutbase, tfilebase,'_spotstats.png');
        saveas(h5,figsavestr);
        close(h5);

        contI_r = contrast_stretch_forced( I_GFP, 4000, 8000 );

        contI_b = uint16(zeros(size(contI_r)));
        plotpart = part_GFP;
        for i = 1 : size(plotpart, 1)
            if and(plotpart(i, 2) == 0, plotpart(i, 1) == 0)
                continue;
            end
            contI_b(plotpart(i, 2), plotpart(i, 1)) = 2^16-1; 
        end
        contI_b = imdilate(contI_b,glob_SE);contI_b = imdilate(contI_b,glob_SE);
        contI_b = imdilate(contI_b,glob_SE);contI_b = imdilate(contI_b,glob_SE);
        contI_btemp = imerode(contI_b,glob_SE);
        contI_b = contI_b - contI_btemp;

        contI_g = uint16(zeros(size(contI_r)));
        plotpart = part_GFP_filt;
        for i = 1 : size(plotpart, 1)
            if and(plotpart(i, 2) == 0, plotpart(i, 1) == 0)
                continue;
            end
            contI_g(plotpart(i, 2), plotpart(i, 1)) = 2^16-1; 
        end
        contI_g = imdilate(contI_g,glob_SE);contI_g = imdilate(contI_g,glob_SE);
        contI_g = imdilate(contI_g,glob_SE);contI_g = imdilate(contI_g,glob_SE);
        contI_gtemp = imerode(contI_g,glob_SE);
        contI_g = contI_g - contI_gtemp;

        trgb_lowres = uint8( cat(3, contI_r+contI_b, contI_g, contI_b) ./ 256 );
        figpeaks = figure;
        figpeaks.Units = 'inches';
        figpeaks.Position = [5 1 8 8];
        figpeaks.PaperPosition = [0 0 8 8];
        subplot('Position', [0 0 1 1]);
        imshow(trgb_lowres);
        figsavestr = strcat(matlaboutbase, tfilebase,'_spots.png');
        saveas(figpeaks,figsavestr);
        close(figpeaks);
       end

       %% Find signals in other channels using 488 centroids
       % [x,y,meanback,sumsig,corrsumsig,SNR,resolvpow,perc_backpix]
       [ part_561_filt ] = one_channel_partID( I_561, 3, corrsumsig_for_selection, part_GFP_filt(:, 1:2) );
       [ part_647_filt ] = one_channel_partID( I_647, 3, corrsumsig_for_selection, part_GFP_filt(:, 1:2) );

       %% Extract signal data [id, x, y, 488, 561, 647]
       acqpartdata = zeros(size(part_GFP_filt,1), 6);
       acqpartdata(:, 1) = (1:size(part_GFP_filt, 1)) + acq*1000;
       acqpartdata(:, 2:3) = part_GFP_filt(:, 1:2);
       acqpartdata(:, 4) = part_GFP_filt(:, 5);
       acqpartdata(:, 5) = part_561_filt(:, 5);
       acqpartdata(:, 6) = part_647_filt(:, 5);


       %% Save datas
       if or(acq == 1, onetime)
           partdata = acqpartdata;
       else
           partdata = vertcat(partdata, acqpartdata);
       end

    end

%     clearvars -except partdata & Mypath & tsamplebase

    sortdata = zeros(size(partdata, 1), 4);
    sortdata(:, 1) = (1:size(partdata, 1)) * 1/size(partdata, 1);
    sortdata(:, 2) = sort(partdata(:, 4));
    sortdata(:, 3) = sort(partdata(:, 5));
    sortdata(:, 4) = sort(partdata(:, 6));
    
    plot(sortdata(:, 2), sortdata(:, 1),'DisplayName','GFP Intensity')
    hold on
    plot(sortdata(:, 4), sortdata(:, 1),'DisplayName','AF647 Intensity')
    hold off
    legend
    title(tsamplebase);


    % gfpsortdata = zeros(size(partdata, 1), 2);
    % gfpsortdata(:, 1) = partdata(:, 4);
    % gfpsortdata(:, 2) = partdata(:, 6);
    % [~,ord] = sort(gfpsortdata(:,1));
    % gfpsortdata = gfpsortdata(ord,:);
    % 
    % numfrac = round(size(partdata, 1) / 4);
    % offset = 4;
    % gfpsort_lo = gfpsortdata(1+offset:numfrac+offset, :);
    % gfpsort_hi = gfpsortdata(end-offset-numfrac+1:end-offset, :);

%     S1.partdata = partdata;
%     T1 = struct2table(S1);
%     out_fn1 = '_partdata.csv';
%     writetable(T1, strcat(Mypath, tsamplebase, out_fn1))
  
%     test =  vertcat(["partid" "x" "y" "488nm" "561nm" "640nm"], partdata);
%     empty = repmat(' ', length(partdata(:,1)), 1);
%     test2 = horzcat(partdata, empty, sortdata);
% 
%     S2.sortdata = sortdata;
%     T2 = struct2table(S2);
%     out_fn2 = '_sortdata.csv';
%     writetable(T2, strcat(Mypath, tsamplebase, out_fn2))

    Dataheader = ["partid","x","y","488nm","561nm","640nm"," ","sort frac","488nm","561nm","640nm"];
    writematrix(partdata, strcat(Mypath, tsamplebase, "raw.xlsx"), 'Range','A2')
    writematrix(sortdata, strcat(Mypath, tsamplebase, "raw.xlsx"), 'Range','H2')
    writematrix(Dataheader, strcat(Mypath, tsamplebase, "raw.xlsx"))

%     %%%%% Matlab 2019b doesn't support "AutoFitWidth" %%%%%%
%     title = ["partid","x","y","488nm","561nm","640nm"," ","sort frac","488nm","561nm","640nm"];
%     writematrix(title, strcat(Mypath, tsamplebase, "test.xlsx"),'AutoFitWidth',false)
%     writematrix(partdata, strcat(Mypath, tsamplebase, "test.xlsx"), 'Range','A2', 'AutoFitWidth', 0)
%     writematrix(sortdata, strcat(Mypath, tsamplebase, "test.xlsx"), 'Range','H2', 'AutoFitWidth', 0)
%     %%%%% Matlab 2019b doesn't support "AutoFitWidth" %%%%%%

end   

