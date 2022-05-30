 % Created July 2012 by Garrett Cobb
% Version 1.7; revised Jan 2013 Garrett Cobb
% Revised May 2021 to import files from Leica LASX - David DeWitt 
%
% This program is designed to fulfill all FCS data processing needs.
% 
% Inputs:
% -Accepts sin files in any mode (single, dual, quad).
% -If batch mode is used, accepts batches of files and automatically parses into datasets.
% -Automatically determines read range of correlation and intensity data,
% and adjusts for correlation measurements of different lengths.
% -Many different fit equations available.
% -Displays both correlation curves and intensity for inspection.
%
% Output files: 
% -ManyFits, which contains the fit parameters from all curves that were
% fit.
% -ManyCorrStdev, which contains each averaged correlation curve followed
% by its stdev, repeated for all of the datasets.
% -ManyStats, which contains intensity, %bleaching, etc.
% -FitData files for each individual dataset, which in addition to the fit
% parameters, contain columns of corrtimes, averaged corves, stdev, and
% fitted curve.

% Notes: 
% -mode 4 (quad mode, no averaging AxB and BxA) is not supported.
% -For future MATLAB releases the batchselect unique function will have to
% be updated.
% -path variable is not passed to import_inspect, which is why it can only
% work in the current path.

% add triplet + background fit equation sometime
% problem sometimes with parsing similar names in batch mode


function fcs_omega

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EDIT HERE to set fitting equation (or 0 for no fitting).

% 0 - No fitting
% 1 - 3D diffusion
% 2 - 3D diffusion fixs
% 3 - 2D diffusion
% 4 - 3D anomalous diffusion fixs
% 5 - 3D twocomponent fixs
% 6 - 3D triplet fixs
% 7 - 2D triplet
% 8 - 3D diffusion fixs background
% 9 - 3D diffusion anomalous fixs background
% 10 - 3D triplet background fixs
% 11 - 3D twocomponent brightnessfactor fixs
% 12 - 3D twocomponent free+anomalous fixs
eqnmode = 9;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EDIT HERE to select appropriate read ranges. For a single component, 
% do ~1.5-2 decades in each direction from expected diffusion time.
% for free fluorophore recommend .001 ms to 10ms
% for labeled protein recommend .01 ms to 100ms
% for in vivo measurements recommend .01 ms to 1000ms

tau_begin =    0.01    ;   % milliseconds
tau_end   =    1000      ;   % milliseconds


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EDIT HERE to set initial parameters, and upper and lower bounds.

%      number  tautriplet  taudiff1  taudiff2    A    alpha    G      s   alpha1   alpha2     
fp = [   20      0.1         0.5     100       0.5    0.5     0     0.17   0.5         0.5;     % initial values
         0       0.001       0.02    1.5       0     0.1   -.003   0.1    0         0.1;      % lower bounds
       10000     1.2         2        500        1      1     .003   0.5    1         5];     % upper bounds% End of recommended user edits.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global s
%s = fp(1,end);
s = fp(1, 8);
disp(s)
fpstr =[{'number   '} 'tautriplet  ' 'taudiff1  ' 'taudiff2  ' '   A    '  ' alpha  '  '   G    ' '   s    ' '   alpha1   ' '   alpha2   ' '  Int   '];       % strings indicating the variables used
op=optimset('MaxFunEvals', 100000, 'MaxIter', 5000, 'TolFun', 1.0000e-018);

% Ask for user input to decide between batch mode and dataset mode
%disp('How many data sets do you have to process?')
%inputdatasets=input('Or enter 0 for batch mode. ');
%if inputdatasets==0
%    [filelistbatch filelist3 filelist4 filelist5 datasets path] = fcs_batchselect; %#ok<*ASGLU>
%elseif inputdatasets > 0
%    datasets=inputdatasets;
%end


% Setup for fitting. Program will change settings based on eqnmode.

% number  tautriplet  taudiff1  taudiff2   A   alpha    G      s    alpha1,alph2    
%   1         2         3          4       5     6      7      8      9,10


if eqnmode == 0
    elseif eqnmode==1
        diffequation = @ FCS_3Ddiffusion;
        disc = [2 4:7 9:10];
        fp(:,disc)=[];fpstr(:,disc)=[];
    elseif eqnmode==2
        diffequation = @ FCS_3Ddiffusion_fixs;
        disc = [2 4:7 9:10];
        fp(:,disc)=[];fpstr(:,disc)=[];
    elseif eqnmode==3
        diffequation = @ FCS_2Ddiffusion;
        disc = [2 4:10];
        fp(:,disc)=[];fpstr(:,disc)=[];
    elseif eqnmode==4
        diffequation = @ FCS_anomalousdiffusion_fixs;
        disc = [2 4 5 7 9:10];
        fp(:,disc)=[];fpstr(:,disc)=[];
    elseif eqnmode==5
        diffequation = @ FCS_3Dtwocomponent_fixs;
        disc = [2 6 7 9:10];
        fp(:,disc)=[];fpstr(:,disc)=[];
    elseif eqnmode==6
        diffequation = @ FCS_3Dtriplet_fixs;
        disc = [4 6 7 9:10];
        fp(:,disc)=[];fpstr(:,disc)=[];
    elseif eqnmode==7
        diffequation = @ FCS_2Dtriplet;
        disc = [3 6:10];
        fp(:,disc)=[];fpstr(:,disc)=[];
    elseif eqnmode==8
        diffequation = @ FCS_3Ddiffusion_fixs_background;
        disc = [2 4:6 9:10];
        fp(:,disc)=[];fpstr(:,disc)=[];
    elseif eqnmode==9
        diffequation = @ FCS_anomalousdiffusion_fixs_background;
        disc = [2 4 5 9:10];
        fp(:,disc)=[];fpstr(:,disc)=[];
    elseif eqnmode==10
        diffequation = @ FCS_3D_triplet_background_fixs;
        disc = [4 6 9:10];
        fp(:,disc)=[];fpstr(:,disc)=[];
    elseif eqnmode==11
        diffequation = @ FCS_3D_twocomponent_brightnessfactor_fixs;
        disc = [2 7 9:10];
        fp(:,disc)=[];fpstr(:,disc)=[];
    elseif eqnmode==12
        diffequation = @ FCS_3Dtwocomponent_fixs_withalpha;
        disc = [2 6 7 9];
        fp(:,disc)=[];fpstr(:,disc)=[];
        
end

    

% declare variables before the first loop
global names1; names1 = []; % label for columns of ManyStats file
global names2; names2 = []; % label for columns of ManyCorrStdev file
global names3; names3 = []; % label for rows of ManyStats file
names4={}; % label for rows of ManyFits file


    
% Assemble the filelist to pass to fcs_import_inspect - either
% automatically from the next dataset in the batch, or by asking the
% user for the next dataset.
%if inputdatasets==0
%    filelist = filelistbatch(:,filelist4(k):filelist5(k));
%elseif inputdatasets > 0
    disp('Press spacebar to select Correlation data'); pause
    filefilter = {'*.xlsx','excel files';'*.csv','csv files';'*.*','All files'};
    [filelist, path]=uigetfile(filefilter,'choose files','Multiselect', 'on' );        
    try
        filelist{1};  %check if the filelist is stored as a cell array  
    catch
        temp=filelist;
        clear filelist
        filelist{1}=temp;  % force filelist to be a cell array
    end
    disp('Press spacebar to select Intensity data'); pause
    filefilter = {'*.xlsx','excel files';'*.csv','csv files';'*.*','All files'};
    [filelistInt, path]=uigetfile(filefilter,'choose files','Multiselect', 'on' );
    try
        filelistInt{1};
    catch
        temp1=filelistInt;
        clear filelistInt
        filelistInt{1}=temp1;
    end
%end
    
    % Call fcs_import_inspect.
    [name corrtimes datasets corrdatafull inttimes intdatafull tracenum filenumber T P Rep Ts Ps Reps] = fcs_import_inspect(filelist,filelistInt,path,tau_begin,tau_end);
    
for k=1:T    
  for k2=1:P
    corrdata=corrdatafull(:,1:Rep);
    intdata=intdatafull(:,1:Rep);
    corrdatafull(:,1:Rep)=[];
    intdatafull(:,1:Rep)=[];
    
    % Plot correlation curves and intensity data.
    set(0,'DefaultAxesLineStyleOrder','-|:|--')
    scrsz = get(0,'ScreenSize');
    topleft = [5 scrsz(4)/20 scrsz(3)*3/4 scrsz(4)*7/8];
    h=figure('OuterPosition',topleft);  % h is the figure handle

    % mode-dependent plotting - plot 1, 2, or 4 correlation functions.
        %hy is the corr subplot handle
        hy=subplot(4,1,[1 3]);semilogx(corrtimes,corrdata(:,:,1))
        legend('show'),ylabel('G(\tau) Channel A'),xlabel('\tau (ms)')  % this subplot is corrdata
        fullinta=reshape(intdata(:,:,1),[],1);  % produce the fullint
        hz=subplot(4,1,4); plot(1/tracenum:1/tracenum:Rep,fullinta)
        %hz is the intensity subplot handle
        ylabel('Intensity, Channel A'),xlabel('FileNumber')
        set(hz,'XGrid','on')


    % Remove aberrant curves from data matrix

    realdiscardlist=[];
    discard1list=[];
    j=1;
    mode=1;
    numcorr=1;
    numint=1;
    for j=1:Rep
        discard1=input('Which curve to discard (0 if none)?');
        while isempty(discard1) || (discard1 < 0) || (discard1 > 10 )
            disp('are you sure?')
            discard1=input('Which curve to discard (0 if none)?');
        end
        figure(h)
        if discard1~=0
            corrdata(:,discard1,:)=[]; % drop the curve from corrdata
            intdata(:,discard1,:)=[]; % drop the curve from intdata
            filelist{discard1}=[]; % drop the filename from the list
            realdiscard=sum(discard1list<=discard1)+discard1;
            discard1list(j)=discard1; %#ok<AGROW>
            while sum(realdiscard==realdiscardlist)>=1;
                realdiscard=realdiscard+1;
            end
            realdiscardlist(j)=realdiscard; %#ok<AGROW>
    
                axes(hy);subplot(4,1,[1 3]);semilogx(corrtimes,corrdata(:,:,1))
                legend('show'),ylabel('G(\tau) Channel A'),xlabel('\tau (ms)')  % this subplot is corrdata
                fullinta=reshape(intdata(:,:,1),[],1);  % produce the fullint
                %axes(hz);subplot(4,1,4);plot(1/tracenum:1/tracenum:(Rep-j),fullinta)
                hz=subplot(4,1,4); plot(1/tracenum:1/tracenum:(Rep-j),fullinta)
                %  hz is the intensity subplot handle
                ylabel('Intensity, Channel A'),xlabel('FileNumber')
                set(hz,'XGrid','on')    
        elseif discard1==0
            break;
        end
    end
    
    
    realdiscardlist=sort(realdiscardlist);
    
    % subtract 1
    corrdatasub=corrdata;
    r=size(corrdatasub,2);
    
    if (P*(k-1)+k2)==1 % If this is the first runthru... 
        % pre-allocate the corr_stdev matrix.
        corr_stdev = zeros(length(corrtimes),datasets*mode*2);
        % remember the mode of the initial files
        allmode = mode;
        % Set names3 to go with the ManyStats file
        if mode==1   
            names3={'Int';'% Int';'Agr1';'Agr2';'Agr3';'Noise';'Discards'};
        elseif mode==2
            names3={'Red Int';'Green Int';'% Red';'% Green';'Red Agr1';'Red Agr2';'Red Agr3';...
                'Red Noise';'Green Agr1';'Green Agr2';'Green Agr3';'Green Noise';'Discards';};
        elseif mode==3
            names3={'Red Int';'Green Int';'% Red';'% Green';'Red Agr1';'Red Agr2';...
                'Red Agr3';'Red Noise';'Green Agr1';'Green Agr2';'Green Agr3';'Green Noise';...
                'XCorr Agr1';'XCorr Agr2';'XCorr Agr3';'XCorr Noise';'Discards';};
        end
    end
    
    if allmode ~=mode
        disp('**Error!** All files analyzed in a single program run must be same mode (i.e. single, dual, or quad)')
        break
    end
    
    % compute avgcorr and sd for the (one or two) autocorrelation functions
    for i=1:mode
        corr_stdev(1:length(corrtimes),(i*2-1+((P*(k-1)+k2)-1)*mode*2) )= mean(corrdatasub(:,:,i), 2);
        if r>1
            stdtemp=std (corrdatasub(:,:,i),0,2);
            stdtemp(stdtemp==0)=min(nonzeros(stdtemp)); % replace any zero stdev with a small non-zero value. zero stdev causes fitting errors later.
        else
            stdtemp=ones(size(corrtimes));
        end  
        corr_stdev(1:length(corrtimes),(i*2  +((P*(k-1)+k2)-1)*mode*2) )= stdtemp;
    end
    
    % now compute avgcorr and sd for the xcorr functions.
    if mode==3
        allxcorr=reshape(corrdatasub(:,:,3:4),length(corrtimes),[],1);
        corr_stdev(1:length(corrtimes),(5+((P*(k-1)+k2)-1)*mode*2) )= mean(allxcorr, 2);
        if r>1
            stdtemp=std (allxcorr,0,2);
            stdtemp(stdtemp==0)=min(nonzeros(stdtemp)); % replace any zero stdev with a small non-zero value.
        else
            stdtemp=ones(size(corrtimes));
        end
        corr_stdev(1:length(corrtimes),(6+((P*(k-1)+k2)-1)*mode*2) )= stdtemp;
    end
    
    % Calculate agreement and noise for each
    for i=1:mode
        agrnoise((i*4-3):i*4) = ...  agreement is StDev between curves, normalized by y.
            ... And averaged over a small window to get a better measure.  Agreement is a function of tau
            [ mean(corr_stdev(6:13,(i*2  +((P*(k-1)+k2)-1)*mode*2)))/mean(corr_stdev(6:13,(i*2-1+((P*(k-1)+k2)-1)*mode*2)))... at intercept (.03-.06ms)
            mean(corr_stdev(44:48,(i*2  +((P*(k-1)+k2)-1)*mode*2)))/mean(corr_stdev(44:48,(i*2-1+((P*(k-1)+k2)-1)*mode*2)))... at 1ms
            mean(corr_stdev(71:74,(i*2  +((P*(k-1)+k2)-1)*mode*2)))/mean(corr_stdev(71:74,(i*2-1+((P*(k-1)+k2)-1)*mode*2)))... at 10ms
            std(corr_stdev(32:42,(i*2-1+((P*(k-1)+k2)-1)*mode*2)))/mean(corr_stdev(32:42,(i*2-1+((P*(k-1)+k2)-1)*mode*2)))];
        % noise is the stdev of the correlation in the intercept region (at t={.3 .75 ms}),
        % normalized by the magnitude
    end

    % Calculate average intensity and intensity stats
    intafterdiscard=reshape(intdata,[],numint);
    avgint  = mean(intafterdiscard); %works for one or two channels of intensity
    for i=1:numint
        normint(:,i) = intafterdiscard(:,i)./ mean(intafterdiscard(1:25,i)); %#ok<*AGROW> %entire series, normalized
    end
    endint  = mean (normint((end-25):end,:)); %fraction of intensity remaining
    comment=[avgint endint agrnoise realdiscardlist]';
    stats(1:length(comment),(P*(k-1)+k2))=comment;
    
    % Plot average curves
    if mode==1
        % Plot final avgcorr and intensity
        hy=subplot(4,1,[1 3]);semilogx(corrtimes,corr_stdev(:,(1*2-1+((P*(k-1)+k2)-1)*mode*2)))
        ylabel('G(\tau) Channel A'),xlabel('\tau (ms)')  % this subplot is corrdata
        fullinta=reshape(intdata(:,:,1),[],1);  % produce the fullint
        hz=subplot(4,1,4); plot(1/tracenum:1/tracenum:(Rep-j+1),fullinta)
        %  hz is the intensity subplot handle
        ylabel('Intensity, Channel A'),xlabel('FileNumber')
        set(hz,'XGrid','on')
    elseif mode==2
        % Plot final avgcorr and intensity
        hy=subplot(4,2,[1 5]);semilogx(corrtimes,corr_stdev(:,(1*2-1+((P*(k-1)+k2)-1)*mode*2)))
        ylabel('G(\tau) Channel A'),xlabel('\tau (ms)')  % this subplot is corrdata
        he=subplot(4,2,[2 6]);semilogx(corrtimes,corr_stdev(:,(2*2-1+((P*(k-1)+k2)-1)*mode*2)))
        ylabel('G(\tau) Channel B'),xlabel('\tau (ms)')  % this subplot is corrdata
        fullinta=reshape(intdata(:,:,1),[],1);  % produce the fullint
        fullintb=reshape(intdata(:,:,2),[],1);
        hz=subplot(4,2,7); plot(1/tracenum:1/tracenum:(Rep-j+1),fullinta)
        %  hz is the intensity subplot handle
        ylabel('Intensity, Channel A'),xlabel('FileNumber')
        set(hz,'XGrid','on')
        hf=subplot(4,2,8); plot(1/tracenum:1/tracenum:(Rep-j+1),fullintb)
        %hz is the intensity subplot handle
        ylabel('Intensity, Channel B'),xlabel('FileNumber')
        set(hf,'XGrid','on')
    elseif mode==3
        % Plot final avgcorr and intensity
        hy=subplot(4,3,[1 7]);semilogx(corrtimes,corr_stdev(:,(1*2-1+((P*(k-1)+k2)-1)*mode*2)))
        ylabel('G(\tau) Channel A'),xlabel('\tau (ms)')  % this subplot is corrdata
        he=subplot(4,3,[2 8]);semilogx(corrtimes,corr_stdev(:,(2*2-1+((P*(k-1)+k2)-1)*mode*2)))
        ylabel('G(\tau) Channel B'),xlabel('\tau (ms)')  % this subplot is corrdata
        fullinta=reshape(intdata(:,:,1),[],1);  % produce the fullint
        fullintb=reshape(intdata(:,:,2),[],1);
        hz=subplot(4,3,10); plot(1/tracenum:1/tracenum:(Rep-j+1),fullinta)
        %  hz is the intensity subplot handle
        ylabel('Intensity, Channel A'),xlabel('FileNumber')
        set(hz,'XGrid','on')
        hf=subplot(4,3,11); plot(1/tracenum:1/tracenum:(Rep-j+1),fullintb)
        %hz is the intensity subplot handle
        ylabel('Intensity, Channel B'),xlabel('FileNumber')
        set(hf,'XGrid','on')
        pa=subplot(4,3,[3 6]);semilogx(corrtimes,corr_stdev(1:length(corrtimes),(5+((P*(k-1)+k2)-1)*mode*2) ))
        pb=subplot(4,3,[9 12]);semilogx(corrtimes,corr_stdev(1:length(corrtimes),(5+((P*(k-1)+k2)-1)*mode*2) ));cla(pb)
    end
  
    % Begin fitting. Change settings based on eqnmode.
    
    if eqnmode==0
        disp('No fitting.')
    elseif eqnmode>0
        disp('Now fitting curves with function:')
        disp(diffequation)
        
        for i=1:mode % loop for fitting mutiple times per dataset (dual or quad mode)
            
            % Define the function to minimize
            weighted_residuals = @(x) (diffequation(x,corrtimes)-corr_stdev(1:length(corrtimes),(i*2-1+((P*(k-1)+k2)-1)*mode*2) ))... this term is avgcorr
                ./corr_stdev(1:length(corrtimes),(i*2  +((P*(k-1)+k2)-1)*mode*2) ); % this term is stdev
            
            % Fit using least squares of weighted_residuals.
            x = lsqnonlin(@(x) weighted_residuals(x), fp(1,:), fp(2,:), fp(3,:), op);
            fitdata = diffequation(x,corrtimes);  % The best fit curve using the optimized parameters in x
            
            % Get results
            resrow=i+((P*(k-1)+k2)-1)*mode;
            result(resrow,1:size(x,2))= x; % x can be of different sizes, depending on the fit equation
            
            if or(i==1,i==2) % Get correct intensity to go with ACF, or none to go with XCF 
                result(resrow,size(x,2)+1)= avgint(i);  
            elseif i==3
                result(resrow,size(x,2)+1)= 0;
            end
            
            % Display data, fitdata & residuals.
            figure
            subplot(4,1,[1 3]), semilogx(corrtimes,corr_stdev(1:length(corrtimes),(i*2-1+((P*(k-1)+k2)-1)*mode*2) ),'.b',corrtimes,fitdata,'-r')
            infobox{1} = cell2mat(fpstr);
            infobox(2) = {num2str(result(resrow,:),'%-12.4g')};
            infobox(3) = {['Fitting from ' name '_T' num2str(Ts(k)) 'P' num2str(Ps(k2))]};
            annotation('textbox', [.5 .6 .4 .3], 'String', infobox,'interpreter','none','FitBoxToText','on');
            ylabel('G(\tau)')
            title(['Correlation ' num2str(i)])
            subplot(4,1,4), semilogx(corrtimes,(corr_stdev(1:length(corrtimes),(i*2-1+((P*(k-1)+k2)-1)*mode*2) )-fitdata)./corr_stdev(1:length(corrtimes),(i*2  +((P*(k-1)+k2)-1)*mode*2) ),'.',corrtimes,0,'-r')
            xlabel('\tau (ms)')
            ylabel('Fitting Residuals (\sigma)')
            
            % Write a FitData file for this interval.
            alldata=[corrtimes corr_stdev(1:length(corrtimes),(i*2-1+((P*(k-1)+k2)-1)*mode*2)) corr_stdev(1:length(corrtimes),(i*2+((P*(k-1)+k2)-1)*mode*2) ) fitdata];
            % Compute chi-square
            chi_square = sum(((alldata(:,4)-alldata(:,2))./alldata(:,3)).^2);
            result(resrow,size(x,2)+2)= chi_square;
            
            data=[path 'FitData_mode_' num2str(eqnmode) '_' name '_T' num2str(Ts(k)) 'P' num2str(Ps(k2)) 'corr' num2str(i) '.txt'];
            dlmwrite(data, result(resrow,:), '\t');
            dlmwrite(data, alldata, '-append','roffset', 1, 'delimiter', '\t')
            
        end

    end % end fitting portion
    clear normint intafterdiscard intdata
    
    str=[ 'Time' num2str(Ts(k)) ' Position' num2str(Ps(k2)) ' completed.'];
    disp(str)
    disp(' ')
    names1{P*(k-1)+k2}=[name '_T' num2str(Ts(k)) 'P' num2str(Ps(k2))];
    for hea=1:mode
        names2{1+((P*(k-1)+k2)-1)*mode*2+(hea-1)*2}=[name 'T' num2str(Ts(k)) 'P' num2str(Ps(k2))];
        names2{2+((P*(k-1)+k2)-1)*mode*2+(hea-1)*2}='sd';
        names4 = [names4 [name '_T' num2str(Ts(k)) 'P' num2str(Ps(k2))]];
    end
  end  
end

% xlswrite is very slow! but seems to be the easiest way to write a combination of strings and numeric.

if eqnmode>0
    try
        xlswrite([path 'ManyFits.xlsx'], ones(40)+NaN); % clear the shreadsheet. Writing NaN writes a blank cell.
        xlswrite([path 'ManyFits.xlsx'], [{''} fpstr]) 
        xlswrite([path 'ManyFits.xlsx'], [names4' num2cell(result)],'Sheet1','A2') 
    catch
        msgbox('Paused for error! Cannot write your data while you have ManyFits.xlsx open! Close the document and press any key to resume','Error','warn')
        disp('Paused for error! Cannot write your data while you have ManyFits.xlsx open!')
        disp('Close the document and press any key to resume')
        pause
        xlswrite([path 'ManyFits.xlsx'], ones(40)+NaN); % clear the shreadsheet. Writing NaN writes a blank cell.
        xlswrite([path 'ManyFits.xlsx'], [{''} fpstr]) 
        xlswrite([path 'ManyFits.xlsx'], [names4' num2cell(result)],'Sheet1','A2') 
    end
end

data=[path 'CorrTime.txt'];
dlmwrite(data, corrtimes, 'delimiter', '\t')
data=[path 'ManyCorrStdev.txt'];
dlmwrite(data, corr_stdev, 'delimiter', '\t')
data=[path 'ManyStats.txt'];
dlmwrite(data, stats, 'delimiter', '\t')


disp(' ')                                             
disp('Alright; finished with these datasets:')
disp(names1')
disp('Goodbye!')
fclose('all');

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fcs_batchselect
function [filelistbatch filelist3 filelist4 filelist5 datasets path] = fcs_batchselect

% import batch and make a list of files
disp('Press spacebar to select the whole batch of files'); pause
filefilter = {'*.sin','correlator files';'*.tsv','flexx files';'*.*','All files'};
[filelistbatch, path]=uigetfile(filefilter,'choose files','Multiselect', 'on' );

if filelistbatch{1}==0
    disp('**Error! No files chosen')
end

% regular expression for the end of the filename (_##.sin)
sinname='(_|-)(\w\w|\w)(.sin|.SIN)';

% get the part of the filename that doesnt include (_##.sin)
filelist1=regexp(filelistbatch,sinname); %filelist1 says what is the last character of the root name
for i=1:length(filelistbatch)
   filelist2{i}=filelistbatch{i}(1:filelist1{i}); %#ok<AGROW> %filelist2 is the list of root names
end

[filelist3,filelist5,null]=unique(filelist2,'legacy');%filelist3 is unique names
filelist5 = filelist5';
% filelist5 is the end index of each unique dataset

datasets=length(filelist3);

% filelist 4 is the start index of each unique dataset
for i=1:datasets
    filelist4(i)=find(strcmp(filelist2,filelist3{i}),1,'first'); %#ok<AGROW>
end

str = ['Here are the ', num2str(datasets), ' datasets:'];
disp(' ')
disp(str)
disp(filelist3')
disp('Press any key to continue'); pause
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%fcs_import_inspect
function [name corrtimes datasets corrdatafull inttimes intdatafull tracenum filenumber T P Rep Ts Ps Reps] = fcs_import_inspect(filelist,filelistInt,path,tau_begin,tau_end)

%[filenumber]=size(filelist,2);
%sinname='(-|_)(\w\w|\w)(.sin|.SIN|.fcs|.xlsx)';
%name1=regexp(filelist{1},sinname);
name=filelist{1};
if isempty(name);
    name=filelist{1}(1:end-5);
end
disp(' ')
str=['Now analyzing ' name];
disp(str)

% data_file = [path filelist{1}];
% corrtimeinit=dlmread(data_file,',', [15 0 190 0]);

data_file = [path filelist{1}];
data_file_int = [path filelistInt{1}];
headertable=readtable(data_file_int, 'VariableNamingRule', 'preserve', 'Range', '1:1');
headervars=headertable.Properties.VariableNames;
%header=char(headervars(length(headervars)-1));
headlength=(length(headervars))/2;
TPRpts=zeros(3,headlength);
for hea=1:headlength
    TPRpoint=sscanf(char(headervars(2*hea-1)),'T%d P%d R%d');
    TPRpts(:,hea)=TPRpoint;
end
Ts=unique(TPRpts(1,:));
Ps=unique(TPRpts(2,:));
Reps=unique(TPRpts(3,:));
disp(' ') 
str2=['Found: ' num2str(length(Ts)) ' Timepoint(s), ' num2str(length(Ps)) ' Position(s), and ' num2str(length(Reps)) ' Repeat(s) at each position'];
disp(str2)
fullsinfile=readtable(data_file, 'VariableNamingRule', 'preserve', 'VariableNamesRange', 'A2');
fullintfile=readtable(data_file_int, 'VariableNamingRule', 'preserve', 'VariableNamesRange', 'A2');
%fclose('all');
% tracenum1 = fullsinfile{1}(193); % this was how it used to be done
% tracenum=str2double(tracenum1{1}(11:end))-1;

%keyphrase = '##NPOINTS='; % phrase preceding the number of points for the correlation data
%skeyphrase = size(keyphrase);
%correlation_numpoints_line=0; % typically 14
%satisfied = 0;
%while (~satisfied)
%    try % just in case our test line is shorter than our keyphrase
%        correlation_numpoints_line=correlation_numpoints_line+1;
%        testline = fullsinfile{1}(correlation_numpoints_line);
%        satisfied = strcmp(testline{1}(1:skeyphrase(2)),keyphrase);
%    catch
%        correlation_numpoints_line=correlation_numpoints_line+1;
%    end
%end
%ncorrpoints = str2double(testline{1}(skeyphrase(2)+1:end)); % typically 176

%keyphrase = '##NPOINTS='; % phrase preceding the number of points for the intensity data
%skeyphrase = size(keyphrase);
%intensity_numpoints_line=correlation_numpoints_line+ncorrpoints; % typically 190
%satisfied = 0;
%while (~satisfied)
%    try % just in case our test line is shorter than our keyphrase
%        intensity_numpoints_line=intensity_numpoints_line+1;
%        testline = fullsinfile{1}(intensity_numpoints_line);
%        satisfied = strcmp(testline{1}(1:skeyphrase(2)),keyphrase);
%    catch
%        intensity_numpoints_line=intensity_numpoints_line+1;
%    end
%end
%nintpoints = str2double(testline{1}(skeyphrase(2)+1:end));


%range1begin = correlation_numpoints_line+1; % +find(corrtimeinit>(tau_begin),1,'first'); % this +1 comes from the successive ##XYPOINTS line
%range1end = range1begin+ncorrpoints-1; % +find(corrtimeinit<(tau_end),1,'last');
%range2begin = intensity_numpoints_line+2; % this +2 comes from the fact that ##DELTAX and ##XYPOINTS lines follow ##NPOINTS
%range2end = range2begin+tracenum-1;

% Get time data from first file.
corrsize=size(fullsinfile,1);
intsize=size(fullintfile,1);
corrtimesfull = table2array(fullsinfile(:,'Time [ms]'));
corrtimetrim = find(corrtimesfull>tau_begin & corrtimesfull<tau_end);
corrtimes=corrtimesfull(corrtimetrim(1):corrtimetrim(end));
inttimes = table2array(fullintfile(1:intsize-50,'Time'));
[columnnumber]=size(fullintfile,2);
[filenumber]=columnnumber/2;
T=length(Ts);
P=length(Ps);
Rep=length(Reps);
datasets=T*P;
tracenum = size(inttimes,1);
% Change time units from seconds to milliseconds.

    
% Get corrdata and intdata from each file.
corrdatafull=zeros(size(corrtimes,1),filenumber,1);
intdatafull=zeros(size(inttimes,1),filenumber,1);
corrdatafull(:,1,:)=table2array(fullsinfile(corrtimetrim(1):corrtimetrim(end),'Correlation Channel 1'));
intdatafull(:,1,:)=table2array(fullintfile(1:intsize-50,'Count Rate Channel 1 [kCounts/s]'));
str1='Correlation Channel 1_';
str2='Count Rate Channel 1 [kCounts/s]_';
for i=1:(filenumber-1)
    try
    str3=append(str1,num2str(i));
    str4=append(str2,num2str(i));
    corrdatafull(:,i+1,:)=table2array(fullsinfile(corrtimetrim(1):corrtimetrim(end),str3));
    intdatafull(:,i+1,:)=table2array(fullintfile(1:intsize-50,str4));    
    end
%    data_file = filelist{i};
%    try
%    corrdata(:,i,:)=dlmread(data_file, ',', [range1begin 1 range1end 1]);
%    intdata(:,i,:)=dlmread(data_file, ',', [range2begin 1 range2end 1]);  
%    catch
%        filelist{i}
%    end
end


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%fitting functions

function F=FCS_3Ddiffusion(x,corrtimes)
number=x(1);
taudiff1=x(2);
s=x(3);

F=((1/number)*(1./(1+corrtimes/taudiff1)).*sqrt(1./(1+s*s*corrtimes/taudiff1)));
%returns F, which is G(tau), using the current values of x
end

function F=FCS_3Ddiffusion_fixs(x,corrtimes)
global s
number=x(1);
taudiff1=x(2);


F=((1/number).*(1./(1+corrtimes./taudiff1)).*sqrt(1./(1+s.*s.*corrtimes./taudiff1)));
%returns F, which is G(tau), using the current values of x
end

function F=FCS_2Ddiffusion(x,corrtimes)
number=x(1);
taudiff1=x(2);

F=((1/number)*(1./(1+corrtimes/taudiff1)));
%returns F, which is G(tau), using the current values of x
end

function F=FCS_anomalousdiffusion_fixs(x,corrtimes)
global s
number=x(1);
taudiff1=x(2);
alpha=x(3);

F=((1/number)*(1./(1+corrtimes/taudiff1).^alpha).*sqrt(1./(1+s*s*(corrtimes/taudiff1).^alpha)));
%returns F, which is G(tau), using the current values of x
end

function F=FCS_3Dtwocomponent_fixs(x,corrtimes)
global s
number=x(1);
taudiff1=x(2);
taudiff2=x(3);
A=x(4);

F=((1/number)*(A*(1./(1+corrtimes/taudiff1)).*sqrt(1./(1+s*s*corrtimes/taudiff1))...
    +((1-A)*1./(1+corrtimes/taudiff2)).*sqrt(1./(1+s*s*corrtimes/taudiff2))));
end

function F=FCS_3Dtwocomponent_fixs_withalpha(x,corrtimes)
global s
number=x(1);
taudiff1=x(2);
taudiff2=x(3);
A=x(4);
alpha2=x(6);

F=((1/number)*(A*(1./(1+corrtimes/taudiff1)).*sqrt(1./(1+s*s*corrtimes/taudiff1))...
    +((1-A)*1./(1+(corrtimes/taudiff2).^alpha2)).*sqrt(1./(1+s*s*(corrtimes/taudiff2).^alpha2))));
end

function F=FCS_3Dtriplet_fixs(x,corrtimes)
global s
number=x(1);
tautriplet=x(2);
taudiff1=x(3);
A=x(4);

F=((1/number)*(1./(1-A)).*(1.-A+A*exp(-corrtimes/tautriplet))...
    .*(1./(1+corrtimes/taudiff1)).*sqrt(1./(1+s*s*corrtimes/taudiff1)));
end

function F=FCS_2Dtriplet(x,corrtimes)

number=x(1);
tautriplet=x(2);
taudiff2=x(3);
A=x(4);

F=((1/number)*(1./(1-A)).*(1.-A+A*exp(-corrtimes/tautriplet)).*(1./(1+corrtimes/taudiff2)));
end

function F=FCS_3Ddiffusion_fixs_background(x,corrtimes)
global s
number=x(1);
taudiff1=x(2);
G=x(3);

F=((1/number).*(1./(1+corrtimes./taudiff1)).*sqrt(1./(1+s.*s.*corrtimes./taudiff1))+G);
%returns F, which is G(tau), using the current values of x
end

function F=FCS_anomalousdiffusion_fixs_background(x,corrtimes)
global s
number=x(1);
taudiff1=x(2);
alpha=x(3);
G=x(4);

F=((1/number)*(1./(1+(corrtimes/taudiff1).^alpha)).*sqrt(1./(1+s*s*(corrtimes/taudiff1).^alpha)))+G;
%returns F, which is G(tau), using the current values of x
end

function F=FCS_3D_triplet_background_fixs(x,corrtimes)
global s
number=x(1);
tautriplet=x(2);
taudiff1=x(3);
A=x(4);
G=x(5);

F=((1/number)*(1./(1-A)).*(1.-A+A*exp(-corrtimes/tautriplet)).*(1./(1+corrtimes/taudiff1)).*sqrt(1./(1+s*s*corrtimes/taudiff1)))+G;
end

function F=FCS_3D_twocomponent_brightnessfactor_fixs(x,corrtimes)
global s
number=x(1);
taudiff1=x(2);
taudiff2=x(3);
A=x(4);
alpha=x(5);

F=((1/number)*(A*(1./(1+corrtimes/taudiff1)).*sqrt(1./(1+s*s*corrtimes/taudiff1))...
    +alpha*((1-A)*1./(1+corrtimes/taudiff2)).*sqrt(1./(1+s*s*corrtimes/taudiff2))));
end




