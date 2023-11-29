function []= forcing_convolution_fun(MemoryLength, NameSuffix, ForcingName, ForcingLength, ForcingDir, SensitivityFieldName, ScaleFactor, WinterSensitivityDir, SpringSensitivityDir, SummerSensitivityDir, FallSensitivityDir)

% Yavor Kostov

% Computes the contribution of a given source of forcing (e.g.
% surface heat flux anomalies) to the variability in the objective function
% (e.g., the volume of Labrador Sea Water)

% Requires the gcmfaces Matlab toolbox
% for the ECCOv4 configuration of the MITgcm


%Input:

% MemoryLength - assumed memory span of the system with respect to past
% forcing

% NameSuffix - a suffix to be added in the name of the output file (e.g.,
% '6YearMemory')

% ForcingLength - total length of the forcing timeseries

% ForcingName - name of the forcing field files (e.g., 'oceQnet')

% ForcingDir - directory where the forcing files are stored

% WinterSensitivityDir - directory storing the sensitivity of the winter
% objective function

% SpringSensitivityDir - directory storing the sensitivity of the spring
% objective function 

% SummerSensitivityDir - directory storing the sensitivity of the summer
% objective function

% FallSensitivityDir - directory storing the sensitivity of the fall
% objective function

% SensitivityFieldName - name of the sensitivity files (e.g., 'ADJhflux')

% ScaleFactor - a constant for rescaling the sensitivity output
%

%Saved output stored in a file called
% strcat('cost_anom_', ForcingName, NameSuffix, '.mat') as follows:

% cost_anom_contrib - the contribution of a given source of forcing
% as a timeseries

% cost_anom_matr - the contribution of a given source of forcing at
% each timestep i to the anomaly in the objective function at timestep
% i+j-1. This matrix can be analyzed to identify preferred lead times.


global mygrid

FileNameSave=strcat('cost_anom_', ForcingName, NameSuffix, '.mat');

myenv.nctilesdir=ForcingDir;
listVars={ForcingName};
for vvv=1:length(listVars)
    vv=listVars{vvv};
    ForcingArray=read_nctiles([myenv.nctilesdir '/' vv],vv);
end

% We obtain the climatological seasonality of the forcing field in order to
% subtract it and compute anomalies relative to the seasonal cycle
for i=1:12
    ForcingArray_MON{i}=repmat(0*mygrid.Depth, [1 1 20]);
    ForcingArray_MON{i}=squeeze(ForcingArray(:,:,i:12:240));
end
for i=1:12
    ForcingArray_MON_clim{i}=nanmean(ForcingArray_MON{i},3);
end

cost_anom_matr=zeros(ForcingLength, ForcingLength);


for i=1:ForcingLength
    
    ForcingField=squeeze(ForcingArray(:,:,i));
    %We check what month of the year corresponds to forcing at time i:
    mon_num=mod(i,12)+12*(mod(i,12)==0);
    
    %We subtract the climatological monthly mean:
    monthly_climatology=squeeze(ForcingArray_MON_clim{mon_num});
    ForcingAnom=ForcingField-monthly_climatology;
    
    for j=1:MemoryLength
        
        if i+j-1<ForcingLength+1
            % We need to know the month of the year corresponding to a
            % lag of j-1 months following the forcing applied at time i:
            MONNUM2=mod(mon_num+(j-1),12)+12*(mod(mon_num+(j-1),12)==0);
            % We then check what season that month corresponds to
            % (April is considered a winter month):
            if (MONNUM2==1 || MONNUM2==2 || MONNUM2==3 || MONNUM2==4)
                SensitivityDirectory=WinterSensitivityDir;
                a=dir(strcat(SensitivityDirectory,'\', SensitivityFieldName, '.*.data'));
            elseif (MONNUM2==5 || MONNUM2==6)
                SensitivityDirectory=SpringSensitivityDir;
                a=dir(strcat(SensitivityDirectory,'\', SensitivityFieldName, '.*.data'));
            elseif (MONNUM2==7 || MONNUM2==8 || MONNUM2==9)
                SensitivityDirectory=SummerSensitivityDir;
                a=dir(strcat(SensitivityDirectory,'\', SensitivityFieldName, '.*.data'));
            else
                SensitivityDirectory=FallSensitivityDir;
                a=dir(strcat(SensitivityDirectory,'\', SensitivityFieldName, '.*.data'));
            end
            ab={a.name};
            abc=ab{length(ab)-j+1}; % the sensitivity output is ordered
                                    % from large lead times 
                                    % to lead time 0
            abcd=abc(10:end-5);
            [SensitivityMap]=ScaleFactor*rdmds2gcmfaces(strcat(SensitivityDirectory,'\', SensitivityFieldName,'.', abcd));
            cost_anom_product=ForcingAnom.*SensitivityMap;
            % The contribution of forcing in month i to the response in
            % month i+j-1 stored in a matrix:
            cost_anom_matr(i,i+j-1)=nansum(cost_anom_product(:))*24*30;
        end
        
    end
    
end

cost_anom_contrib=nansum(cost_anom_matr,1); % the contribution of the
                                            % forcing anomalies to the
                                            % objective function

save(FileNameSave, 'cost_anom_contrib', 'cost_anom_matr');