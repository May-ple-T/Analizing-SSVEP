function makeSSVEP

% define data size
channels = 8;
max_t_length = 128 * 5;
sets = 20;
targetNum = 4;


data = [load('yoshimoto_4.csv'); load('yoshimoto_5.csv')];

SSVEPdata = zeros(channels, max_t_length, sets, targetNum);
rawDataTemp = zeros(max_t_length,1);

for channel = 1:channels
    for set = 1:sets
        % セット数の条件を見たすデータのインデックスを列挙
        eachSetData = find(data(:,1)==(set-1));
        setData = data(eachSetData,:);
        for targetId = 1:targetNum
            eachTargetdata = find(setData(:,12)==(targetId-1));
            targetData = setData(eachTargetdata,:);
%             if (channel == 1) && (targetId == 1)
%                 figure;
%                 [numRows,~] = size(targetData);
%                 x=1:numRows;
%                 plot(x,targetData(:,4),x,targetData(:,5),x,targetData(:,6),x,targetData(:,7),x,targetData(:,8),x,targetData(:,9),x,targetData(:,10),x,targetData(:,11))
%                 legend("P3","P7","O1","Pz","Oz","O2","P8","P4");
%             end
            for sumpleID = 1:max_t_length
                rawDataTemp(sumpleID) = targetData((sumpleID+32),(channel+3));
%                 view = rawDataTemp(sumpleID)
            end
            %bandpassDataTemp = bandpass(rawDataTemp, [0.1 60], 128);
            %zRawDataTemp = zscore(bandpassDataTemp);
            zRawDataTemp = zscore(rawDataTemp);
%             m = mean(zRawDataTemp)
%             s = std(zRawDataTemp)
            for sumpleID = 1:max_t_length
                SSVEPdata(channel, sumpleID, set, targetId) = zRawDataTemp(sumpleID);
            end
        end
    end
end

a = SSVEPdata(1,:,1,1);

save("yoshiInew.mat", "SSVEPdata");