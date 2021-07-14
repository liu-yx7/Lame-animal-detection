clear all
clc

q = 10; %quantity of the sample

for c1 = 1:q

%% Ancortek reading part
[filename,pathname] = uigetfile('*.dat');
fileID = fopen(filename, 'r');
dataArray = textscan(fileID, '%f');
fclose(fileID);
radarData = dataArray{1};
clearvars fileID dataArray ans;
fc = radarData(1); % Center frequency
Tsweep = radarData(2); % Sweep time in ms
Tsweep=Tsweep/1000; %then in sec
NTS = radarData(3); % Number of time samples per sweep
Bw = radarData(4); % FMCW Bandwidth. For FSK, it is frequency step;
% For CW, it is 0.
Data = radarData(5:end); % raw data in I+j*Q format

fs=NTS/Tsweep; % sampling frequency ADC
record_length=length(Data)/NTS*Tsweep; % length of recording in s
nc=record_length/Tsweep; % number of chirps

%% Reshape data into chirps and do range FFT (1st FFT)
Data_time=reshape(Data, [NTS nc]);

% Data_time=Data_time(:,20001:end);

%Hamming window prior to FFT may reduce the sidelobes in range but not
%relevant for Doppler (will also need to discard 2 raneg bins not 1 later
%on!!!)
% win = repmat(hamming(NTS),1,size(Data_time,2));
win = ones(NTS,size(Data_time,2));

%Part taken from Ancortek code for FFT and IIR filtering
tmp = fftshift(fft(Data_time.*win),1);
Data_range(1:NTS/2,:) = tmp(NTS/2+1:NTS,:);
% IIR Notch filter
ns = oddnumber(size(Data_range,2))-1;
Data_range_MTI = zeros(size(Data_range,1),ns);
[b,a] = butter(4, 0.0075, 'high');
[h, f1] = freqz(b, a, ns);
for k=1:size(Data_range,1)
  Data_range_MTI(k,1:ns) = filter(b,a,Data_range(k,1:ns));
end
freq =(0:ns-1)*fs/(2*ns); 
range_axis=(freq*3e8*Tsweep)/(2*Bw);
%Need to remove the first range bin as it has got some strong residual from
%filtering?
Data_range_MTI=Data_range_MTI(2:size(Data_range_MTI,1),:);
Data_range=Data_range(2:size(Data_range,1),:);

% figure
% colormap(jet)
% % imagesc([1:10000],range_axis,20*log10(abs(Data_range_MTI)))
% imagesc(20*log10(abs(Data_range_MTI)))
% xlabel('No. of Sweeps')
% ylabel('Range bins')
% title('Range Profiles after MTI filter')
% clim = get(gca,'CLim'); axis xy; ylim([1 100])
% set(gca, 'CLim', clim(2)+[-60,0]);
% drawnow

%% Spectrogram processing for 2nd FFT to get Doppler
% This selects the range bins where we want to calculate the spectrogram
bin_indl = 3;
bin_indu = 60;
%Parameters for spectrograms
MD.PRF=1/Tsweep;
MD.TimeWindowLength = 600;
MD.OverlapFactor = 0.5;
MD.OverlapLength = round(MD.TimeWindowLength*MD.OverlapFactor);
MD.Pad_Factor = 4;
MD.FFTPoints = MD.Pad_Factor*MD.TimeWindowLength;
MD.DopplerBin=MD.PRF/(MD.FFTPoints);
MD.DopplerAxis=-MD.PRF/2:MD.DopplerBin:MD.PRF/2-MD.DopplerBin;
MD.WholeDuration=size(Data_range_MTI,2)/MD.PRF;
MD.NumSegments=floor((size(Data_range_MTI,2)-MD.TimeWindowLength)/floor(MD.TimeWindowLength*(1-MD.OverlapFactor)));
    
%Method 1 - COHERENT SUM
% myvec_MTI=sum(Data_range_MTI(bin_indl:bin_indu,:));
% myvec=sum(Data_range(bin_indl:bin_indu,:));
% Data_spec_MTI=fftshift(spectrogram(myvec_MTI,MD.TimeWindowLength,MD.OverlapLength,MD.FFTPoints),1);
% Data_spec=fftshift(spectrogram(myvec,MD.TimeWindowLength,MD.OverlapLength,MD.FFTPoints),1);

%Method 2 - SUM OF RANGE BINS
Data_spec_MTI2=0;
Data_spec2=0;
for RBin=bin_indl:1:bin_indu
    Data_MTI_temp = fftshift(spectrogram(Data_range_MTI(RBin,:),MD.TimeWindowLength,MD.OverlapLength,MD.FFTPoints),1);
    Data_spec_MTI2=Data_spec_MTI2+abs(Data_MTI_temp);                                
    Data_temp = fftshift(spectrogram(Data_range(RBin,:),MD.TimeWindowLength,MD.OverlapLength,MD.FFTPoints),1);
    Data_spec2=Data_spec2+abs(Data_temp);
end

%Method 3 - AVERAGE OF RANGE BINS
% myvec_MTI3=mean(Data_range_MTI(bin_indl:bin_indu,:));
% myvec3=mean(Data_range(bin_indl:bin_indu,:));
% Data_spec_MTI3=fftshift(spectrogram(myvec_MTI3,MD.TimeWindowLength,MD.OverlapLength,MD.FFTPoints),1);
% Data_spec3=fftshift(spectrogram(myvec3,MD.TimeWindowLength,MD.OverlapLength,MD.FFTPoints),1);

MD.TimeAxis=linspace(0,MD.WholeDuration,size(Data_spec_MTI2,2));

% Normalise and plot micro-Doppler
% Data_spec=flipud(Data_spec./max(Data_spec(:)));
% Data_spec_MTI=flipud(Data_spec_MTI./max(Data_spec_MTI(:)));
% Data_spec2=flipud(Data_spec2./max(Data_spec2(:)));
% Data_spec_MTI2=flipud(Data_spec_MTI2./max(Data_spec_MTI2(:)));
Data_spec_MTI2=flipud(Data_spec_MTI2);
% Data_spec3=flipud(Data_spec3./max(Data_spec3(:)));
% Data_spec_MTI3=flipud(Data_spec_MTI3./max(Data_spec_MTI3(:)));

%% time segment 
% warmup (2 min), walk(3 min) walk fast( 3 min)  and cool down(2 min)
% cen_mean = zeros(3,1);
% bdw_mean = zeros(3,1);
% cen_std = zeros(3,1);
% bdw_std = zeros(3,1);
feature = zeros(3,4);

for c2 = 1:3
    col_num = size(Data_spec_MTI2,2);
    switch c2
        case 1
%             cut_between = [1,round(col_num*0.2)];
%             Data_spec_MTI2=Data_spec_MTI2(:,1:round(col_num*0.2));
%             MD.TimeAxis=MD.TimeAxis(1:round(col_num*0.2));
            Data_spec_MTI2_div=Data_spec_MTI2(:,1:round(col_num/3));
%             MD.TimeAxis=MD.TimeAxis(1:round(col_num/3));

        case  2
%             cut_between = [round(col_num*0.2),round(col_num*0.5)];
            Data_spec_MTI2_div=Data_spec_MTI2(:,round(col_num/3):round(col_num*2/3));
%             MD.TimeAxis_mod = MD.TimeAxis(round(col_num/3):round(col_num*2/3));

        case 3
%             cut_between = [round(col_num*0.5),round(col_num*0.8)];
            Data_spec_MTI2_div=Data_spec_MTI2(:,round(col_num*2/3):col_num);
%             MD.TimeAxis_mod=MD.TimeAxis(round(col_num*2/3):col_num);
        case 4
%             cut_between = [round(col_num*0.8),col_num];
            Data_spec_MTI2_div=Data_spec_MTI2(:,round(col_num*0.8),col_num);
%             MD.TimeAxis_mod=MD.TimeAxis(round(col_num*0.8),col_num);
    end

        
%% Cut from 1 to 9 s
% Data_spec_MTI2=Data_spec_MTI2(:,100:895);
% MD.TimeAxis=MD.TimeAxis(100:895);

%%
% centriod
column_add = sum(Data_spec_MTI2_div,1);
[m,n] = size(Data_spec_MTI2_div);

centriod = zeros(1,n);
for i = 1:n
    centriod(1,i) = (1:m)*Data_spec_MTI2_div(:,i)/column_add(1,i);
end

%bandwidth
bandwidth = zeros(1,n);
for j = 1:n
    bandwidth(1,j) = ((((1:m)- repmat(centriod(1,j),1,m)).^2)* Data_spec_MTI2_div(:,j)/column_add(1,j)).^0.5;
end

%%

% data  extraction (static Mean & Standard deviation)

% cen_mean(c2,1) = mean(centriod,2);
% bdw_mean(c2,1) = mean(bandwidth,2);
% cen_std(c2,1)= std(centriod,0,2);
% bdw_std(c2,1) = std(bandwidth,0,2);

feature(c2,1) = mean(centriod,2);
feature(c2,2) = mean(bandwidth,2);
feature(c2,3)= std(centriod,0,2);
feature(c2,4) = std(bandwidth,0,2);

% normalization
Mapped_feature = mapminmax(feature, -1, 1);
end
% Writing data into excel
timeData={datestr(now,30)};
if~exist('myData1.xlsx','file')
xlswrite('myData1.xlsx',Mapped_feature,1,'A1');
% xlswrite('myData1.xlsx',bdw_mean,1,'B1');
% xlswrite('myData1.xlsx',cen_std,1,'C1');
% xlswrite('myData1.xlsx',bdw_std,1,'D1');
% xlswrite('myData1.xlsx',timeData,1,'E1');
else

 [tmp1,tmp2,tmpRaw]=xlsread('myData1.xlsx');

if size(tmp1,1)==0&&size(tmp2,1)==0

mRowRange='1';

else

mRowRange=num2str(size(tmpRaw,1)+1);

end


xlswrite('myData1.xlsx',Mapped_feature,1,['A',mRowRange]);
% xlswrite('myData1.xlsx',bdw_mean,1,['B',mRowRange]);
% xlswrite('myData1.xlsx',cen_std,1,['C',mRowRange]);
% xlswrite('myData1.xlsx',bdw_std,1,['D',mRowRange]);
% xlswrite('myData1.xlsx',timeData,1,['E',mRowRange]);
end
%%
%Initialization variable
clear Data_range 

end
% draw table 
% f ={'Non Lame';'Non Lame';'Non Lame';'Non Lame';'Non Lame';'Non Lame';'Lame';'Lame';'Lame';'Lame';'Lame';'Lame';'Lame';'Lame';'Non Lame';'Non Lame';'Non Lame';'Non Lame';'Non Lame';};
% 
% Lable = char(f);
% feature = table(cen_mean,bdw_mean,cen_std,bdw_std,Lable);
% 
% filename = 'Warwick feature_win600overlap0.5.txt';
% writetable(feature,filename);
%%
% figure
% imagesc(MD.TimeAxis,MD.DopplerAxis,20*log10(abs(Data_spec)),[-50 0]); colormap('jet'); axis xy
% figure
% imagesc(MD.TimeAxis,MD.DopplerAxis,20*log10(abs(Data_spec_MTI))); colormap('jet'); axis xy
% ylim([-250 250])
% % figure
% imagesc(MD.TimeAxis,MD.DopplerAxis,20*log10(abs(Data_spec2)),[-50 0]); colormap('jet'); axis xy
% ylim([-250 250])


% figure
% imagesc(MD.TimeAxis,MD.DopplerAxis.*3e8/2/5.8e9,20*log10(abs(Data_spec_MTI2))); colormap('jet'); axis xy
% ylim([-6 6]); colorbar
% colormap; %xlim([1 9])
% clim = get(gca,'CLim');
% set(gca, 'CLim', clim(2)+[-40,0]);
% xlabel('Time[s]', 'FontSize',16);
% ylabel('Velocity [m/s]','FontSize',16)
% set(gca, 'FontSize',16)
% title(filename)
% 
% figure
% imagesc(MD.TimeAxis,MD.DopplerAxis,20*log10(abs(Data_spec_MTI2))); colormap('jet'); axis xy
% ylim([-500 500]); colorbar;
% 
% figure
% plot(centriod); 
% ylim([1 800]);
% 
% figure
% plot(bandwidth); 





% figure
% imagesc(MD.TimeAxis,MD.DopplerAxis,20*log10(abs(Data_spec_MTI2))); colormap('jet'); axis xy
% ylim([-250 250]); xlabel('Time [s]'), ylabel('Doppler [Hz]'), colorbar

% figure
% imagesc(MD.TimeAxis,MD.DopplerAxis,20*log10(abs(Data_spec3)),[-50 0]); colormap('jet'); axis xy
% figure
% imagesc(MD.TimeAxis,MD.DopplerAxis,20*log10(abs(Data_spec_MTI3)),[-50 0]); colormap('jet'); axis xy
% ylim([-250 250])

