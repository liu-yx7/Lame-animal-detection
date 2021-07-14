clear all
clc
% close all
%% Read the data from the file you select - the script must be in the same folder as the data
fileFolder=fullfile('H:\新建文件夹\data');
dirOutput=dir(fullfile(fileFolder,'*.dat'));
fileNames={dirOutput.name};
for fnum=1:length(fileNames)
filename = char(fileNames(fnum));
myname=filename(1:end-4)
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
record_length=size(Data,1)/fs; % length of recording in s
nc=record_length/Tsweep; % number of chirps

%% Reshape data into chirps and do 1st FFT to get Range-Time
%Note that also an IIR notch filter is applied to the data - this is done
%to remove the 0Hz components in the frequency domain that would be
%associated to static (0Hz Doppler shift) targets
Data_time=reshape(Data, [NTS nc]);
win = ones(NTS,size(Data_time,2));
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

Data_range_MTI=Data_range_MTI(2:size(Data_range_MTI,1),:);
Data_range=Data_range(2:size(Data_range,1),:);

%Plot of the Range-Time-Intensity - note that the plot is normalised at 40 dB and
%in logarithmic scale
figure
colormap(jet)
imagesc(20*log10(abs(Data_range_MTI)))
xlabel('No. of Sweeps')
ylabel('Range bins')
title('Range Profiles after MTI filter')
clim = get(gca,'CLim'); axis xy;
set(gca, 'CLim', clim(2)+[-40,0]);
drawnow
colorbar
im1=getframe;
Myname1 = strcat(myname,'RCS.png');
imwrite(im1.cdata,Myname1);
%% Spectrogram processing - This is done through STFT Short Time Fourier Transform
%MATLAB implements STFT with the spectrogram() function

% This selects the range bins where we want to calculate the spectrogram,
% basically where you see the target signature in the RTI
bin_indl = 5
bin_indu = 25
%Parameters for spectrograms
MD.PRF=1/Tsweep;
MD.TimeWindowLength = 200;
MD.OverlapFactor = 0.95;
MD.OverlapLength = round(MD.TimeWindowLength*MD.OverlapFactor);
MD.Pad_Factor = 4;
MD.FFTPoints = MD.Pad_Factor*MD.TimeWindowLength;
MD.DopplerBin=MD.PRF/(MD.FFTPoints);
MD.DopplerAxis=-MD.PRF/2:MD.DopplerBin:MD.PRF/2-MD.DopplerBin;
MD.WholeDuration=size(Data_range_MTI,2)/MD.PRF;
MD.NumSegments=floor((size(Data_range_MTI,2)-MD.TimeWindowLength)/floor(MD.TimeWindowLength*(1-MD.OverlapFactor)));

%There are three methods listed here to generate the spectrograms - we
%normally use the second one but it is not set in stone
%Method 1 - COHERENT SUM
myvec_MTI=sum(Data_range_MTI(bin_indl:bin_indu,:));
myvec=sum(Data_range(bin_indl:bin_indu,:));
Data_spec_MTI=fftshift(spectrogram(myvec_MTI,MD.TimeWindowLength,MD.OverlapLength,MD.FFTPoints),1);
Data_spec=fftshift(spectrogram(myvec,MD.TimeWindowLength,MD.OverlapLength,MD.FFTPoints),1);

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
myvec_MTI3=mean(Data_range_MTI(bin_indl:bin_indu,:));
myvec3=mean(Data_range(bin_indl:bin_indu,:));
Data_spec_MTI3=fftshift(spectrogram(myvec_MTI3,MD.TimeWindowLength,MD.OverlapLength,MD.FFTPoints),1);
Data_spec3=fftshift(spectrogram(myvec3,MD.TimeWindowLength,MD.OverlapLength,MD.FFTPoints),1);

MD.TimeAxis=linspace(0,MD.WholeDuration,size(Data_spec_MTI,2));

% Normalise and plot micro-Doppler
Data_spec=flipud(Data_spec./max(Data_spec(:)));
Data_spec_MTI=flipud(Data_spec_MTI./max(Data_spec_MTI(:)));
Data_spec2=flipud(Data_spec2./max(Data_spec2(:)));
Data_spec_MTI2=flipud(Data_spec_MTI2./max(Data_spec_MTI2(:)));
Data_spec3=flipud(Data_spec3./max(Data_spec3(:)));
Data_spec_MTI3=flipud(Data_spec_MTI3./max(Data_spec_MTI3(:)));

%% Part to plot the spectrograms
% figure
% imagesc(MD.TimeAxis,MD.DopplerAxis,20*log10(abs(Data_spec)),[-50 0]); colormap('jet'); axis xy
% figure
% imagesc(MD.TimeAxis,MD.DopplerAxis,20*log10(abs(Data_spec_MTI)),[-50 0]); colormap('jet'); axis xy
% figure
% imagesc(MD.TimeAxis,MD.DopplerAxis,20*log10(abs(Data_spec2)),[-50 0]); colormap('jet'); axis xy
% ylim([-250 250])
figure
imagesc(MD.TimeAxis,MD.DopplerAxis,20*log10(abs(Data_spec_MTI2)),[-45 0]); colormap('jet'); axis xy
ylim([-150 150]); xlabel('Time [s]'), ylabel('Doppler [Hz]'), colorbar
im2=getframe;
Myname = strcat(myname,'FT.png');
imwrite(im2.cdata,Myname);
close all;
clearvars -except fnum dirOutput fileNames fileFolder
end
% saveas(gcf,strcat(myname,'.png'));
% figure
% imagesc(MD.TimeAxis,MD.DopplerAxis,20*log10(abs(Data_spec3)),[-50 0]); colormap('jet'); axis xy
% figure
% imagesc(MD.TimeAxis,MD.DopplerAxis,20*log10(abs(Data_spec_MTI3)),[-50 0]); colormap('jet'); axis xy

