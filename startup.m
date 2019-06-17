
clear all; clc; close all;

addpath(genpath('C:\Users\s152891\OneDrive - TU Eindhoven\jaar 4, 2018-2019\Kwartiel 3\Klinische module, neurologie\Project\final scripts'));
addpath(genpath('C:\Users\s152891\Documents\MATLAB\RIPPLELAB-master'));
addpath(genpath('C:\Users\s152891\Documents\MATLAB\STFT'));

%% Initial detection

%open('p_RippleLab.m')

%% Load raw and Ripplelab data

% Ripplelab
raw_Data_1 = load('originele_data_deel1.mat');                                           %filename raw data
%Raw_Data_2 = load('originele_data_deel2.mat');                                          %filename raw data
%Raw_Data_3 = load('originele_data_deel3.mat');                                          %filename raw data
%Raw_Data_Total = [Raw_Data_1.record_deel1;Raw_Data_2.record_deel2;Raw_Data_3.record_deel3];

RippleLab_Data = importdata('MARIO07_200116_GRMES_1-8.rhfe');                            %filename Ripplelabdata

%% Time-frequency analysis

RippleLab_Data1_32 = importdata('1-32final.rhfe');
nbrElectrode = 32;
names = fieldnames(RippleLab_Data1_32);

for i = 2:33
    x{i-1} = cell(size(RippleLab_Data1_32.(names{i}).st_HFOInfo.m_EvtLims,1),1); %initialise sizes 
    for j =  1:size(x{i-1},1)
        event{j,:} = RippleLab_Data1_32.(names{i}).st_HFOInfo.m_EvtLims(j,1) : RippleLab_Data1_32.(names{i}).st_HFOInfo.m_EvtLims(j,2);
        %event(j:length(event) ) = [];
        %x{i-1} = event;
    end
    event(j+1:length(event) ) = [];
    x{i-1} = event;
end

[ S, f, t, X ] = time_frequency_analysis(RippleLab_Data, Unfiltered_HFO_Data, HFOnr);

%% Denoising

S_abs        = abs(S);                                          % S is the STFT matrix
S_abs        = S_abs.^2;
S_struct.(genvarname('S')) = S;

S_1 = S_abs;
S_2 = S_abs;

denoise_factor_1 = 0.90;                                        % Conservation of 90 percent of the energy
[ S_struct ] = denoising(S_1, S, S_struct, denoise_factor_1);

denoise_factor_2 = 0.97;                                        % Conservation of 97 percent of the energy
[ S_struct ] = denoising(S_2, S, S_struct, denoise_factor_2);

% Compute the spectrogram
spect_list_S = zeros(length(S_struct.S),1);
spect_list_S90 = zeros(length(S_struct.S),1);
spect_list_S97 = zeros(length(S_struct.S),1);

for ii = 1 : size(S_struct.S,1)
    tot_kol_S = sum(S_struct.S(ii,:));
    tot_kol_S90 = sum(S_struct.S_90(ii,:));
    tot_kol_S97 = sum(S_struct.S_97(ii,:));
    spect_list_S(ii) = tot_kol_S;
    spect_list_S90(ii) = tot_kol_S90;
    spect_list_S97(ii) = tot_kol_S97;
end

spect_list_S(spect_list_S < 0) = 0;
spect_list_S90(spect_list_S90 < 0) = 0;
spect_list_S97(spect_list_S97 < 0) = 0;

%% Plot the original t,f-spectrum against the 90% and 97% energy preserved
% figure(2);
% title('Amplitude spectrogram of the signal')
% 
% subplot(1,3,1);
% surf(t, f, S_struct.S);
% title('100% energy preserved')
% subplot(1,3,2);
% surf(t, f, S_struct.S_97);
% title('97% energy preserved')
% subplot(1,3,3);
% surf(t, f, S_struct.S_90);
% title('90% energy preserved')
% 
% for ii = 1 : 3
%     subplot(1,3,ii);
%     colormap(jet);
%     shading interp
%     axis tight
%     view(0, 90)
%     set(gca, 'FontName', 'Times New Roman', 'FontSize', 14)
%     xlabel('Time, s')
%     ylabel('Frequency, Hz') 
% 
%     hcol = colorbar;
%     set(hcol, 'FontName', 'Times New Roman', 'FontSize', 14)
%     ylabel(hcol, 'Magnitude, dB')
% 
% end
% 
% %% Plot the denoised spectra
% figure(3);
% subplot(1,3,1);
% plot(f,spect_list_S);
% title('100% energy preserved')
% subplot(1,3,2);
% plot(f,spect_list_S97);
% title('97% energy preserved')
% subplot(1,3,3);
% plot(f,spect_list_S90);
% title('90% energy preserved')
% 
% for ii = 1 : 3
%     subplot(1,3,ii);
%     set(gca, 'FontName', 'Times New Roman', 'FontSize', 14)
%     xlabel('Frequency, Hz')
%     ylabel('Magnitude, -') 
% end

%% Plot the results in one figure
plotResults(Unfiltered_HFO_Data, RippleLab_Data, S_struct, t, f, spect_list_S, spect_list_S97, spect_list_S90);

%% Feature extraction
% For further research
%[ fm ] = peakToNotch(S);
%[ ] = spectralEntropy(S);
%[ Sub_band_power_ratio ] = subBandPowerRatio(S);
