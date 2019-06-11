
function [ S, f, t, X ] = time_frequency_analysis(RippleLab, unfiltered, HFOnr)
%%% Input:
    % A: imported EEG data from the indicated electrodes.
    % HFOnr: number of the initially detected HFO.

%%% Output:
    % S: STFT-matrix (only unique points, time across columns, frequency across rows)
    % f: frequency vector, Hz
    % t: time vector, s
    % C: coherent amplification of the window
    % x: imported EEG data from the HFOnr detected on electrode GRMES2

% Intervals and sample rate 
fs      = RippleLab.EEG_GRMES_5.st_HFOInfo.s_Sampling;    % Sample rate 

%maxInt = A.EEG_GRMES_1.st_HFOInfo.m_IntervLims(:,2);
%minInt = A.EEG_GRMES_1.st_HFOInfo.m_IntervLims(:,1);
%midden_van_HFO = (maxInt + minInt)/2;
%x = (maxInt + minInt)/2;

% Define analysis parameters
wlen = 256;                                             % window length (recomended to be power of 2)
hop  = wlen/4;                                          % hop size (recomended to be power of 2)
nfft = 4096;                                            % number of fft points (recomended to be power of 2)

% Perform STFT
win = hanning(wlen, 'periodic');
%win2 = blackman(wlen, 'periodic');
%unfiltered = unfiltered - mean(unfiltered);
unfiltered = unfiltered - unfiltered(1);
[S, f, t, X] = stft(unfiltered, win, hop, nfft, fs);

% calculate the coherent amplification of the window
C = sum(win)/wlen;

% take the amplitude of fft(x) and scale it, so not to be a
% function of the length of the window and its coherent amplification
S = abs(S)/wlen/C;

% correction of the DC & Nyquist component
if rem(nfft, 2)                     % odd nfft excludes Nyquist point
    S(2:end, :) = S(2:end, :).*2;
else                                % even nfft includes Nyquist point
    S(2:end-1, :) = S(2:end-1, :).*2;
end

% convert amplitude spectrum to dB (min = -120 dB)
S = 20*log10(S + 1e-6);

% plot the spectrogram
% figure(1)
% surf(t, f, S)
% shading interp
% axis tight
% view(0, 90)
% set(gca, 'FontName', 'Times New Roman', 'FontSize', 14)
% xlabel('Time, s')
% ylabel('Frequency, Hz')
% title('Amplitude spectrogram of the signal')
% 
% hcol = colorbar;
% set(hcol, 'FontName', 'Times New Roman', 'FontSize', 14)
% ylabel(hcol, 'Magnitude, dB')