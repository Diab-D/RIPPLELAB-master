function  [ ] = plotResults( inputRaw, RippleLab_Data, S_struct, t, f, spect_list_S, spect_list_S97, spect_list_S90 )
    % Plot the results for the 
    %       - Raw iEEG data
    %       - Filtered iEEG data
    %       - T-f map
    %       - 0% denoised spectrum
    %       - 97% denoised spectrum
    
    %% Plot of raw data
    figure(1); plot(inputRaw);
    ylabel('Amplitude (mV)');
    
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 20)
    
    savefig('raw_data_HFOnr1_GRMES1');
    
    %% Plot of the filtered data
    figure(2); plot(RippleLab_Data.EEG_GRMES_1.v_Intervals{1,1});
    ylabel('Amplitude (mV)');
    
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 20)
    
    savefig('filtered_data_HFOnr1_GRMES1');
    %% T-f map
    % Plot the original t,f-spectrum against the 90% and 97% energy preserved
    figure(3);
    title('Amplitude spectrogram of the signal')

    subplot(1,3,1);
    surf(t, f, S_struct.S);
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 20)
    title('100% energy preserved')
    subplot(1,3,2);
    surf(t, f, S_struct.S_97);
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 20)
    title('97% energy preserved')
    subplot(1,3,3);
    surf(t, f, S_struct.S_90);
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 20)
    title('90% energy preserved')

    for ii = 1 : 3
        subplot(1,3,ii);
        colormap(jet);
        shading interp
        axis tight
        view(0, 90)
        set(gca, 'FontName', 'Times New Roman', 'FontSize', 20)
        if ii == 1
            ylabel('Frequency, [Hz]') 
        elseif ii == 2
            xlabel('Time, [s]')
        elseif ii == 3
            hcol = colorbar;
            set(hcol, 'FontName', 'Times New Roman', 'FontSize', 20)
            ylabel(hcol, 'Magnitude, [dB]')
        end
    end
    
    savefig('t_f_maps_HFOnr1_GRMES1');
    
    %% Denoising spectra
    figure(4);
    subplot(1,3,1);
    plot(f,spect_list_S);
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 20)
    title('100% energy preserved')
    subplot(1,3,2);
    plot(f,spect_list_S97);
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 20)
    title('97% energy preserved')
    subplot(1,3,3);
    plot(f,spect_list_S90);
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 20)
    title('90% energy preserved')

    for ii = 1 : 3
        subplot(1,3,ii);
        set(gca, 'FontName', 'Times New Roman', 'FontSize', 20)
        if ii == 1
            set(gca, 'FontName', 'Times New Roman', 'FontSize', 20)
            ylabel('Magnitude, [-]');
            set(gca, 'FontName', 'Times New Roman', 'FontSize', 20)
        elseif ii == 2
            xlabel('Frequency, [Hz]')
        end
    end
    
    savefig('denoised_spectra_HFOnr1_GRMES1');
    
    figure(5);
    plot(f,spect_list_S);
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 20)
    ylabel('Magnitude, [-]');
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 20)
    xlabel('Frequency, [Hz]')
    savefig('spectrum_100_HFOnr1_GRMES1')
    xlim([0 500]);
    
    figure(6);
    plot(f,spect_list_S97);
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 20)
    ylabel('Magnitude, [-]');
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 20)
    xlabel('Frequency, [Hz]')
    savefig('spectrum_97_HFOnr1_GRMES1');
    xlim([0 500]);
    
end