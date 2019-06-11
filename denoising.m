function [ S_struct, spect_list_S, spect_list_S90, spect_list_S97 ] = denoising(S, S_original, S_struct, denoise_factor)
    % 1. Compute cumulative energy of the time-freq map using expansion
    % coefficients for each event  --> POWER SPECTRUM
    % 2. Retain only the components that made up the majority of the total energy
    % and omit the rest
    % Signalen die niet significant bijdragen aan de cumulative energy uit het 
    % power spectrum filteren.

    %B = cumsum(S);
    
    %%  Sort all values in the time-frequency spectrum
    %S        = abs(S);      % Alle waarden in 1 kolom zetten
    S        = S.^2; 
    S1       = S(:);
    S1_sorted = sort(S1, 'descend');
    
    %% Define threshold
    threshold = sum(S1) * denoise_factor;           % Dit nog uitbreiden naar ook 97%
    
    %% Define what values have to be deleted from S
    tot = 0;
    thres_val = 0;
    
    for ii = 1 : length(S1_sorted)      % Loopen door alle waardes die in S1 staan
        if tot >= threshold             % Als het optellen groter of gelijk is aan de threshold waarde,
            thres_val = S1_sorted(ii);  % Definieer de waarde die moet worden verwijderd uit S en de 
                                        % waardes die kleiner zijn dan thres_val
            break                       % Verlaat de for-loop
        end
        tot = tot + S1_sorted(ii);
    end  
    
    S(S <= thres_val) = 0;
   
    for jj = 1 : size(S_original,2)
        for ii = 1 : length(S_original)
            if S(ii,jj) == 0
                S_original(ii,jj) = 0; % min(min(S_original));
            end
        end
    end
    S_struct.(genvarname(['S_', num2str(denoise_factor*100)])) = S_original;
    
end


