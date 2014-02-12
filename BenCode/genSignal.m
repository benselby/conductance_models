function [x,A] = genSignal(T,dt,rms,bandwidth,randomSeed,type)

    if nargin < 6
        type = 0;
    end

    %Seed future random number generations
    randn('seed',randomSeed);

    %Get the number of time steps
    num_time_steps = 2*floor(T/(2*dt)) + 1;
    num_osc_freq = (num_time_steps-1)/2; %Because of symmetry
    
    %Get angular frequencies for each oscillation frequency step
    w = 2*pi/T * [0:1:num_osc_freq];

    %Initialize amplitudes for each angular freq
    A = zeros(1, length(w));
    
    %Generate random amplitudes for all values
    A = randn(1, length(A));
    
    %Do clipping based on requested power spectrum
    if(type == 1)
        
        %Get center of bandwidth
        bw_c = bandwidth(1) + (bandwidth(2) - bandwidth(1))/2;
        
        %Get variance
        bw_var = (bandwidth(2) - bandwidth(1))/2;
        
        %Create a gaussian centered around the bandwith
        gauss = gaussmf(w, [bw_var bw_c]);
        
        %Multiply the amplitudes by the gaussian
        A = A .* gauss;
        
    else
        
        %Multiply everything in A that isn't the bandwidth by 0
        A = A .* (abs(w) >= bandwidth(1));
        A = A .* (abs(w) <= bandwidth(2));
           
    end
    
    
    %Flip these amplitudes because of the symetry (don't repeat the zero)
    flipped_A = fliplr(conj(A(2:end)));
    
    %Add flipped values to original
    A = [A, flipped_A];
    
    %Do an inverse FFT to get the signal
    x = ifft(A);
    
    %Remove unreal components
    x = real(x);
    
    %Get the signals RMS
    signal_rms = sqrt( sum(x.^2)/length(x) );
    
    %Apply the requested RMS
    x = x * (rms/signal_rms);
    A = A * (rms/signal_rms);

end

