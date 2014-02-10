function [onOffCount, tOnOff] = genOnOffLIFSpikes(x, dt, params)
    
    if nargin < 3

        %Set neuron parameters
        tau_ref = 0.002;
        tau_rc = 0.02;
        R = 1;
        v_thresh = 1;
        
        %Make the on neuron
        aon_x0 = 40;
        aon_x1 = 150;
        j_bias_on = (1-exp( (1-aon_x0*tau_ref)/((-1)*aon_x0*tau_rc)  ))^(-1);
        alpha_on = (1-exp( (1-aon_x1*tau_ref)/((-1)*aon_x1*tau_rc)  ))^(-1) - j_bias_on;
    
    else
       
        %Get parameters
        tau_ref = params(1);
        tau_rc = params(2);
        v_thresh = params(3);
        j_bias_on = params(4);
        alpha_on = params(5);
        R = 1;
        
    end
    
    %Initialize outputs
    tOnOff = zeros(2, 2*floor(dt*length(x)/tau_ref));
    onOffCount = zeros(1,2);
    
    %Make the off neuron
    j_bias_off = j_bias_on; 
    alpha_off = (-1) * alpha_on;
    
    %Initialize spike trains
    on_spikes = [];
    off_spikes = [];
    
    %Initialize voltages
    on_v = 0;
    off_v = 0;
    
    %%% Get spikes for on neuron %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    i = 1;
    while i < length(x)
       
        %Calc dv
        on_dv = (-1/tau_rc)*(on_v-(alpha_on*x(i) + j_bias_on)*R) * dt;
        
        %Update voltage
        on_v = on_v + on_dv;
        
        %Create a spike
        if( on_v >= v_thresh)
           on_spikes = [on_spikes i*dt]; 
           on_v = 0;
           
           %Wait tau_ref before initiating a new spike
           i = i + ceil(tau_ref/dt);
        end
        
        i = i+1;
    end
    
    %%% Get spikes for off neuron %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    i = 1;
    while i < length(x)
       
        %Calc dv
        off_dv = (-1/tau_rc)*(off_v-(alpha_off*x(i) + j_bias_off)*R) * dt;
        
        %Update voltage
        off_v = off_v + off_dv;
        
        %Create a spike
        if( off_v >= v_thresh)
           off_spikes = [off_spikes i*dt]; 
           off_v = 0;
           
           %Wait tau_ref before initiating a new spike
           i = i + ceil(tau_ref/dt);
        end
        
        i = i+1;
    end
    
    %Get the number of spikes
    onOffCount(1) = length(on_spikes);
    onOffCount(2) = length(off_spikes);
    
    %Store the spikes
    tOnOff(1, 1:onOffCount(1)) = on_spikes;
    tOnOff(2, 1:onOffCount(2)) = off_spikes;
    
end