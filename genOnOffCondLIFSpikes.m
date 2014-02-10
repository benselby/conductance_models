function [OnOffSpikes OnOffCount] = genOnOffCondLIFSpikes(dt, ga, gb, g_params, neuron_params)
    
    reset_v = -0.07;
    v_thresh = -0.05;
    R = 0.1;
    
    if nargin < 5

        %Set default neuron parameters
        tau_ref = 0.002;
        tau_rc = 0.02;
        
        if nargin < 4
           g_bias = 0; %Default g_bias 
           g_gain = 1;
        else
            g_bias = g_params(1);
            g_gain = g_params(2);
        end
     
    else
       
        %Get parameters
        tau_ref = neuron_params(1);
        tau_rc = neuron_params(2);
        
    end
    
    %Default voltage 
    on_v = reset_v;
    off_v = reset_v;
    
    %Initialize outputs
    OnOffSpikes = zeros(2, length(ga));
    OnOffCount = zeros(2,1);
    
    i = 1;
    while i < length(ga)
               
        %Get the current
        current = getCurrent(ga(i), gb(i), on_v, [g_bias g_gain]);
                
        %Calc dv
        on_dv = (1/tau_rc)*(on_v+(current)*R) * dt;
        
        %Update voltage
        on_v = on_v + on_dv;
        
        if(on_v < reset_v)
           on_v = reset_v; 
        end
        
        %Create a spike
        if( on_v >= v_thresh)
            
            %Find first nonzero spike time
            index = find(OnOffSpikes(1,:) == 0);
            OnOffSpikes(1,index(1)) = i*dt;
            
            %Keep count of the spike
            OnOffCount(1) = OnOffCount(1) + 1;
            
           %Wait tau_ref before initiating a new spike
           on_v = reset_v;
           i = i + ceil(tau_ref/dt);
        end
        
        i = i+1;
    end    
    
    
    i = 1;
    while i < length(gb)
               
        %Get the current
        current = getCurrent(gb(i), ga(i), off_v, [g_bias g_gain]);
                
        %Calc dv
        off_dv = (1/tau_rc)*(off_v+(current)*R) * dt;
        
        %Update voltage
        off_v = off_v + off_dv;
        
        if(off_v < reset_v)
           off_v = reset_v; 
        end
        
        %Create a spike
        if( off_v >= v_thresh)
            
            %Find first nonzero spike time
            index = find(OnOffSpikes(2,:) == 0);
            OnOffSpikes(2,index(1)) = i*dt;
            
            %Keep count of the spike
            OnOffCount(2) = OnOffCount(2) + 1;
            
           %Wait tau_ref before initiating a new spike
           off_v = reset_v;
           i = i + ceil(tau_ref/dt);
        end
        
        i = i+1;
    end
    
end


