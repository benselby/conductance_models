function [spikes v cur] = genLIFSpikes(x, dt, ga, gb, g_params, neuron_params)
    
    reset_v = -0.07;
    v_thresh = -0.05;
    R = 0.1;%8000000;
    

    if nargin < 6

        %Set default neuron parameters
        tau_ref = 0.002;
        tau_rc = 0.02;
        
        if nargin < 5
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
    on_spikes = [];
    v = zeros(1,length(x));
    cur = zeros(1,length(x));
    
    i = 1;
    while i < length(x)
       
        v(i) = on_v;
        
        %Get the current
        current = getCurrent(ga, gb, on_v, [g_bias g_gain]);
        
        cur(i) = current;
        
        %Calc dv
        on_dv = (1/tau_rc)*(on_v+(current)*R) * dt;
        
        %Update voltage
        on_v = on_v + on_dv;
        
        if(on_v < reset_v)
           on_v = reset_v; 
        end
        
        %Create a spike
        if( on_v >= v_thresh)
           on_spikes = [on_spikes i*dt]; 
           on_v = reset_v;
           
           %Wait tau_ref before initiating a new spike
           i = i + ceil(tau_ref/dt);
        end
        
        i = i+1;
    end
    
    spikes = on_spikes;
    
    %Copy final values
    cur(length(x)) = current;
    v(length(x)) = on_v;

end

