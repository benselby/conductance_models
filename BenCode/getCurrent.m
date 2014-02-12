function J = getCurrent(g_excitatory, g_inhibitory, membrane_v, neuron_parms)

    if nargin < 4
        
       g_bias = 0;
       g_gain = 1;
       
    else
        
        g_bias = neuron_parms(1);
        g_gain = neuron_parms(2);
        
    end

    %Nernst potential
    nernst_excititory = 0.005;
    nernst_inhibitory = -0.07;
    
    J = (g_gain*g_excitatory).*(nernst_excititory - membrane_v) - (g_inhibitory + g_bias).*(membrane_v - nernst_inhibitory);

end