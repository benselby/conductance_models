clear all;
close all;

dt = 0.0001;
max_t = 1;
t = [0:dt:max_t];

dga = 1;
max_ga = 100;
ga = [0:dga:max_ga];

dgb = 1;
max_gb = 100;
gb = [0:dga:max_ga];

%Create a 1 second signal
x = ones(1,length(t));

%Create a random g_bias
g_bias = 0;%rand(1) * 100;
g_gain = 5;%rand(1); %Lower gain means weaker excitation

%Plot a set of spikes, voltages, and currents
[spikes v curr] = genLIFSpikes(x, dt, 100, 0, [g_bias g_gain]);
figure(1);
hold on;
spike_plot = plot([spikes;spikes],[2*ones(size(spikes));2.5*ones(size(spikes))], 'k-');
plot(t, v, 'b-', t, curr, 'r-');
spike_group = hggroup;
set(spike_plot,'Parent',spike_group);
set(get(get(spike_group,'Annotation'),'LegendInformation'),'IconDisplayStyle','on');
legend('voltage', 'current', 'spikes');
xlabel('time (s)');

%Make a 3d plot of the firing rate as conductances change
a = zeros(length(ga), length(gb));
for(i = 1:length(ga))
    
    for(j = 1:length(gb))
        
        spikes = genLIFSpikes(x, dt, ga(i), gb(j), [g_bias g_gain]);
        
        spike_rate = length(spikes)/max_t;
        
        a(i,j) = spike_rate;
        
    end
     
end

%Plot the plot
X = ga;
Y = gb;
[X,Y] = meshgrid(X,Y);

figure(2);
surf(X,Y,a');
xlabel('Excititory Conductance');
ylabel('Inhibitory Conductance');
zlabel('Firing Rate (Hz)');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% POPULATION DECODING %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = 50;
dt = 0.001;
T = 1;

%Neuron parms
tau_rc = 0.02;
tau_ref = 0.002;
Jth = 1;
sigma = 0.1;

%Create a post-synaptic conductance (psg)
RandomSeed = 10;
tau_psg = 0.005;
n = 0;
range = 0:dt:0.15;
psg_kernel = exp(-range/tau_psg);
psg_kernel = .1*psg_kernel/sum(dt*psg_kernel); % normalize filter

% Generate input signal
upperBandLimit = 2*pi*5; % Upper frequency limit
lowerBandLimit = 2*pi*0;  % Lower frequency limit
rms = 1;           % RMS signal level
bandwidth = [lowerBandLimit upperBandLimit];
[S,Amps] = genSignal(T,dt,rms,bandwidth,RandomSeed);

% Standard LIF population for encoding arbitrary signal
x_ints = 4*rand(1,N) - 2;
rates = 100*rand(1,N)+200;
alpha = zeros(size(x_ints));

for i=1:N
    U = -1*(1/rates(i) - tau_ref)/tau_rc;
    alpha(i) = ((1 - exp(U))^-1 - 1)/(2-x_ints(i));
    J_bias(i) = 1 - alpha(i)*x_ints(i);
end

Non = 0; % 'on' neuron spike count
Noff = 0; % 'off' neuron spike count
TonOff = [0;0];

for i=1:N
    LIFparms = [tau_ref,tau_rc,Jth,J_bias(i),alpha(i)];    
    % generate spikes for 2 neurons
    [onOffCnt spike_times] = genOnOffLIFSpikes(S,dt,LIFparms);      
    TonOff = [TonOff spike_times];
    Non=Non+onOffCnt(1);     
    Noff = Noff+onOffCnt(2);  
end

on_times = TonOff(1, find(TonOff(1,:)));
off_times = TonOff(2, find(TonOff(2,:)));

on_spike_times = floor(on_times/dt);
off_spike_times = floor(off_times/dt);
on_spikes = zeros(1, T/dt);
off_spikes = zeros(1, T/dt);
on_spikes(on_spike_times) = 1;
off_spikes(off_spike_times) = 1;

%Get the inhibitory and excitatory conductances:
excitatory_g = conv(psg_kernel, on_spikes, 'full'); 
inhibitory_g = conv(psg_kernel, off_spikes, 'full'); 

%Clip the conductances (they're too long due to full convolution)
excitatory_g = excitatory_g(1:T/dt);
inhibitory_g = inhibitory_g(1:T/dt);

%%% Create 25 post synaptic Neurons %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
g_bias = -100 + rand(1,N/2) * 200;
g_gain = rand(1, N/2);

%Copy the conductances
ga = excitatory_g;
gb = inhibitory_g;

spikes = zeros(N, length(ga));
spikeCount = zeros(N,1);

for i=1:N/2
   
    %Get the spikes
    [s c] = genOnOffCondLIFSpikes(dt, ga, gb, [g_bias(i) g_gain(i)]);
    
    %Save spikes
    spikes(i, :) = spikes(i, :) + s(1,:);
    spikes(N/2+i, :) = spikes(N/2+i, :) + s(2,:);
    
    %Save spike count
    spikeCount(i) = c(1);
    spikeCount(N/2+i) = c(2);
end

%Convert spikes to firing rate for each neuron
bin = 0.01; %10ms bins
firing_rates = zeros(N, length(ga));
for i=1:T/bin
   
    for j=1:N
        
        %Find spikes
        index1 = find( spikes(j,:) > (i-1)*bin);
        index2 = find( spikes(j,index1) < (i)*bin);
       
        %Put in firing rates
        start_of_bin = (i-1)*bin/dt + 1;
        firing_rates(j, start_of_bin:start_of_bin+bin/dt) = length(index2)/bin;
       
    end
    
end


