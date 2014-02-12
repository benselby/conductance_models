% Dynamics of an integrate-and-fire conductance model

inLevels = 0:5:400;
exLevels = 0:5:100;

in = repmat(inLevels', 1, length(exLevels));
ex = repmat(exLevels, length(inLevels), 1);

dt = .0001;
T = 2;
time = dt:dt:T;

Vres = -65;
Vth = -50;
EEx = 0;
EIn = -65;
Eleak = -65;
gMaxEx = 1;
gMaxIn = 1;
gLeak = 0;
Tref = .002;
C = 1; % we set this to 1 and omit it from the loop 

%TODO: add leak -- is this an OK bias term?

V = zeros(size(in));
spikeCounts = zeros(size(in));
lastSpikeTime = -ones(size(in));

for i = 1:length(time)
    dVdt = -gMaxEx*ex.*(V-EEx) - gMaxIn*in.*(V-EIn) - gLeak*(V-Eleak);
    refractory = find(time(i)-lastSpikeTime < Tref);
    dVdt(refractory) = 0;
    V = V + dt*dVdt;
    spikeInd = find(V>=Vth);
    lastSpikeTime(spikeInd) = time(i);% - dVdt(spikeInd);
    V(spikeInd) = Vres;
    spikeCounts(spikeInd) = spikeCounts(spikeInd) + 1;
end

mesh(ex, in, spikeCounts/T);
set(gca, 'FontSize', 20)
xlabel('excitation', 'FontSize', 20)
ylabel('inhibition', 'FontSize', 20)
zlabel('rate', 'FontSize', 20)

% n = 200;
% sg = LIFSpikeGenerator(.001, .002, .02, -1+2*rand(1,n), 100+100*rand(1,n), 0);
% cp = CosinePopulation([1], sg, 'test');  
% do = addOrigin(cp, 'x', @(x) sin(1*pi*x));
% % do = addOrigin(cp, 'x', @(x) x);
% 
% x = -1:.1:1;
% rates = getRates(cp, x, 0, 0);
% sumEx = max(0, do.decoders) * rates;
% sumIn = min(0, do.decoders) * rates;
% 
% figure, plot(sumEx, -sumIn, 'k')
% axis equal
% set(gca, 'XLim', [0 max(sumEx)]), set(gca, 'YLim', [0 max(sumEx)])
% set(gca, 'FontSize', 20)
% xlabel('Excitation', 'FontSize', 20)
% ylabel('Inhibition', 'FontSize', 20)
% 
% figure, hold on
% plot(x, sumEx, 'r')
% plot(x, sumIn, 'b')
