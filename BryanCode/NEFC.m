% This code demonstrates tweaking of NEF decoding so that the spike rate of 
% a post-synaptic conductance neuron is closer to being a function of the 
% ideal decoded function. 
% 
% To clean up: 
% - something seems wrong in gradient calculation (see comments in code)
% - the conductance isn't in sensible units
classdef NEFC < handle 
    
    methods (Static)
        
        % Run this method for a demo with plots. Set fun in the code to
        % change the decoded function. 
        function runMe() 
            fun = @(x) sin(pi*x);
%             fun = @(x) sin(2*pi*x);
%             fun = @(x) x.^2;
            [x, preRates] = NEFC.setup(300);
            NEFC.decode(x, preRates, fun);
        end
        
        % Discretizes presynaptic variable x and samples spike rates of a 
        % presynaptic LIF ensemble. 
        % 
        % n: # of neurons in presynaptic population 
        % 
        % x: discrete values of represented variable
        % preRates: spikes rates of neurons in an example presynaptic
        %   population at above x values
        function [x, preRates] = setup(n)
            x = (-1:.025:1);
                        
            intercepts = -1 + 2*rand(1, n);
            maxRates = 80 + 40*rand(1, n);
            tauRef = .005;
            tauRC = .01;
            encoders = -1 + 2*(rand(1, n)>0.5);
            preRates = NEFC.getLIFRatesX(intercepts, maxRates, tauRef, tauRC, encoders, x);
        end
        
        % Find excitatory and inhibitory decoding weights for a single
        % post-synaptic neuron that represents the given function with
        % encoder 1. 
        % 
        % Implementation: Initial decoders are found using using NEF
        % methods. If fun is non-monotonic, rates may not be consistent for 
        % different occurances of the same fun(x), so gradient descent is 
        % used to make rates more consistent as a function of fun. The
        % signs are not changed during gradient descent (e.g. inhibitory
        % synapses stay inhibitory). 
        % 
        % x: discrete values of variable represented by presynaptic neurons
        % preRates: spike rates of presynaptic neurons at x
        % fun: function of x to decode
        % 
        % decoders: -ve values are inhibitory; +ve excitatory 
        function decoders = decode(x, preRates, fun)
            fsn = 'FontSize'; 
            fs = 18;
            
            decoders = NEFC.getNEFDecoders(x, preRates, fun);
            posDec = decoders >= 0;
            negDec = decoders < 0;
            decoders = NEFC.scaleDecodersForConductance(preRates, decoders, 80);
                                    
            [exTable, inTable, rateTable, dRdExTable, dRdInTable] = NEFC.getRateTable();
            figure(1), hold on
            mesh(exTable, inTable, rateTable), set(gca, fsn, fs)
            xlabel('excitatory conductance', fsn, fs), ylabel('inhibitory conductance', fsn, fs), zlabel('spike rate', fsn, fs)            
            
            f = fun(x);
            [exCond, inCond] = NEFC.getConductance(preRates, decoders);

            initialPostRate = NEFC.getLIFRateCond(exCond, inCond); 
            err = NEFC.getSimilarityError(initialPostRate, f); 
            refRate = initialPostRate - err;
            meanRefRate = mean(refRate);
            plot3(exCond, inCond, initialPostRate, 'k')

            kappa = .0005;
            nIter = 400; 
            for i = 1:nIter
                pr = max(0, preRates + (i<nIter)*0*randn(size(preRates))); %we can add a little noise here but we drop it in last iteration for clean plot
                [exCond, inCond] = NEFC.getConductance(pr, decoders);
                postRate = interp2(exTable, inTable, rateTable, exCond, inCond);
                
                % recalculate reference but keep mean rate constant ...
                err = NEFC.getSimilarityError(postRate, f); 
                refRate = postRate - err;
                refRate = refRate * meanRefRate / mean(refRate);
                % ... end recalculate reference
                
                err = postRate - refRate;
                
                figure(2)
                subplot(1,2,1), plot(f, postRate, 'r', f, initialPostRate, 'k', f, refRate, 'gx-')
                set(gca, fsn, fs)
                xlabel('ideal decoded function', fsn, fs), ylabel('spike rate (red) & ideal rate (green)', fsn, fs)
                subplot(1,2,2), plot(x, err, 'r')
                set(gca, fsn, fs)
                xlabel('presynaptic variable x', fsn, fs), ylabel('postsynaptic spike rate error', fsn, fs)
                
                dRdEx = interp2(exTable, inTable, dRdExTable, exCond, inCond);
                dRdIn = interp2(exTable, inTable, dRdInTable, exCond, inCond);
                dExdD = pr .* repmat(posDec, 1, length(x)); % derivative of excitatory conductance wrt decoders (#neurons x #samples)
                dIcdD = pr.* repmat(negDec, 1, length(x)); % derivative of inhibitory conductance wrt decoders (#neurons x #samples)

                dRdD = dExdD .* repmat(dRdEx, size(preRates,1), 1) + dIcdD .* repmat(dRdIn, size(preRates,1), 1);

                % this was a mistake but works pretty well ... 
                dRdDScaled = dRdD ./ repmat(sum(dRdD.^2,1).^.5, size(preRates,1), 1) / size(preRates,1);
                dD = -kappa * (dRdDScaled*err'); 

                % this seems right to me but doesn't work as well as the above mistake ... 
%                 dD = - maxRateChange / size(preRates,1) * (dRdD ./ repmat(sum(dRdD.^2,2).^.5, 1, size(preRates,2))) * (err' / norm(err)) ./ mean(dRdD,2);

                decoders = decoders + dD;
                decoders(decoders < 0 & posDec) = 0; 
                decoders(decoders > 0 & negDec) = 0; 
            end
            figure(1), plot3(exCond, inCond, postRate, 'r')
        end

        % Scales NEF decoders so that they make sense in a conductance context. 
        % 
        % preRates: rates of presynaptic neurons at various x
        % decoders: weights of synapses onto a postsynaptic neuron 
        %   (+ve: excitatory; -ve: inhibitory)
        % maxExCond: maximum excitatory conductance
        % 
        % decoders: decoders scaled so that e.g. excitatory conductance = 
        %   preRates' * (decoders .* decoders>0)
        function decoders = scaleDecodersForConductance(preRates, decoders, maxExCond)
            [exCond, inCond] = NEFC.getConductance(preRates, decoders);

            exScale = maxExCond / max(exCond);
            inScale = 50/15 * exScale; % hard-coded from conductance model: (Eex-Eth) / (Ein-Eth) 
            
            posDec = decoders >= 0;
            negDec = decoders < 0;
            
            decoders(posDec) = decoders(posDec) * exScale;
            decoders(negDec) = decoders(negDec) * inScale;
        end
        
        % exCond: Excitatory conductance (can be a scalar, vector, or
        %   matrix)
        % inCond: Inhibitory conductance
        % 
        % dRdEx: Derivative of spike rate wrt excitatory conductance
        % dRdIn: Derivative of spike rate wrt inhibitory conductance
        function [dRdEx, dRdIn] = getLIFRateDerivative(exCond, inCond)
            dCond = 1; % conductance step with which to estimate derivative
            
            refRate = NEFC.getLIFRateCond(exCond, inCond);
            exRate = NEFC.getLIFRateCond(exCond+dCond, inCond);
            inRate = NEFC.getLIFRateCond(exCond, inCond+dCond);
            
            dRdEx = exRate - refRate;
            dRdEx(dRdEx == 0) = 1e-3;
            
            dRdIn = inRate - refRate;
            dRdIn(dRdIn== 0) = -1e-3;
        end 
        
        % Calculates a relative decoding error in terms of the differences
        % between firing rates of a post-syaptic neuron at different x for
        % which the decoded function is supposed to be the same (so the
        % firing rate should be the same). This doesn't impose a particular
        % decoded function -> firing rate relationship, so the tuning curve
        % could potentially look strange with zero error, as long as it is
        % internally consistent. 
        % 
        % r: postsynaptic spike rates (#samples)
        % f: ideal values of function to be decoded (#samples)
        % 
        % error: difference between r and some compromize that is
        %   a consistent function of f (e.g. mean of r over instances of
        %   same f)
        function error = getSimilarityError(r, f)
            error = zeros(size(f));
            for i = 1:length(f)
                lower = f < f(i);
                crossingInd = find(lower(1:end-1) ~= lower(2:end));
                if abs(f(1)-f(i))<1e-10, crossingInd = union(crossingInd, 1); end % whether it is already there depends on slope
                if abs(f(end)-f(i))<1e-10, crossingInd = union(crossingInd, length(f)-1); end
                if isempty(crossingInd), crossingInd = i; end

                n = length(crossingInd);
                rCrossing = zeros(1, n);
                for j = 1:n  % interpolate rate at each point f crossing f(i)
                    cj = crossingInd(j);
                    rCrossing(j) = r(cj) + (r(cj+1)-r(cj)) * (f(i)-f(cj)) / (f(cj+1)-f(cj));
                end
                meanr = mean(rCrossing);
                error(i) = r(i) - meanr;
            end
        end
        
        % x: discrete values of presynaptic represented variable
        % rates: spike rates of presynaptic neurons at x
        % fun: function of x to decode
        % 
        % decoders: NEF decoders
        function decoders = getNEFDecoders(x, rates, fun)
            f = fun(x)';
            gamma = rates * rates'; 
            decoders = pinv(gamma + .05*mean(diag(gamma))*eye(size(gamma))) * (rates*f);
        end

        % Returns spike rates of LIF neurons given the usual parameters, as
        % a function of represented variable x. 
        function rates = getLIFRatesX(intercepts, maxRates, tauRef, tauRC, encoders, x)
            scales = (1 ./ (1 - exp( (tauRef - (1 ./ maxRates)) / tauRC)) - 1) ./ (1 - intercepts);
            biases = 1 - scales .* intercepts;

            drive = encoders' * x;
            current = repmat(scales', 1, length(x)) .* drive + repmat(biases', 1, length(x)); 
            rates = zeros(size(current));
            activeIndices = find(current > 1);
            rates(activeIndices) = 1 ./ (tauRef - tauRC * log(1 - 1./current(activeIndices)));  
        end
        
        % Returns lookup tables for conductance-neuron spike rates and
        % derivatives. 
        % 
        % Use like this: rate = interp2(exTable, inTable, rateTable, ex, in)
        function [exTable, inTable, rateTable, dRdExTable, dRdInTable] = getRateTable()
            dEx = 2;
            dIn = 8;
            exLevels = 0:dEx:100;
            inLevels = 0:dIn:400;
            exTable = repmat(exLevels, length(inLevels), 1);
            inTable = repmat(inLevels', 1, length(exLevels));
            rateTable = NEFC.getLIFRateCond(exTable, inTable); 
            dRdExTable = [rateTable(:,2)-rateTable(:,1) rateTable(:,3:end)-rateTable(:,1:end-2) rateTable(:,end)-rateTable(:,end-1)] / dEx; 
            dRdExTable(dRdExTable <= 0) = 0.1;
            dRdInTable = [rateTable(2,:)-rateTable(1,:); rateTable(3:end,:)-rateTable(1:end-2,:); rateTable(end,:)-rateTable(end-1,:)] / dIn; 
            dRdInTable(dRdInTable >= 0) = -0.1;
        end
        
        % Returns spike rate of single LIF neuron as function of excitatory
        % and inhibitory conductance. 
        % 
        % exCond: excitatory conductance (can be vector or matrix) 
        % inCond: ihibitory conductance (can be vector or matrix) 
        % 
        % rate: spike rate
        function rate = getLIFRateCond(exCond, inCond)
            % we find this by simulation, because the current for given
            % conductance depends on the membrane potential, which varies
            % over time (we have to integrate, and simulation is a good
            % way to do so)
            
            dt = .0001; %small step needed for fine rate variations
            T = 3;
            time = dt:dt:T;

            Vres = -65; %potential to which cell is reset after spike
            Vth = -50; %spike threhold potential 
            EEx = 0; %excitatory reversal potential
            EIn = -65; %inhibitory reversal potential
            Eleak = -65; %leak current reversal potential
            gLeak = 0; %constant leak conductance
            Tref = .002; %spike refractory time
%             C = 1; %we set membrane capacitance=1 so it doesn't show up in the equations

            V = zeros(size(exCond));
            spikeCounts = zeros(size(exCond));
            lastSpikeTime = -ones(size(exCond));

            for i = 1:length(time)
                dVdt = -exCond.*(V-EEx) - inCond.*(V-EIn) - gLeak*(V-Eleak);
                refractory = find(time(i)-lastSpikeTime < Tref); 
                dVdt(refractory) = 0;
                V = V + dt*dVdt;
                spikeInd = find(V>=Vth);
                lastSpikeTime(spikeInd) = time(i);
                V(spikeInd) = Vres;
                spikeCounts(spikeInd) = spikeCounts(spikeInd) + 1;
            end
            
            rate = spikeCounts / T;
        end
        
        % rates: spike rates of presynaptic neurons
        % decoders: conductance decoders
        % 
        % exCond: excitatory conductance
        % inCond: inhibitory conductance
        function [exCond, inCond] = getConductance(rates, decoders)
            exCond = (rates' * (decoders .* (decoders>0)))';
            inCond = -(rates' * (decoders .* (decoders<0)))';
        end
        
    end
end