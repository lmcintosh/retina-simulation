function [spikeTrain,nonlinearOutput,spikeTimestamps,time] = lnp(input, varargin)
% lnp is a linear-nonlinear-poisson cascade model used to model retinal ganglion neurons.
% REQUIRED INPUTS: stimulus (vector)
% OPTIONAL INPUTS: 'binLength', 'gain', 'threshold', 'peakFiringRate', 'meanFiringRate', and 'plots' flag
% OUTPUTS: spikeTrain, stimulus, nonlinearOutput

p = inputParser;
addRequired(p, 'input',@isnumeric);
addParamValue(p, 'binLength',1,@isnumeric);
addParamValue(p, 'gain',1,@isnumeric);
addParamValue(p, 'threshold',1,@isnumeric);
addParamValue(p, 'peakFiringRate',1,@isnumeric);
addParamValue(p, 'meanFiringRate',1.1,@isnumeric);
addParamValue(p, 'plots',1,@isnumeric);
parse(p, input, varargin{:});

% create the linear filter
% 1, -0.5, 1, 30 are parameters for a biphasic ON cell
%kernel = linearKernel(1,-0.5,1,100); % inputs: freq, phase, var, resolution
kernel = linearKernel(1,-0.5,1,30);

% pass the stimulus through the linear filter
linearOutput = conv(p.Results.input,kernel,'same'); % should this be 'full' to maintain causality, and then snip the end?

% pass the linear filter's output through the nonlinearity/threshold
nonlinearOutput = sigmoid(linearOutput,'gain',p.Results.gain,'threshold',p.Results.threshold,'maximum',p.Results.peakFiringRate);
nonlinearOutput = col(nonlinearOutput);
normConstant    = sum(nonlinearOutput);
nonlinearOutput = p.Results.meanFiringRate*nonlinearOutput/normConstant;

% use poisson probability of spiking to generate a spike train
spikeTrain = poissrnd(length(p.Results.input)*p.Results.binLength*nonlinearOutput);


time   = linspace(0,p.Results.binLength*length(input),length(nonlinearOutput));
ts     = 0:length(spikeTrain)-1;
spikes = ts(find(spikeTrain));

spikeTimestamps = time(spikes);


if p.Results.plots ~= 0
    figure; subplot(4,1,1), plot(time,p.Results.input,'k'), title('Stimulus'),
    subplot(4,1,2), plot(time,linearOutput,'k','LineWidth',1), title('Linear Output'),
    subplot(4,1,3), plot(time,nonlinearOutput,'k','LineWidth',1), title('Nonlinear Output'),
    subplot(4,1,4), ...
        for i = 1:length(spikes)
            line([spikes(i),spikes(i)],[0,1],'Color','k');
        end
end

%rasters(spikeTrain)
