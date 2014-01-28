%function [retinaOutput] = retina(movie, gain, threshold)
% linearInput, nonlinearOutput, spikes
% INPUT: number of ON parvocellular, number of OFF parvocellular, number ON
% magnocellular, number OFF magnocellular, whether to use grating or
% spatial white noise
% OUTPUT: cell array of spike times

numOnP = 32;
numOffP = 32;
numOnM = 32;
numOffM = 32;

% params
resolution = 32;
peakFiringRate = 10;

image_height = size(movie,1);
image_width  = size(movie,2);
max_rf_radius = 10;
center_rate_density = 0.5;
T = size(movie,3);
Ts = 2.5e-3;
stimulus = zeros(image_height,image_width,T/Ts);


for i = 1:T/Ts
    stimulus(:,:,i) = movie(:,:,ceil(i*Ts)); 
end

% Make Center-surround filters
filterHeightP = 7;
filterWidthP = 7;
centerRadiusP = 2;
surroundRadiusP = 4;
on_centerP = center_surround(filterHeightP,filterWidthP,centerRadiusP,surroundRadiusP,center_rate_density);
off_centerP = center_surround(filterHeightP,filterWidthP,centerRadiusP,surroundRadiusP,-center_rate_density);

filterHeightM = 15;
filterWidthM = 15;
centerRadiusM = 3;
surroundRadiusM = 7;
on_centerM = center_surround(filterHeightM,filterWidthM,centerRadiusM,surroundRadiusM,center_rate_density);
off_centerM = center_surround(filterHeightM,filterWidthM,centerRadiusM,surroundRadiusM,-center_rate_density);

% Make biphasic in-time filter
kernel = linearKernel(1,-0.5,1,resolution); % inputs: freq, phase, var, resolution

retinaOutput = struct;
retinaOutput.time = 0:Ts:T;
retinaOutput.stimulus = stimulus;

% Large ON Cells
outputMon = cell(numOnM,1);
locations = zeros(2,numOnM);
for i = 1:numOnM
    spatialFiltered = zeros(image_height,image_width,T/Ts);
    output = zeros(T/Ts,1);
    %  choose one particular segment of the image (x;y)
    location = floor(rand(2,1).*([image_height image_width]'-2*max_rf_radius)) + max_rf_radius;
    filterRadius = floor(filterHeightM/2);
    weights = zeros(image_height+filterHeightM,image_width+filterWidthM);
    weights(location(1):location(1)+filterHeightM-1,location(2):location(2)+filterWidthM-1) = on_centerM;
    weights = weights(filterRadius:image_height+filterRadius-1,filterRadius:image_width+filterRadius-1);
    
    
    for t = 1:T/Ts
        spatialFiltered(:,:,t) = stimulus(:,:,t).*weights;
        tmp = spatialFiltered(:,:,t);
        output(t) = sum(tmp(:));
    end
    
    
    % filter in time
    linearOutput = conv(output,kernel,'same'); % should this be 'full' to maintain causality, and then snip the end?
    nonlinearOutput = sigmoid(linearOutput,'gain',gain,'threshold',threshold,'maximum',peakFiringRate)
    nonlinearOutput = col(nonlinearOutput);
    
    
    % Poisson process
    spikeTrain = poissrnd(peakFiringRate*(Ts/T)*nonlinearOutput);
    moreSpikes = find(spikeTrain>1);
    spikeTrain(moreSpikes) = 1;
    
    outputMon{i} = spikeTrain;
    locations(:,i) = location;
    
    i
end

retinaOutput.Mon = outputMon;
retinaOutput.Mon_loc = locations;

% Large OFF Cells
outputMoff = cell(numOffM,1);
locations = zeros(2,numOffM);
for i = 1:numOnM
    spatialFiltered = zeros(image_height,image_width,T/Ts);
    output = zeros(T/Ts,1);
    %  choose one particular segment of the image (x;y)
    location = floor(rand(2,1).*([image_height image_width]'-2*max_rf_radius)) + max_rf_radius;
    filterRadius = floor(filterHeightM/2);
    weights = zeros(image_height+filterHeightM,image_width+filterWidthM);
    weights(location(1):location(1)+filterHeightM-1,location(2):location(2)+filterWidthM-1) = off_centerM;
    weights = weights(filterRadius:100+filterRadius-1,filterRadius:100+filterRadius-1);
    
    for t = 1:T/Ts
        spatialFiltered(:,:,t) = stimulus(:,:,t).*weights;
        tmp = spatialFiltered(:,:,t);
        output(t) = sum(tmp(:));
    end
    
    % filter in time
    linearOutput = conv(output,kernel,'same'); % should this be 'full' to maintain causality, and then snip the end?
    nonlinearOutput = threshold(linearOutput,point,slope);
    nonlinearOutput = col(nonlinearOutput);
    
    % Poisson process
    spikeTrain = poissrnd(meanFiringRate*(Ts/T)*nonlinearOutput);
    moreSpikes = find(spikeTrain>1);
    spikeTrain(moreSpikes) = 1;
    
    outputMoff{i} = spikeTrain;
    locations(:,i) = location;
    
    i
end

retinaOutput.Moff = outputMoff;
retinaOutput.Moff_loc = locations;

% Small ON Cells
outputPon = cell(numOnM,1);
locations = zeros(2,numOnP);
for i = 1:numOnM
    spatialFiltered = zeros(image_height,image_width,T/Ts);
    output = zeros(T/Ts,1);
    %  choose one particular segment of the image (x;y)
    location = floor(rand(2,1).*([image_height image_width]'-2*max_rf_radius)) + max_rf_radius;
    filterRadius = floor(filterHeightP/2);
    weights = zeros(image_height+filterHeightP,image_width+filterWidthP);
    weights(location(1):location(1)+filterHeightP-1,location(2):location(2)+filterWidthP-1) = on_centerP;
    weights = weights(filterRadius:100+filterRadius-1,filterRadius:100+filterRadius-1);
    
    for t = 1:T/Ts
        spatialFiltered(:,:,t) = stimulus(:,:,t).*weights;
        tmp = spatialFiltered(:,:,t);
        output(t) = sum(tmp(:));
    end
    
    % filter in time
    linearOutput = conv(output,kernel,'same'); % should this be 'full' to maintain causality, and then snip the end?
    nonlinearOutput = threshold(linearOutput,point,slope);
    nonlinearOutput = col(nonlinearOutput);
    
    % Poisson process
    spikeTrain = poissrnd(meanFiringRate*(Ts/T)*nonlinearOutput);
    moreSpikes = find(spikeTrain>1);
    spikeTrain(moreSpikes) = 1;
    
    outputPon{i} = spikeTrain;
    locations(:,i) = location;
    
    i
end

retinaOutput.Pon = outputPon;
retinaOutput.Pon_loc = locations;

% Small OFF Cells
outputPoff = cell(numOnM,1);
locations = zeros(2,numOffP);
for i = 1:numOnM
    spatialFiltered = zeros(image_height,image_width,T/Ts);
    output = zeros(T/Ts,1);
    %  choose one particular segment of the image (x;y)
    location = floor(rand(2,1).*([image_height image_width]'-2*max_rf_radius)) + max_rf_radius;
    filterRadius = floor(filterHeightP/2);
    weights = zeros(image_height+filterHeightP,image_width+filterWidthP);
    weights(location(1):location(1)+filterHeightP-1,location(2):location(2)+filterWidthP-1) = off_centerP;
    weights = weights(filterRadius:100+filterRadius-1,filterRadius:100+filterRadius-1);
    
    for t = 1:T/Ts
        spatialFiltered(:,:,t) = stimulus(:,:,t).*weights;
        tmp = spatialFiltered(:,:,t);
        output(t) = sum(tmp(:));
    end
    
    % filter in time
    linearOutput = conv(output,kernel,'same'); % should this be 'full' to maintain causality, and then snip the end?
    nonlinearOutput = threshold(linearOutput,point,slope);
    nonlinearOutput = col(nonlinearOutput);
    
    % Poisson process
    spikeTrain = poissrnd(meanFiringRate*(Ts/T)*nonlinearOutput);
    moreSpikes = find(spikeTrain>1);
    spikeTrain(moreSpikes) = 1;
    
    outputPoff{i} = spikeTrain;
    locations(:,i) = location;
    
    i
end

retinaOutput.Poff = outputPoff;
retinaOutput.Poff_loc = locations;

save retinaOutput.mat retinaOutput

%{
%% calculate sample ganglion cell response
%generate spike trains using Poisson process with rate which is a function
%of:
%   type - midget/parasol (or simplified model)
%   on/off center
%   location
%   size
% subplot(122)
r_max = 20;
g_1 = 0.5;
L_half = 0;

on_gang_rates = imfilter(stimulus,on_centerP);
off_gang_rates = imfilter(stimulus,off_centerP);
on_gang_rates_nl = non_linearity(r_max,g_1,L_half,on_gang_rates);
off_gang_rates_nl = non_linearity(r_max,g_1,L_half,off_gang_rates);
% imshow(on_gang_rates)

N = T/Ts;
rgc_response_ON = zeros(image_height, image_width, N);
rgc_response_OFF = zeros(image_height, image_width, N);

for i = 1:image_height
    for j = 1:image_width
        rgc_response_ON (i,j,:) = rate_encode(on_gang_rates_nl(i,j),T,Ts);
        rgc_response_OFF(i,j,:) = rate_encode(off_gang_rates_nl(i,j),T,Ts);
    end
end


% imshow((result-min(min(result)))/(max(max(result))-min(min(result))));


% connect ganglion cells together



%% calculate sample lgn response
%magnocellular/parvocellular (or simplified model)
%for now simply relay
lgn_response_ON = rgc_response_ON;
lgn_response_OFF = rgc_response_OFF;


%% calculate sample V1 response


% plot(v1)
%% 
%}
