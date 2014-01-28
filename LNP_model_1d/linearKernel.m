function [kernel] = linearKernel(frequency,phase,variance,resolution)
% linearKernel generates the linear filter of an LN model typical of the retina
% This filter is the dot product of a sinusoid and a gaussian pdf.
% INPUTS: frequency (of the sinusoid), phase (of the sinusoid), variance (of the gaussian), resolution (how many samples do you want?)
% OUTPUTS: row vector of the kernel

x = linspace(0,2*pi,resolution);

sinusoid = sin(frequency*x + phase);
y = gaussian(0,variance,10e3);
gauss = hist(y,resolution);


kernel = col(sinusoid).*col(gauss);

kernel = kernel/max(abs(kernel));
