% Sigmoid nonlinearity function
% Author: Niru Maheswaranathan
% Edited: Lane McIntosh
% Thu Jul 26 12:27:22 2012
% f = sigmoid(x, varargin)

function f = sigmoid(input, varargin)

p = inputParser;
addRequired(p,'input',@isnumeric);
addParamValue(p,'gain',1,@isnumeric);
addParamValue(p,'threshold',0,@isnumeric);
addParamValue(p,'maximum',1,@isnumeric);
parse(p,input,varargin{:});

f = p.Results.maximum./(1+exp(-p.Results.gain*(p.Results.input-p.Results.threshold)));
