function [result] = threshold(x,point,slope)
% threshold.m thresholds your vector x with a piecewise function that is nonlinear at the point "point".
% INPUTS: x, point (the x value below which the result is zero), slope (the slope of the line after point)
% OUTPUTS: result (col vector)

result = zeros(length(x),1);

for i = 1:length(x)
    if x(i) >= point
        result(i) = slope*(x(i) - point);
    end
end
