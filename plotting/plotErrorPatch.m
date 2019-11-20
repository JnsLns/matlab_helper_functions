function hPatch = plotErrorPatch(hLine,errorVals,errorDim)
%
% Plot polygon around 2d curve, representing extent of error (similar to
% errorbars). Patch properties such as transparency (alpha) can later be
% adjusted using the returned patch handle.
%
% hLine: Handle to a line object existing in the workspace around which the
% errorbar should be plotted.
%
% errorVals: vector holding an error value for each data point in the line
% (i.e., in line's 'xdata' and 'ydata' properties). The patch will extend
% from the corresponding data points along dimension errorDim by the value
% specified in errorVals *in both directions* (i.e., total extent will be
% two times the value in errorVals).
%
% errorDim (optional; default 1, i.e. y dimension): define within which
% dimension the errors live, that is, which dimension "error bars" will
% span (the different "error bars" will be spread across the other dim.).
% 1 for errors spanning y-dimension (i.e., error values apply to the data 
% in the line object's 'YData' property).
% 2 for errors spanning x-dimension (i.e., error values apply to the data 
% in the line object's 'XData' property).


if nargin < 3
    errorDim = 1;
end

% make sure errorVals is column vector
errorVals = errorVals(:);

% Line data
x_data = get(hLine,'XData')';
y_data = get(hLine,'YData')';

% Axes of target line
hAxes = get(hLine,'parent');

if errorDim == 1 % errors span y dimension (different error bars along x)
    errBarPositions = x_data;
    errBarOrigins = y_data;
elseif errorDim == 2 % errors span x dimension (different error bars along y)
    errBarPositions = y_data;
    errBarOrigins = x_data;
end

% Remove nans from error values
errorVals(isnan(errorVals)) = 0;

x_vertexCoords = [errBarPositions;flipud(errBarPositions)];
y_vertexCoords = [errBarOrigins+errorVals;flipud(errBarOrigins)-flipud(errorVals)];

axes(hAxes); 
wasHold = ishold;
if ~wasHold
    hold(hAxes,'on');
end

% Plot patch  
hPatch = fill(x_vertexCoords, y_vertexCoords, get(hLine,'color'),'facealpha',0.07,'linestyle','none','parent',hAxes);
% Disable legend entry for patch
set(get(get(hPatch(end),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

if ~wasHold
    hold(hAxes,'off');
end

