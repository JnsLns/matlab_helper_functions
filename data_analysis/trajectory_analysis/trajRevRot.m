
% function transformReversed = trajRevRot(dataPoints,usedRotateToAxis,coordCols,reverseRefDirs,trajStart)
%
% Reverses the transformation performed in a call of trajRot, which may
% include a translation and a rotation. This reverse transformation can be
% applied to either the trajectories themselves or any other data points,
% enabling, e.g., to first derive new data points from trajectories recti-
% fied with trajRot, and then transform the new points back into the original
% coordinate space of the raw trajectories. The output of this function is
% a cell array with each element holding a matrix of reverse-transformed
% data points (either trajectories or other), each cell corresponding to
% one of the trajectories supplied to trajRot in the intial transformation.
% Note that the output includes only the transformed data from the columns 
% given in coordCols, not other columns of the argument dataPoints.
%
% dataPoints: n-by-m cell array holding in each element one set of data
% points (or a single point), e.g., trajectories (max. three dimensions),
% in the form of an n-by-m matrix: In this matrix, n = data points, and
% m = point coordinates. dataPoints may also be supplied as an n-by-m
% matrix of that form holding only one set of points. The order of columns
% (x-,y-,z-values) must conform to that specified in coordCols (and with
% that used in the original trajectories).
%
% usedRotateToAxis: scalar value denoting the axis to which trajectories
% were rotated (i.e., on which their first and last points ended up) in the
% call to trajRot that produced the data in reverseRefDirs and trajStart
% (see below). The value refers to the respective element of coordCols (i.e.
% to the axis along which the values live that are supplied in the column
% coordCols(usedRotateAxis) of the dataPoint matrices).
%
% coordCols: Two- or three-element row vector, giving the numbers of the
% columns of the matrices within dataPoints that hold the coordinates of the
% points that should be reverse-transformed. For valid reversal, this mapping
% must be consistent with that used in the initial transform using trajRot.
% (i.e., if the numbers of x,y,z-columns in the input matrices are the same
% as in the trajectories passed to trajRot, then the same coordCols can be
% used).
%
% reverseRefDirs: matrix supplied by trajRot as output of the same name (in
% each row holding one reference direction that when being "coupled" to the
% corresponding rotated trajectory and then rotated to the axis used in the
% original transformation would result in the original unrotated trajectory.
% Here it is used to rotate arbitratry data points to the coordinate frame
% of the original trajectory)
%
% trajStart (optional, by default no translation of coordinates is
% is done, which is appropriate in case shiftToOrigin was set to 0 in the
% call to trajRot that produced the data supplied in the other arguments of
% the current function call) is a cell array with each element holding one row
% vector with numel(coordCols) elements, which provides the original starting
% point of one trajectory before it was transformed using trajRot. Based on
% this value, the set of dataPoints in the corresponding element of
% dataPoints is translated equivalent to a reversal of the coordinate shift
% performed by trajRot (in case shiftToOrigin was 1). Usually, the
% corresponding output provided by trajRot can simply be passed on to
% trajRevRot.


function transformReversed = trajRevRot(dataPoints,usedRotateToAxis,coordCols,reverseRefDirs,trajStart)

if isa(dataPoints,'double')
    dataPoints = {dataPoints};
end

% ensure that dataPoints has only one column (avoids problems with cellfun
% outputs, permute, and in case of interaction with other functions)
dataPoints = reshape(dataPoints,numel(dataPoints),1);

% default for trajStart is list of points [0,0,...]
if nargin < 5
    trajStart = mat2cell(repmat(zeros(1,numel(coordCols)),numel(dataPoints),1),ones(1,numel(dataPoints)));
end

rotationReversed = trajRot(dataPoints,usedRotateToAxis,coordCols,reverseRefDirs,0);

transformReversed = cellfun(@(dat,start) dat(:,coordCols) + start(ones(size(dat,1),1),:), ...
    rotationReversed, trajStart, 'uniformoutput', 0);


end