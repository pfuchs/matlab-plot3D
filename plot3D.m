function im = plot3D(Data, Mask, options)
% Plotting function inspired by "im" from Jeff Fesslers MIRT toolbox
% automatically crops the image domain to the mask and plots mosaic or
% midplane orientations depending on the inputs given.
%
% Also allows maximum intensity projections through the plot type "mip3"
%
% Example
%   [X,Y,Z] = ndgrid(-10:10);
%   R = (X.^2 + Y.^2 + (Z+4).^2);
%   sp = R < 100;
%   im = sp.*(rand(size(sp))/5+exp(-R/25));
%   figure(); subplot(211)
%   plot3D(im); % plot without masking
%   subplot(212);
%   plot3D(im, sp); % plot with masking
%
% see also IMAGE, IMAGESC, REGIONPROPS3, ALPHA

%% Input validation and defaults
arguments
    Data (:,:,:,:) double {mustBeNumeric}
    Mask (:,:,:,:) double {mustBeNumeric} = ones(size(Data))
    options.Type (1,1) string {mustBeMember(options.Type, ...
              ["mosaic","square","line","stack","mid3","mip3"])} = "mosaic"
    options.Crop (1,1) logical = true
    options.FigHandle (1,1) = 1
end

figure(FigHandle)
%% Preprocessing
if options.Crop
    bbox = floor(getfield(regionprops3(Mask ~=0 ,'BoundingBox'),'BoundingBox'));
    bboxIdx = substruct('()',{bbox(1,2)+(1:bbox(1,5)), ...
                              bbox(1,1)+(1:bbox(1,4)), ...
                              bbox(1,3)+(1:bbox(1,6))});
    
    Data = subsref(Data, bboxIdx);
    Mask = subsref(Mask, bboxIdx);
end

%% Extract image mosaic
switch options.Type
    case {"line", "stack", "square", "mid3"}
        imMat = mid3(Data, options.Type);
        mskMat = mid3(Mask, options.Type);
    case {"mosaic"}
        imMat = mosaic(Data);
        mskMat = mosaic(Mask);
    case {"mip3"}
        imMat = mip3(Data);
        mskMat = mip3(Mask);
    otherwise
        error('No valid image type given.')
end

im = image(imMat, 'CDataMapping', 'scaled'); 

im.AlphaData = mskMat;
im.Parent.Color = 'none';
im.Parent.Box = 'off';
 
axis image off;
colormap(gray);

colorbar;
% clim = min(abs(caxis));
% caxis([-clim, clim]);

title(inputname(1))

end

function mosaic = mid3(Data, varargin)
% Create mosaic of midplanes in xy/xz/zy axis.

type = 'square';
if nargin>1
    type = varargin{1};
end

sz = size(Data);

% Extract midplanes in the right orientation
xy = flipud(permute(Data(:,:,ceil(end/2)), [2 1 3]));
xz = flipud(permute(Data(:,ceil(end/2),:), [3 1 2]));
zy = flipud(permute(Data(ceil(end/2),:,:), [3 2 1]));

% Combine them into a tiled mosaic
switch lower(type)
    case {'line'}
        ndz = sz(2) - sz(3);
        nxy = sz(1) + sz(2);
        mosaic = [xy, [nan(floor(ndz/2),nxy); ...
                       [xz, zy]; ...
                       nan(ceil(ndz/2),nxy)]];
    case {'square', 'mid3'}
        mosaic = [[xy; xz], [nan(sz(2),sz(2)); zy]];
    case {'stack'}
        ndx = sz(2) - sz(1);
        nyz = sz(2) + sz(3);
        mosaic = [[nan(nyz,floor(ndx/2)),[xy;xz], nan(nyz,ceil(ndx/2))];...
                       zy];
    otherwise
        error('No valid mosaic type given, valid types are "line" and "square"');
end

end

function mosaic = mosaic(Data, varargin)
% Inspired by "montager" from Jeff Fesslers MIRT toolbox
% Creates a mosaic from 3 or 4 dimensional data, slicing axially along the
% outermost dimension optional input limits, gives column and row limits to
% reshape the mosaic

sz = size(Data);
nz = prod(sz(3:end));

aspectRatio = 1.4;
if sz(1) == sz(2) && nz == round(sqrt(nz))^2 % perfect square
			lims(1) = round(sqrt(nz));
else
    lims(1) = ceil(sqrt(nz * sz(1) / sz(2) * aspectRatio));
end
lims(2) = ceil(nz / lims(1));

if nargin>1
    lims = varargin{1};
end

if ndims(Data) > 3
	Data = reshape(Data, [sz(1) sz(2) nz]);
end

mosaic = zeros(sz(1) * lims(1), sz(2) * lims(2));
for iz=0:(nz-1)
	iy = floor(iz / lims(1));
	ix = iz - iy * lims(1);
	tmp = Data(:,:,iz+1);

	mosaic((1:sz(1))+ix*sz(1), (1:sz(2))+iy*sz(2)) = fliplr(tmp);
end

mosaic = mosaic';
end

function mosaic = mip3(Data)
% Create mosaic of the maximum intensity projections of the data in
% xy/xz/zy planes

sz = size(Data);

% Extract maximum intensity projections in the right orientation
xy = flipud(permute(max(Data,[],3), [2 1 3]));
xz = flipud(permute(max(Data,[],2), [3 1 2]));
zy = flipud(permute(max(Data,[],1), [3 2 1]));

% Combine them into a tiled mosaic
mosaic = [[xy; xz], [nan(sz(2),sz(2)); zy]];

end