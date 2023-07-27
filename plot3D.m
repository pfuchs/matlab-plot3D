function [fig, im] = plot3D(Data, Mask, options)
% Plotting function inspired by "im" from Jeff Fesslers MIRT toolbox
% automatically crops the image domain to the mask and plots mosaic or
% midplane orientations depending on the inputs given.
%
% PLOT3D(Volume) plots the 3D data in volume in a mosaic with a colorbar.
% The plot will be titled by the variable name of "Volume"
%
% PLOT3D(Volume, mask) trims the content of volume to the mask boundaries, 
% and plots the data in a mosaic, with transparent background.
%
% Options can be supplied in parameter/value pairs as
%
% PLOT3D(...,'PropertyName',PropertyValue,...) where valid property names
% are type, crop, OrientationLabels, fighandle, caxis, location, label
%
% Also allows maximum intensity projections through the plot type "mip3"
%
% To access the plotted image data matrix and alpha map use the returned
% Image handle's "CData" and "AlphaMap" fields. This will be an array of
% Images for 4D datasets.
%
% Code has been verified to work with MATLAB versions
%       '9.14.0.2206163 (R2023a)' & '9.10.0.2015706 (R2021a) Update 7'
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
arguments %(Input)
    Data (:,:,:,:) double {mustBeNumeric}
    Mask (:,:,:,:) double {mustBeNumeric} = ones(size(Data))
    options.Type (1,1) string {mustBeMember(options.Type, ...
              ["mosaic","square","line","stack","mid3","mip3"])} = "mosaic"
    options.Crop (1,1) logical = true
    options.OrientationLabels (1,1) logical = false
    options.FigHandle (1,1) = gcf
    options.caxis = 'auto'
    options.Location (1,:) char {mustBeMember(options.Location, ...
                                 ['north','south','east','west'])} = 'east'
    options.Label (1,1) string = ""
    options.Slices (1,3) int8 = ceil(size(Data)/2)
end

% arguments (Output) % Doesn't work with older versiions of MATLAB
%     fig matlab.ui.Figure
%     im (1,:) matlab.graphics.primitive.Image
% end

im = repmat(image(0),1,size(Data,4));

if isa(options.FigHandle,'matlab.graphics.axis.Axes')
    fig = options.FigHandle; 
    tcl = fig;
else
    fig = figure(options.FigHandle);
    tcl = tiledlayout(fig,'flow','TileSpacing','compact', 'Padding', 'loose');
end

for tt = 1:size(Data,4)
%% Preprocessing
if options.Crop
    bbox = floor(getfield(regionprops3(Mask(:,:,:,min(tt,size(Mask,4))) ~=0 ,'BoundingBox'),'BoundingBox'));
    bboxIdx = substruct('()',{bbox(1,2)+(1:bbox(1,5)), ...
                              bbox(1,1)+(1:bbox(1,4)), ...
                              bbox(1,3)+(1:bbox(1,6)), ...
                              tt});
    
    PlotData = subsref(Data, bboxIdx);
    
    bboxIdx = substruct('()',{bbox(1,2)+(1:bbox(1,5)), ...
                              bbox(1,1)+(1:bbox(1,4)), ...
                              bbox(1,3)+(1:bbox(1,6)), ...
                              min(tt,size(Mask,4))});
    PlotMask = subsref(Mask, bboxIdx);
else
    PlotData = Data(:,:,:,tt);
    PlotMask = Mask(:,:,:,min(tt,size(Mask,4)));
end

%% Extract image mosaic
switch options.Type
    case {"line", "stack", "square", "mid3"}
        imMat = mid3(PlotData, options.Type, options.Slices);
        mskMat = mid3(PlotMask, options.Type, options.Slices);
    case {"mosaic"}
        imMat = mosaic(PlotData);
        mskMat = mosaic(PlotMask);
    case {"mip3"}
        imMat = mip3(PlotData);
        mskMat = mip3(PlotMask);
    otherwise
        error('No valid image type given.')
end

%% Plot
if ~isa(fig,'matlab.graphics.axis.Axes'); nexttile; end
im(tt) = image(imMat, 'CDataMapping', 'scaled',...
                      'AlphaData', mskMat); 

im(tt).Parent.Color = 'none';
im(tt).Parent.Box = 'off';


%% Add labels
if options.OrientationLabels
    sz = size(PlotData);
    pltsz = size(imMat);
    switch lower(options.Type)
        case {"square", "mid3"}
            xshift = ceil(pltsz(1)/40);
            yshift = ceil(pltsz(2)/40);
            % axial
            text(xshift,         sz(2)/2,            'R','Color','r','FontWeight','bold')
            text(sz(1) - xshift, sz(2)/2,            'L','Color','r','FontWeight','bold')
            text(sz(1)/2,        yshift,             'A','Color','r','FontWeight','bold')
            text(sz(1)/2,        sz(2)-yshift,       'P','Color','r','FontWeight','bold')
            % coronal
            text(xshift,         sz(3)/2+sz(2),      'R','Color','r','FontWeight','bold')
            text(sz(1) - xshift, sz(3)/2+sz(2),      'L','Color','r','FontWeight','bold')
            text(sz(1)/2,        sz(2)+yshift,       'S','Color','r','FontWeight','bold')
            text(sz(1)/2,        sz(2)+sz(3)-yshift, 'I','Color','r','FontWeight','bold')
            % sagittal
            text(sz(1)+sz(2)-xshift, sz(3)/2+sz(2), 'A','Color','r','FontWeight','bold')
            text(sz(1)+xshift,       sz(3)/2+sz(2), 'P','Color','r','FontWeight','bold')
            text(sz(2)/2+sz(1),      sz(2)+yshift,       'S','Color','r','FontWeight','bold')
            text(sz(2)/2+sz(1),      sz(2)+sz(3)-yshift, 'I','Color','r','FontWeight','bold')
        case {'line'}
            xshift = ceil(pltsz(1)/40);
            yshift = ceil(pltsz(2)/60);
            % axial
            text(xshift,         sz(2)/2,            'R','Color','r','FontWeight','bold')
            text(sz(1) - xshift, sz(2)/2,            'L','Color','r','FontWeight','bold')
            text(sz(1)/2,        yshift,             'A','Color','r','FontWeight','bold')
            text(sz(1)/2,        sz(2)-yshift,       'P','Color','r','FontWeight','bold')
            % coronal
            text(sz(1) + xshift,  sz(2)/2,            'R','Color','r','FontWeight','bold')
            text(2*sz(1) - xshift,sz(2)/2,            'L','Color','r','FontWeight','bold')
            text(sz(1)/2 + sz(1), yshift,             'S','Color','r','FontWeight','bold')
            text(sz(1)/2 + sz(1), sz(2)-yshift,       'I','Color','r','FontWeight','bold')
            % sagittal
            text(2*sz(1) + xshift,  sz(2)/2,            'P','Color','r','FontWeight','bold')
            text(2*sz(1) + sz(2) - xshift,sz(2)/2,      'A','Color','r','FontWeight','bold')
            text(sz(1)/2 + 2*sz(1), yshift,             'S','Color','r','FontWeight','bold')
            text(sz(1)/2 + 2*sz(1), sz(2)-yshift,       'I','Color','r','FontWeight','bold')      
        otherwise
            warning("Labels have not been defined for this projection (yet).")
    end
end

axis image off;
colormap(gray);
if ndims(Data)>3; title(sprintf('Echo %i',tt)); end

c = colorbar(im(tt).Parent, 'Location', [options.Location 'outside']);
c.Label.String = options.Label;

caxis(im(tt).Parent, options.caxis); %#ok<CAXIS> % Kept as caxis for backwards compatibility
end

title(tcl, inputname(1))

end

function mosaic = mid3(Data, varargin)
% Create mosaic of midplanes in xy/xz/zy axis.

type = 'square';
slices = ceil(size(Data)/2);
if nargin>1
    type = varargin{1};
end
if nargin>2
    slices = varargin{2};
end

sz = size(Data);

% Extract midplanes in the right orientation
xy = flipud(permute(Data(:,:,slices(3)), [2 1 3]));
xz = flipud(permute(Data(:,slices(2),:), [3 1 2]));
zy = flipud(permute(Data(slices(1),:,:), [3 2 1]));

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