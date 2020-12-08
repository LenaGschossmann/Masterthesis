function [image, imInfo] = load_bf_file(filepointer, readmetadata)
%% Load bioformat file and read metadata
% Sources:
% Melissa Linkert, Curtis T. Rueden, Chris Allan, Jean-Marie Burel, Will
% Moore, Andrew Patterson, Brian Loranger, Josh Moore, Carlos Neves,
% Donald MacDonald, Aleksandra Tarkowska, Caitlin Sticco, Emma Hill,
% Mike Rossner, Kevin W. Eliceiri, and Jason R. Swedlow (2010) Metadata
% matters: access to image data in the real world. The Journal of Cell
% Biology 189(5), 777-782. doi: 10.1083/jcb.201004104
%
% Swedlow JR, Goldberg I, Brauner E, Sorger PK (2003) Informatics and
% quantitative analysis in biological imaging. Science 300(5616), 100-2. Published 4 April 2003


% addpath('bfmatlab');

% Read image file
data = bfopen(filepointer);

imInfo = struct();
imInfo.width_px = size(data{1,1}{1,1},2);
imInfo.height_px = size(data{1,1}{1,1},1);
imInfo.planes = size(data{1,1},1);

% Create image stack
image = zeros(imInfo.height_px,imInfo.width_px,imInfo.planes,'uint16');
for iSeq = 1:imInfo.planes
    image(:,:,iSeq) = data{1,1}{iSeq,1};
end

% Read metadata
if readmetadata
    metadata = data{1, 2};
    metadataKeys = metadata.keySet().iterator();
    num_metadata = metadata.size();
%     information = {'',''};
%     num_metadata = size(information,1);
    metainfo = cell(num_metadata,1);
    for i=1:num_metadata
        key = metadataKeys.nextElement();
%         key = information{i,1};
        value = metadata.get(key);
        metainfo{i,1} = sprintf('%s = %s\n', key, value);
    end
    imInfo.metadata = metainfo;
    
    % OME metadata
    omeMeta = data{1, 4};
%     imInfo.voxelSizeXdefaultValue = omeMeta.getPixelsPhysicalSizeX(0).value();           % returns value in default unit
    imInfo.vx_scale_unit = omeMeta.getPixelsPhysicalSizeX(0).unit().getSymbol(); % returns the default unit type
    vx_x_scaled = omeMeta.getPixelsPhysicalSizeX(0).value(ome.units.UNITS.MICROMETER); % in µm
    imInfo.vx_x_scaled = vx_x_scaled.doubleValue();                                  % The numeric value represented by this object after conversion to type double
    vx_y_scaled = omeMeta.getPixelsPhysicalSizeY(0).value(ome.units.UNITS.MICROMETER); % in µm
    imInfo.vx_y_scaled = vx_y_scaled.doubleValue();                                  % The numeric value represented by this object after conversion to type double
    try
        voxelSizeZ = omeMeta.getPixelsPhysicalSizeZ(0).value(ome.units.UNITS.MICROMETER); % in µm
        voxelSizeZdouble = voxelSizeZ.doubleValue();                                  % The numeric value represented by this object after conversion to type double
    end
end

end