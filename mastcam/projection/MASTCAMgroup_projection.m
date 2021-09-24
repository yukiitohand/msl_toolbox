classdef MASTCAMgroup_projection < handle
    % MASTCAMgroup_projection
    %  Storing the information of the projection results
    %  Properties
    %    XYZ : ENVIRasterMultBand class obj
    %      double image [L_mastcam x S_mastcam x 3]
    %      XYZ-coordinate of the center of pixel of the MASTCAM image in
    %      the defined geographical coordinate system.
    %      1st layer: X coordinate, 2nd layer: Y coordinate, 3rd layer: Z
    %      coordinate
    %    NE  : ENVIRasterMultBand class obj
    %      double image [L_mastcam x S_mastcam x 2]
    %      northing and easting of the center of pixel of the MASTCAM image
    %      in the defined geographical coordinate system.
    %      1st layer: northing, 2nd layer: easting.
    %    latlon : ENVIRasterMultBand class obj
    %      double image [L_mastcam x S_mastcam x 2]
    %      latitude and longitude of the center of pixel of the MASTCAM image
    %      in the defined geographical coordinate system.
    %      1st layer: latitude, 2nd layer: longitude in degree
    %    zenith : ENVIRasterSingleLayer class obj
    %      double image [L_mastcam x S_mastcam]
    %      The coordinate value in the zenith direction. Elevation or
    %      radius.
    %    range : ENVIRasterSingleLayer class obj
    %      double image [L_mastcam x S_mastcam]
    %      The range of the center of the pixels. in meter.
    %    surfplnc: ENVIRasterMultBand class obj
    %      double image [L_mastcam x S_mastcam x 4]
    %      parameters for the surface plane equation. Surface plane is
    %      defined as the triangle that intersects with central pixel pointing
    %      vectors. 1st, 2nd, 3rd layers: x,y,z component of the normal
    %      vector. 4th layer: plane constant
    %    emiang : ENVIRasterSingleLayer class obj
    %      double image [L_mastcam x S_mastcam]
    %      The emission angle of the center of the pixels. in degree.
    %    nn_msldem  : ENVIRasterMultBand class obj
    %      int16 image [L_mastcam x S_mastcam x 2]
    %      nearest neighbor pixels in the base MSLDEM image. The first page
    %      is the sample indices and the second is the line indexes.
    %    ref_msldem : ENVIRasterMultBand class obj
    %      indicate which triangle in the MSLDEMdata, the pixel is located.
    %      pages 1,2,3 are sample, line, and j. j is either 1 or 2, indicating
    %      which triangle at the (sample, line).
    properties
        MASTCAMdataObj
        MSLDEMCdata
        XYZ
        NE %
        latlon
        zenith
        range
        surfplnc
        emiang
        nn_msldem  % nearest neighbor pixel in MSLDEM
        ref_msldem % ref triangles in the 
    end
    
    methods
        function obj = MASTCAMgroup_projection()
        end
        
    end
end

