function [bandpass] = mastcamBandpassRead(filter_name,varargin)
% [bandpass] = mastcamBandpassRead(filter_name)
%   Read bandpass function. Bandpass could be either lens or filter
%  In case of "Lens" and "L#", the bandpass has two fields: "wavelength" 
%  and "transmission". In case of "Red", "Green", or "Blue", the bandpass 
%  has two fileds: "wavelength" and "qe".
%
% Usage:
%  > [bandpass] = mastcamBandpassRead('Lens','R')
%  > [bandpass] = mastcamBandpassRead('R0')
%  > [bandpass] = mastcamBandpassRead('L1')
%  > [bandpass] = mastcamBandpassRead('Red')
%  > [bandpass] = mastcamBandpassRead('Green')
%  > [bandpass] = mastcamBandpassRead('Blue')

switch upper(filter_name)
    case 'LENS'
        cam_eye = varargin{1};
        fname = sprintf('M%1s_LENS_BANDPASS.txt',cam_eye);
    case {'L0','L1','L2','L3','L4','L5','L6','L7','R0','R1','R2','R3','R4','R5','R6','R7'}
        cam_eye = filter_name(1);
        filter_id = filter_name(2);
        fname = sprintf('M%1s_FILTER_%1s_BANDPASS.txt',cam_eye,filter_id);
    case {'RED','GREEN','BLUE'}
        fname = sprintf('M_FILTER_%s_BANDPASS.txt',filter_name);
    otherwise
        error('The input is invalid. Refer USAGE in the help');
end

fid = fopen(fname,'r');
column_header = fgetl(fid);
column_headers = regexp(column_header,'\s+','split');
data = textscan(fid,'%f %f');
fclose(fid);

bandpass.(column_headers{1}) = data{1};
bandpass.(column_headers{2}) = data{2};

switch upper(filter_name)
    case 'LENS'
        if ~strcmpi(column_headers{2},'transmission')
            error('Check the file, the column is not labeled transmission.');
        end
    case {'L0','L1','L2','L3','L4','L5','L6','L7','R0','R1','R2','R3','R4','R5','R6','R7'}
        bandpass.(column_headers{2}) = bandpass.(column_headers{2})/100;
        if ~strcmpi(column_headers{2},'transmission')
            error('Check the file, the column is not labeled transmission.');
        end
    case {'RED','GREEN','BLUE'}
        if ~strcmpi(column_headers{2},'qe')
            error('Check the file, the column is not labeled qe.');
        end
    otherwise
        
end



end

