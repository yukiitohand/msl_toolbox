function [prop] = getProp_basenameMASTCAM(basename_mastcam,varargin)
% [prop] = getProp_basenameMASTCAM(basename_mastcam,varargin)
%   Get properties from the basename of MASTCAM data
%  Input Parameters
%   basename_mastcam: string, like
%    SSSSIIFFFFFFLLLXXCCCCCPGV_DXXX
%     SSSS: Four-digit sol number after landing day (which was defined as sol 0).
%     II: Two-digit camera code: ?ML? = Mastcam Left (M-34) and ?MR? = Mastcam Right (M-100).
%     FFFFFF: Six-digit sequence number identifier.
%     LLL: Three-digit command number within the sequence that corresponds to this image.
%     XX: Two-digit Camera Data Product Identifier (CDPID) counter that records the number of times this CDPID has been used over the lifetime of the mission.
%     CCCCC: Five-digit CDPID value, uniquely assigned by the camera to an image product. 
%     P: One-letter product type (see Table 3).
%     G: One-letter Group of Pictures (GOP) hexadecimal counter, for video sequences.
%     V: One-digit version number.
%     DXXX: Four-letter data processing code (see Table 13).
%     ZZZ: Three-letter file extension (typically, ?DAT? or ?IMG?).
%  Output Parameters
%   Output
%   prop: struct storing properties
%    'sol'               
%    'cam_code'         
%    'seq_id'           
%    'command_num'      
%    'cdpid_counter'
%    'unique_cdpid'
%    'product_type'     
%    'gop'              
%    'data_proc_code'   
%    'version'

[ prop_ori ] = create_propMASTCAMbasename();
[basenameptrn] = get_basenameMASTCAM_fromProp(prop_ori);

prop = regexpi(basename_mastcam,basenameptrn,'names');

if ~isempty(prop)
    prop.sol = str2num(prop.sol);
    prop.seq_id = str2num(prop.seq_id);
    prop.command_num = str2num(prop.command_num);
    prop.cdpid_counter = str2num(prop.cdpid_counter);
    prop.unique_cdpid = str2num(prop.unique_cdpid);
    prop.version = str2num(prop.version);
end


end