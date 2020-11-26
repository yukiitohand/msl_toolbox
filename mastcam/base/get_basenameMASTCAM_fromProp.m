function [basename_mastcam] = get_basenameMASTCAM_fromProp(prop)
% [basename_mastcam] = get_basenameMASTCAM_fromProp(prop)
%  Input Parameters
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
%  Output Parameters
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

sol            = prop.sol;
cam_code       = prop.cam_code;
seq_id         = prop.seq_id;
command_num    = prop.command_num;
cdpid_counter  = prop.cdpid_counter;
unique_cdpid   = prop.unique_cdpid;
product_type   = prop.product_type;
gop            = prop.gop;
data_proc_code = prop.data_proc_code;
vr             = prop.version;

if isnumeric(sol)
    sol = sprintf('%04d',sol);
end

if isnumeric(seq_id)
    seq_id = sprintf('%06d',seq_id);
end
if isnumeric(command_num)
    command_num = sprintf('%03d',command_num);
end
if isnumeric(cdpid_counter)
    cdpid_counter = sprintf('%02d',cdpid_counter);
end
if isnumeric(unique_cdpid)
    unique_cdpid = sprintf('%05d',unique_cdpid);
end
if isnumeric(gop)
    gop = sprintf('%1d',gop);
end

if isnumeric(vr)
    vr = sprintf('%1d',vr);
end


basename_mastcam = sprintf('%s%s%s%s%s%s%s%s%s_%s',sol,cam_code,seq_id,...
    command_num,cdpid_counter,unique_cdpid,product_type,gop,vr,...
    data_proc_code);

end

