function [ prop ] = create_propMASTCAMbasename( varargin )
% [ prop ] = create_propMASTCAMbasename( varargin )
%   return struct of MASTCAM product property
% 
%   Output
%   prop: struct storing properties
%    'sol'               : (default) '(?<sol>([\d]{4}|DEV_|TVC_|CAL_|DEL_|ATL_CRU_))'
%    'cam_code'          : (default) '(?<cam_code>M[RLHD]{1}'
%    'seq_id'            : (default) '(?<seq_id>[\d]{6})'
%    'command_num'       : (default) '(?<command_num>[\d]{3})'
%    'cdpid_counter'     : (default) '(?<cdpid_counter>[\d]{2})'
%    'unique_cdpid'      : (default) '(?<unique_cdpid>[\d]{5})'
%    'product_type'      : (default) '(?<product_type>[a-uA-U]{1})'
%    'gop'               : (default) '(?<gop>[\da-fA-F]{1})'
%    'data_proc_code'    : (default) '(?<data_proc_code>[DRCLX]{4})'
%    'version'           : (default) '(?<version>[\d]{1})'
%
%   Optional Parameters
%    'SOL', 'CAM_CODE', 'SEQ_ID', 'COMMAND_NUM', 'CDPID_COUNTER',
%    'UNIQUE_CDPID', 'PRODUCT_TYPE','GOP','DATA_PROC_CODE','VERSION'
% 
%  REFERENCE
%   SSSSIIFFFFFFLLLXXCCCCCPGV_DXXX.ZZZ defined as follows:
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

sol            = '(?<sol>([\d]{4}|DEV_|TVC_|CAL_|DEL_|ATL_CRU_))';
cam_code       = '(?<cam_code>M[RLHD]{1})';
seq_id         = '(?<seq_id>[\d]{6})';
command_num    = '(?<command_num>[\d]{3})';
cdpid_counter  = '(?<cdpid_counter>[\d]{2})';
unique_cdpid   = '(?<unique_cdpid>[\d]{5})';
product_type   = '(?<product_type>[a-uA-U]{1})';
gop            = '(?<gop>[\da-fA-F]{1})';
data_proc_code = '(?<data_proc_code>[DRCLX]{4})';
vr             = '(?<version>[\d]{1})';

if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        if ~isempty(varargin{i+1})
            switch upper(varargin{i})
                case 'SOL'
                    sol = varargin{i+1};
                case 'CAM_CODE'
                    cam_code = varargin{i+1};
                case 'SEQ_ID'
                    seq_id = varargin{i+1};
                case 'COMMAND_NUM'
                    command_num = varargin{i+1};
                case 'CDPID_COUNTER'
                    cdpid_counter = varargin{i+1};
                case 'UNIQUE_CDPID'
                    unique_cdpid = varargin{i+1};
                case 'PRODUCT_TYPE'
                    product_type = varargin{i+1};
                case 'GOP'
                    gop = varargin{i+1};
                case 'DATA_PROC_CODE'
                    data_proc_code = varargin{i+1};
                case 'VERSION'
                    vr = varargin{i+1};
                otherwise
                    error('Unrecognized option: %s', varargin{i});   
            end
        end
    end
end

prop = [];
prop.sol = sol;
prop.cam_code = cam_code;
prop.seq_id = seq_id;
prop.command_num = command_num;
prop.cdpid_counter = cdpid_counter;
prop.unique_cdpid = unique_cdpid;
prop.product_type = product_type;
prop.gop = gop;
prop.data_proc_code = data_proc_code;
prop.version = vr;

end