function [ data,colinfo2,colinfo_names ] = msl_telemetryCSVread(fpath,lbl_telemetry)
% [ data,colinfo2,colinfo_names ] = msl_telemetryCSVread(fpath,lbl_telemetry)
%   Read MSL telemetry csv file
%    Input Parameters
%    fpath: file path to the csv file '*.csv'
%    lbl_telemtry: pds lbl file of the telemetry file
%    Output Parameters
%     data: struct, field names are defined in 'NAME' in each of struct in 
%           lbl_telemetry.OBJECT_SPREADSHEET.OBJECT_FIELD
%     colinfo: struct, field are same as the each element in
%              lbl_telemetry.OBJECT_SPREADSHEET.OBJECT_FIELD
%     colinfo_names: struct, field are same as the each element in
%              lbl_telemetry.OBJECT_SPREADSHEET.OBJECT_FIELD, field names are colinfo(i).NAME for
%              easy access with their names.


obj_csv = lbl_telemetry.OBJECT_SPREADSHEET;

nCols = obj_csv.FIELDS;
nRows = obj_csv.ROWS;
colinfo = obj_csv.OBJECT_FIELD;

if length(colinfo)==1
    colinfo = {colinfo};
end

nameList = cell(1,length(colinfo));

[name] = mod_fieldname(colinfo{1}.NAME);
nameList{1} = name;
data = struct(name,cell(nRows,1));

for c=2:nCols
    [name] = mod_fieldname(colinfo{c}.NAME);
    nameList{c} = name;
    [data.(name)] = deal([]);
end

fp = fopen(fpath,'r');
% fseek(fp,skip_byte,0);
fgets(fp);

% create format spec
formatSpec = '';
for c=1:nCols
    if c>1
        formatSpec = [formatSpec ''];
    end
    switch colinfo{c}.DATA_TYPE
        case 'CHARACTER'
            formatSpec_c = ['%s'];
            formatSpec = [formatSpec formatSpec_c];
        case 'ASCII_REAL'
            formatSpec_c = lower(['%f']);
            formatSpec = [formatSpec formatSpec_c];

        case {'INTEGER','ASCII_INTEGER'}
            formatSpec_c = ['%d'];
            formatSpec = [formatSpec formatSpec_c];
        otherwise
            error('c=%d,DATA_TYPE %s is not defined',c,colinfo{c}.DATA_TYPE);
    end
end

data_raw = textscan(fp,formatSpec,nRows,'delimiter',',');
%for i=1:size(data,2)
%    spc(:,i) = data{i};
%end
for c=1:nCols
    switch colinfo{c}.DATA_TYPE
        case 'CHARACTER'
            data_rawc = data_raw{c};
        case {'ASCII_REAL','INTEGER','ASCII_INTEGER'}
            data_rawc = num2cell(data_raw{c});
        otherwise
            error('c=%d,DATA_TYPE %s is not defined',c,colinfo{c}.DATA_TYPE);
    end
    [data.(nameList{c})] = data_rawc{:};
end

fclose(fp);


if length(colinfo)==1
    colinfo2 = colinfo;
else
    colinfo2 = merge_struct(colinfo{:});
end

colinfo_names = [];
for i=1:length(colinfo2)
    colinfo_names.(nameList{i}) = colinfo2(i);
end

end

function [name_m] = mod_fieldname(name)
if ~isempty(regexpi(name,'^[\d]+.*','ONCE'))
    name = ['COLNAME_' name];
end
name = replace(name,{',',' ',':','(',')'},'_');
name_m = replace(name,{';','^','/'},'');
end