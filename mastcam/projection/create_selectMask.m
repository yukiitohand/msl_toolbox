function [selectMask,srange,lrange] = create_selectMask(x,y)
% [selectMask,srange,lrange] = create_selectMask(x,y)

x_min = min(x);
x_max = max(x);
y_min = min(y);
y_max = max(y);
samplesc = x_max - x_min + 1;
linesc   = y_max - y_min + 1;
sample_offset = x_min-1;
line_offset   = y_min-1;
srange = [x_min x_max];
lrange = [y_min y_max];

x_ofst = x - sample_offset;
y_ofst = y - line_offset;

selectMask = zeros(linesc,samplesc,'uint8');
Npts = length(x_ofst);
for ni=1:Npts
    selectMask(y_ofst(ni),x_ofst(ni)) = 1;
end

if linesc==1 && samplesc==1
    selectMask = zeros([3,3],'uint8');
    selectMask(2,2) = 1;
    srange = [x_min-1 x_min x_min+1];
    lrange = [y_min-1 y_min y_min+1];
elseif linesc==1
    selectMask_pd = zeros([3,samplesc],'uint8');
    selectMask_pd(2,:) = selectMask;
    selectMask = selectMask_pd;
    lrange = [y_min-1 y_min y_min+1];
elseif samplesc==1
    selectMask_pd = zeros([linesc,3],'uint8');
    selectMask_pd(:,2) = selectMask;
    selectMask = selectMask_pd;
    srange = [x_min-1 x_min x_min+1];
end

end