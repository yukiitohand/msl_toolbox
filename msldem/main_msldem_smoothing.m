MSLDEMdata = HSI('MSL_Gale_DEM_Mosaic_1m_v3','/Users/yukiitoh/data/');
wndw_size = 5;
tic; msldem_numeqnb = msldem_get_numeqnb_mex(MSLDEMdata.imgpath,MSLDEMdata.hdr,wndw_size); toc;
% remove the pixels that cannot be traced back to the pixel that has the
% value larger than 4.
eval_num = 4; eval_num_min = 2;
tic; msldem_numeqnbcl = msldem_cleanup_numeqnb_step1_mex(MSLDEMdata.hdr,msldem_numeqnb,eval_num,eval_num_min); toc;
% some pixels that cannot be detected but has the same value with low
% resolution pixels are compensated.
tic; msldem_numeqnbcl2 = msldem_cleanup_numeqnb_step2_mex(...
    MSLDEMdata.imgpath,MSLDEMdata.hdr,msldem_numeqnb,msldem_numeqnbcl,5); toc;

% remove pixels that are not directly connected to the pixel with the value
% larger than 7.
tic; msldem_numeqnbcl2p = msldem_cleanup_numeqnb_step2p_mex(...
    MSLDEMdata.hdr,msldem_numeqnbcl2,7); toc;

% detect low resolution pixel with the value of 1.
tic; msldem_numeqnbcl3 = msldem_cleanup_numeqnb_step3_mex(...
    MSLDEMdata.imgpath,MSLDEMdata.hdr,msldem_numeqnb,msldem_numeqnbcl2p,5,1); toc;

rc_exclude = [...
    7944 26658;...
    7945 26658];

rc_include = [...
    39041 19780;...
    39042 19780;...
    19379 7121;...
    19380 7121];

for i=1:size(rc_exclude,1)
    msldem_numeqnbcl3(rc_exclude(i,1),rc_exclude(i,2)) = 0;
end
for i=1:size(rc_include,1)
    msldem_numeqnbcl3(rc_include(i,1),rc_include(i,2)) = 1;
end

tic; msldem_numeqnbcl3p = msldem_cleanup_numeqnb_step3p_mex(...
    MSLDEMdata.imgpath,MSLDEMdata.hdr,msldem_numeqnb,msldem_numeqnbcl3,5,0); toc;

% correspond to 6 12 15 19 20 25 26 28 30 36 37?
rc_exclude_0 = [...
       19455        7132;...
       18574       13102;...
       13444       15253;...
       14148       18541;...
       46425       19262;...
       39870       19884;...
       33522       20089;...
        9684       23497;...
       10237       23570;...
       11905       25386;...
       34677       29284];
   
for i=1:size(rc_exclude_0,1)
    msldem_numeqnbcl3p(rc_exclude_0(i,1),rc_exclude_0(i,2)) = 0;
end


%%
%
% help detecting false positives in the last step. by changing the last
% argument with 2 and 3, I created the following list of pixels that are
% falsely detected as low resolution pixels.
% tic; msldem_numeqnbcl3pp = msldem_cleanup_numeqnb_step3pp_mex(...
%     MSLDEMdata.imgpath,MSLDEMdata.hdr,msldem_numeqnb,msldem_numeqnbcl3p,21,2); toc;

rc_exclude_3pp = [...
       33558        4952;...
       33559        4952;...
       33560        4952;...
       30363        4547;...
       30364        4547;...
       30363        4548;...
       31874        4735;...
       31874        4736;...
       31874        4737;...
       34023        4744;...
       34024        4744;...
       34025        4744;...
       22524        5504;...
       22524        5506;...
       22524        5507;...
       17672        6919;...
       17673        6919;...
       17673        6920;...
       17821        6936;...
       17822        6936;...
       17822        6937;...
       22260        7465;...
       22260        7466;...
       22260        7467;...
       16452        9601;...
       16452        9603;...
       16454        9603;...
       16302       10232;...
       16303       10232;...
       16302       10233;...
       16402       12323;...
       16401       12324;...
       16402       12324;...
       19981       12907;...
       19982       12907;...
       19981       12909;...
       17208       12937;...
       17209       12937;...
       17210       12937;...
       13403       16023;...
       13405       16023;...
       13403       16024;...
       14232       17991;...
       14234       17991;...
       14236       17991;...
        7819       18235;...
        7819       18236;...
        7820       18236;...
       14145       18541;...
       14146       18541;...
       14147       18541;...
       43193       20155;...
       43195       20155;...
       43196       20155;...
        7620       20277;...
        7620       20278;...
        7621       20278;...
        9072       23417;...
        9071       23418;...
        9072       23418];

for i=1:size(rc_exclude_3pp,1)
    msldem_numeqnbcl3p(rc_exclude_3pp(i,1),rc_exclude_3pp(i,2)) = 0;
end



%%

bname_out = [MSLDEMdata.basename '_dave'];
fpath_out_img = joinPath(MSLDEMdata.dirpath,[bname_out '.img']);
msldem_average_20x20_lowrespixels_double_direct_mex(MSLDEMdata.imgpath,...
    MSLDEMdata.hdr,fpath_out_img,msldem_numeqnbcl3p);
% 
fpath_out_hdr = joinPath(MSLDEMdata.dirpath,[bname_out '.hdr']);
copyfile(MSLDEMdata.hdrpath,fpath_out_hdr);

%%
% x = round(ax.XLim);
% y = round(ax.YLim);

%%

%for i=12:length(row)
% i=22;
%rowc = row(i); colc = col(i);
% rowc = 23460; colc = 3681;
% rowc = 13681; colc = 13302;
% rowc = 25937; colc = 3996;
% rowc = 7944; colc = 26658;
% rowc = 19379; colc = 7121;
% rowc = 39042; colc = 19780;
rowc = 14148; colc =18541;
%rowc = 25300; colc = 3908;
%rowc = 23148; colc = 3641;
y = [(rowc-300),(rowc+300)];
x = [(colc-300),(colc+300)];
%
tic; msldemc_img = msldem_lazyenvireadRect(MSLDEMdata,...
    x(1)-1,y(1)-1,x(2)-x(1)+1,y(2)-y(1)+1,'precision','single'); toc; 

MSLDEMdata_ave = HSI('MSL_Gale_DEM_Mosaic_1m_v3_ave','/Users/yukiitoh/data/');
MSLDEMdata_dave = HSI('MSL_Gale_DEM_Mosaic_1m_v3_dave','/Users/yukiitoh/data/');

tic; msldemcave_img = msldem_lazyenvireadRect(MSLDEMdata_ave,...
    x(1)-1,y(1)-1,x(2)-x(1)+1,y(2)-y(1)+1,'precision','single'); toc; 
tic; msldemcdave_img = msldem_lazyenvireadRect(MSLDEMdata_dave,...
    x(1)-1,y(1)-1,x(2)-x(1)+1,y(2)-y(1)+1,'precision','single'); toc; 

%
a = ImageStackView({},'Ydir','reverse');
a.add_layer(x,y,msldemc_img,'name','dem');
a.Update_ImageAxes_LimHomeAuto();
a.Update_ImageAxes_LimHome();
a.Update_axim_aspectR();
a.Restore_ImageAxes2LimHome();
% a.add_layer(x,y,msldemcave_img,'name','dem_ave');
a.add_layer(x,y,msldemcdave_img,'name','dem_dave');
% a.add_layer(x,y,msldem_numeqnb(y(1):y(2),x(1):x(2)),'name','idx');
% a.add_layer(x,y,msldem_numeqnbcl(y(1):y(2),x(1):x(2)),'name','idxc');
% a.add_layer(x,y,msldem_numeqnbcl2(y(1):y(2),x(1):x(2)),'name','idxc2');
% a.add_layer(x,y,msldem_numeqnbcl3(y(1):y(2),x(1):x(2)),'name','idxc3');
% a.add_layer(x,y,msldem_numeqnbcl3p(y(1):y(2),x(1):x(2)),'name','idxc3p');
% a.add_layer(x,y,msldem_numeqnbcl3pp(y(1):y(2),x(1):x(2)),'name','idxc3pp');
% a.add_layer(x,y,msldem_numeqnbcl3p(y(1):y(2),x(1):x(2)),'name','idxc3p');
% a.add_layer(x,y,msldem_numeqnbcl4(y(1):y(2),x(1):x(2)),'name','idxc4');
% a.add_layer(x,y,msldem_numeqnbcl4(y(1):y(2),x(1):x(2)),'name','idxc4');
%pause;
%delete(a);
%end