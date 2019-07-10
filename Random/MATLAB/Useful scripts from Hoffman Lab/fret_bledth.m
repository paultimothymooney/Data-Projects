function [abt,dbt] = fret_bledth(doa,dod,dof,varargin)
% files = fret_bledth(doaf,dodf,dofr,param)
% files = fret_bledth(doaf,dodf,dofr,baf,bdf,brd,param)
% files = fret_bledth(doaf,dodf,dofr,baf,bdf,brd,aoaf,aodf,aofr,param)

% PURPOSE: Calculate the bleed-throughs and cross-talk associated with FRET
% imaging. Can do linear or non-linear bleedthroughs and cross-talks.
% Bleed-throughs are defined as either fluorophore bleeding into the FRET
% channel and cross-talk is defined as one fluorophore crossing over into
% another's channel. Linear bleed-throughs assume that these percentages
% are constant as a function of brightness. Non-linear bleed-throughs do
% not. In this program, images of single labeled samples and backgrounds
% are read in and the percentages that the fluorophore emit in each channel
% are calculated. Files representing the non-linear bleed-throughs are
% written out and linear estimates are printed to the screen.

% Background subtraction is handled one of two ways. If background images
% are included then these images are averaged and subtracted from the
% single and double labeled images.

%--------------------------------------------------------------------------

% Created 9/5/12 by Katheryn Rothenberg
% Updated 9/12/12 by Wes Maloney - Translated original body of lister
%       subfunction
% Updated 9/14/12 by Wes Maloney - Translated original body of smooth_table
%       and fit_samp subfunctions
% Updated 9/14/12 by Katheryn Rothenberg - Translated original body of main
%       function with adjustments to file loading
% Updated 9/20/12 by Katheryn Rothenberg - Debugged and functioning
% Updated 11/6/12 by Katheryn Rothenberg - Fixed a bug in the fit_samp
%       function with calculation of the mean while binning the data
% Updated 11/9/12 by Katheryn Rothenberg - Adjusted image reading, output
%       plot matches sample data plot, any further changes will be a new
%       version

%--------------------------------------------------------------------------

% INPUTS: Filenames must be passed in either a group of 3, 6, or 9. The
% last input will be a structure with parameters.

% If you input three filenames, the program will calculate the bleedthrough
% assuming it is a set of donor only images. If it is really a set of
% acceptor images, switch the position of the donor and acceptor channel.
% doaf - base name for donor only acceptor channel image files
% dodf - base name for donor only donor channel image files
% dofr - base name for donor only FRET channel image files

% If you input six filenames, the program assumes you have given it a three
% channel set of donor only and background images and will calculate the
% bleedthrough of the donor images with background flat-fielding. Use the
% no background keyword if images are actually a set of three channel
% images for the donor and acceptor
% baf - base name for background images from the acceptor channel
% bdf - base name for background images from the donor channel
% brd - base name for background images from the FRET channel

% If you input nine filenames (RECOMMENDED), the program will calculate
% cross-talks and bleedthroughs for all channels and perform flat-fielding
% calculations.
% aoaf - base name for acceptor only acceptor channel image files
% aodf - base name for acceptor only donor channel image files
% aofr - base name for acceptor only FRET channel image files

% params will be a structure containing a field for each of the following
% parameters:
% sourcefolder - specifies the folder containing all the images being
% analyzed
% destfolder - specifies the destination folder for the output files and
% images
% nobkg - to get naming right if you don't have background images
% bit - specifies bit of the image
% dthres - threshold in the donor only donor channel image (only pixels
% above this intensity are used in the calculation
% athres - threshold in the acceptor only acceptor channel image
% outname - to specify the name of the non-linear correction files
% width - set the number of points over which the non-linear correction
% curves should be smoothed. Large values cause errors at intensity extrema
% ocimg - set to output the images
% npnts - sets the number of points in the smoothed curve
% avg - will average results over block size given, useful if some
% brightness regions are not sufficiently covered
% pf - set to write out files with only prefix names

%--------------------------------------------------------------------------


% initialize/get variables
param = varargin{end};
if ~isfield(param,'athres')
    param.athres = 1;
end
if ~isfield(param,'dthres')
    param.dthres = 1;
end
if ~isfield(param,'bit')
    param.bit = 12;
end
if ~isfield(param,'width')
    param.width = 1500;
end
if ~isfield(param,'npnts')
    param.npnts = 1000;
end
if ~isfield(param,'avg')
    param.avg = 1;
end

param.npnts = double(param.npnts);
daflag = 0;
nps = nargin-1;
param.bin = 1;

doaf = file_search(doa,param.sourcefolder);
dodf = file_search(dod,param.sourcefolder);
dofr = file_search(dof,param.sourcefolder);

if param.nobkgd && nps > 3
    aoaf = file_search(varargin{1},param.sourcefolder);
    aodf = file_search(varargin{2},param.sourcefolder);
    aofr = file_search(varargin{3},param.sourcefolder);
end

% load files
filecell = {'doafn','dodfn','dofrn','bafn','bdfn','bfrn','aoafn','aodfn','aofrn'};
filenames = {doaf,dodf,dofr};
ndoa = length(doaf);
ndod = length(dodf);
ndof = length(dofr);

if param.nobkgd && nps > 3
    filenames{4} = aoaf;
    filenames{5} = aodf;
    filenames{6} = aofr;
    naoa = length(aoaf);
end

for i = 1:3
    for j = 1:ndoa
        eval(sprintf('%s{j} = double(imread(fullfile(''%s'',''%s'')));',...
            filecell{i},param.sourcefolder,filenames{i}{j}));
    end
end


% do background averages
if ~param.nobkgd && nps >= 6
    for i = 1:3 % read in 4-6 as background
        for j = 1:length(varargin{1})
            eval(sprintf('%s{j} = double(imread(''%s''));',filecell{i+3},filenames{i+3}{j}));
        end
    end
    tot = 0;
    tot2 = 0;
    tot3 = 0;
    nba = length(varargin{1});
    nbd = length(varargin{2});
    nbf = length(varargin{3});
    for i=1:nba
        tot = tot+bafn{i};
        tot2 = tot2+bdfn{i};
        tot3 = tot3+bfrn{i};
    end
    axamb = tot/nba;
    dxdmb = tot2/nbd;
    dxamb = tot3/nbf;
elseif nps >= 6 && param.nobkgd
    for i = 1:3 % read in 4-6 as ao
        for j = 1:length(aoaf)
            eval(sprintf('%s{j} = double(imread(fullfile(''%s'',''%s'')));',...
                filecell{i+6},param.sourcefolder,filenames{i+3}{j}))
        end
    end
    daflag = 1;
end

if nps == 9
    for i = 1:3 % read in 7-9 as ao
        for j = 1:length(varargin{4})
            eval(sprintf('%s{j} = double(imread(''%s''));',filecell{i+6},filenames{i+6}{j}))
        end
    end
    daflag = 1;
end

doaxam = cell(1,ndoa);
dodxam = cell(1,ndoa);
dodxdm = cell(1,ndoa);
for i=1:ndoa
    [doaxam{i}, dodxdm{i}, dodxam{i}] = deover(doafn{i},dodfn{i},dofrn{i},param.bit);
    
    if nps >= 6 && ~param.nobkgd
        dodxdm{i} = bs_ff(dodxdm{i},dxdmb,param);
        dodxam{i} = bs_ff(dodxam{i},dxamb,param);
        doaxam{i} = bs_ff(doaxam{i},axamb,param);
    else
        dodxdm{i} = bs_ff(dodxdm{i},param);
        dodxam{i} = bs_ff(dodxam{i},param);
        doaxam{i} = bs_ff(doaxam{i},param);
    end
    
    % get the mask
    wdo = find(dodxdm{i} > param.dthres);
    
    if i ==1
        ddd = lister(dodxdm{i},wdo,i);
        dda = lister(dodxam{i},wdo,i);
        daa = lister(doaxam{i},wdo,i);
    else
        ddd = lister(dodxdm{i},wdo,i,ddd);
        dda = lister(dodxam{i},wdo,i,dda);
        daa = lister(doaxam{i},wdo,i,daa);
    end
    % Structure for tiff writing
    % write files if requested
    if param.ocimg
        if param.nobkgd
            imwrite2tif(dodxdm,[],fullfile(pwd,param.destfolder,['bs' dodfn{i}]),'single')
            imwrite2tif(dodxdm,[],fullfile(pwd,param.destfolder,['bs' doafn{i}]),'single')
            imwrite2tif(dodxdm,[],fullfile(pwd,param.destfolder,['bs' dofrn{i}]),'single')
        else
            imwrite2tif(dodxdm,[],fullfile(pwd,param.destfolder,['bsff' dodfn{i}]),'single')
            imwrite2tif(dodxdm,[],fullfile(pwd,param.destfolder,['bsff' doafn{i}]),'single')
            imwrite2tif(dodxdm,[],fullfile(pwd,param.destfolder,['bsff' dofrn{i}]),'single')
        end
        
    end
end

aoaxam = cell(1,naoa);
aodxam = cell(1,naoa);
aodxdm = cell(1,naoa);
if daflag
    for i = 1:length(aoaxam)
        % do flatfielding and background subtraction
        [aoaxam{i},aodxdm{i},aodxam{i}] = deover(aoafn{i},aodfn{i},aofrn{i},param.bit);
        if nps == 9 && ~param.nobkgd
            aodxdm{i} = bs_ff(aodxdm{i},dxdmb,param);
            aodxam{i} = bs_ff(aodxam{i},dxamb,param);
            aoaxam{i} = bs_ff(aoaxam{i},axamb,param);
        else
            aodxdm{i} = bs_ff(aodxdm{i},param);
            aodxam{i} = bs_ff(aodxam{i},param);
            aoaxam{i} = bs_ff(aoaxam{i},param);
        end
        
        % get the mast
        wao = find(aoaxam{i} > param.athres);
        
        if i ==1
            aaa = lister(aoaxam{i},wao,i);
            ada = lister(aodxam{i},wao,i);
            add = lister(aodxdm{i},wao,i);
        else
            aaa = lister(aoaxam{i},wao,i,aaa);
            ada = lister(aodxam{i},wao,i,ada);
            add = lister(aodxdm{i},wao,i,add);
        end
        
        if param.ocimg
            if param.nobkgd
                imwrite2tif(aodxdm,[],fullfile(pwd,param.destfolder,['bs' aodfn{i}]),'single')
                imwrite2tif(aoaxdm,[],fullfile(pwd,param.destfolder,['bs' aoafn{i}]),'single')
                imwrite2tif(aodxam,[],fullfile(pwd,param.destfolder,['bs' aofrn{i}]),'single')
            else
                imwrite2tif(aodxdm,[],fullfile(pwd,param.destfolder,['bsff' aodfn{i}]),'single')
                imwrite2tif(aoaxdm,[],fullfile(pwd,param.destfolder,['bsff' aodfn{i}]),'single')
                imwrite2tif(aodxam,[],fullfile(pwd,param.destfolder,['bsff' aodfn{i}]),'single')
            end
        end
        
    end
else
    ada = [0 0];
    aaa = [0 0];
    add = [0 0];
end

[ldbt,mddd,mdda] = fit_samp(ddd,dda,param.avg);
[labt,maaa,mada] = fit_samp(aaa,ada,param.avg);
[ldct,mddd,mdaa] = fit_samp(ddd,daa,param.avg);
[lact,maaa,madd] = fit_samp(aaa,add,param.avg);

xmax = max([mddd , maaa]);
plot(mddd,mdda./mddd,'o') % donor bleed-through
axis([0 xmax 0 1.2])
hold on
plot(maaa,mada./maaa,'s') % acceptor bleed-through
plot(mddd,mdaa./mddd,'m^') % donor cross-talk
plot(maaa,madd./maaa,'mx') % acceptor cross-talk

plot([0 2^param.bit],[ldbt(1) ldbt(1)],'r')
plot([0 2^param.bit],[labt(1) labt(1)],'r')
plot([0 2^param.bit],[ldct(1) ldct(1)],'r')
plot([0 2^param.bit],[lact(1) lact(1)],'r')

if max(mddd) > 0 && max(mdda) > 0
    nldbt = smooth_table(mddd,mdda,param.width,param.npnts);
else
    nldbt = zeros(5,2);
end
if max(maaa) > 0 && max(mada) > 0
    nlabt = smooth_table(maaa,mada,param.width,param.npnts);
else
    nlabt = zeros(5,2);
end
if max(mddd) > 0 && max(mdaa) > 0
    nldct = smooth_table(mddd,mdaa,param.width,param.npnts);
else
    nldct = zeros(5,2);
end
if max(maaa) > 0 && max(madd) > 0
    nlact = smooth_table(maaa,madd,param.width,param.npnts);
else
    nlact = zeros(5,2);
end

plot(nldbt(:,1),nldbt(:,2),'g',nlabt(:,1),nlabt(:,2),'g',nldct(:,1),nldct(:,2),'g',...
    nlact(:,1),nlact(:,2),'g')

if nps <= 6 && ~daflag
    fprintf('The bleed-through in the FRET channel is: %.2f\n',ldbt(1))
    fprintf('The cross-talk into the non-FRET channel is: %.2f\n', ldct(1))
end
if daflag
    fprintf('The donor bleed-through into the FRET channel is: %.2f\n',ldbt(1))
    dbt = ldbt(1);
    fprintf('The acceptor bleed-through into the FRET channel is: %.2f\n',labt(1))
    abt = labt(1);
    fprintf('The donor crosstalk into acceptor channel is: %.2f\n',ldct(1))
    fprintf('The acceptor crosstalk into donor channel is: %.2f\n',lact(1))
end

if ~isfield(param,'pf')
    save(fullfile(param.destfolder,'nlabt.dat'),'nlabt','-ascii')
    save(fullfile(param.destfolder,'nldbt.dat'),'nldbt','-ascii')
    save(fullfile(param.destfolder,'nlact.dat'),'nlact','-ascii')
    save(fullfile(param.destfolder,'nldct.dat'),'nldct','-ascii')
else
    if ~isfield(param,'outname')
        param.outname = ' ';
        nele = length(param.pf);
        if nele >= 1
            save(fullfile(pwd,param.destfodler,[param.pf(1) '_' param.outname]),'nldbt','-ascii')
        end
        if nele >= 2
            save(fullfile(pwd,param.destfodler,[param.pf(2) '_' param.outname]),'nldct','-ascii')
        end
        if nele >= 3
            save(fullfile(pwd,param.destfodler,[param.pf(3) '_' param.outname]),'nlabt','-ascii')
        end
        if nele >= 4
            save(fullfile(pwd,param.destfodler,[param.pf(4) '_' param.outname]),'nlact','-ascii')
        end
    end
end

end

%--------------------------------------------------------------------------
% SUBFUNCTIONS

function list = lister(img,w,i,list)
% change the data structure from images type (1024x1024) to list type
% (npixels), also removes zeroes. This means w must be the same in all
% calls

if min(w) > -1
    if i > 1
        list=[list;img(w)];
    else
        list=img(w);
    end
else
    if i == 1
        list=[0,0];
    end
end

end

function [res,xm,ym] = fit_samp(x,y,avg)
% program to optionally bit, fit, and smooth data


if sum(x+y) > 0
    %This is just a quick way to bin by size avg
    x2=round(x./avg);
    [un,s] = unique(x2);
    %     xm = x2(s).*avg;
    %     ym = y(s);
    s = [];
    for i = 1:numel(un)
        s = [s; find(un(i) == x2)];
    end
    
    x2 = x2(s)*avg;
    y2 = y(s);
    [trash,u]=unique(x2);
    nele=length(u);
    xm=zeros(1,nele);
    ym=zeros(1,nele);
    for i=1:nele
        if i == 1
            n=u(1);
            if n >= 3
                xm(i)=mean(x2(1:u(1)));
                ym(i)=mean(y2(1:u(1)));
            end
            if n == 1
                xm(1)=x2(1);
                ym(1)=y2(1);
            end
            if n == 2
                xm(i)=sum(x2(1:2))/double(2);
                ym(i)=sum(y2(1:2))/double(2);
            end
        else
            n=u(i)-(u(i-1)+1);
            xm(i)=mean(x2(u(i-1)+1:u(i)));
            ym(i)=mean(y2(u(i-1)+1:u(i)));
        end
    end
    res=polyfit(xm,ym,1);
else
    res=[0,0];
    xm=res;
    ym=res;
end
end

function res = smooth_table(tx,ty,width,npnts)
% program to smooth data
[x,i]=sort(tx);
y=ty(i);

w=find(x ~= 0 & y ~= 0);
x=x(w);
y=y(w);
f=y./x;

[spx,spf] = spline_p_k(x,f);
sf = smooth(spf,width);

nele=length(sf);
del=(nele/npnts);
if del >= 1
    vec=floor((1:npnts)*del);
    res = [spx(vec) ;sf(vec)']';
else
end
end