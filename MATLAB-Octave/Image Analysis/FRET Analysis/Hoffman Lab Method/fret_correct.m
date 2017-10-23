function fret_correct(aexp,dexp,fexp,abtfn,dbtfn,varargin)
% outputs = fret_correct(af,df,fr,abtfn,dbtfn,params)
% outputs = fret_correct(af,df,fr,abtfn,dbtfn,baf,bdf,bfr,params)

% PURPOSE: A program to remove the intensity in a set of FRET images due to
% bleed-throughs or cross-talks. Can do linear or non-linear bleedthroughs
% and cross-talks. Linear bleed-throughs assume that these percentages are
% constant as a function of brightness. Non-linear bleed-throughs do not.
% If all four corrections are used, the corrections scheme becomes a
% non-linear set of coupled equations. In this case, Newton's method is
% used to solve the set of equations. This is not fully validated and
% probably should have some sort of regularization to help enure that the
% solution is a global minimum and not a local one

% Non-linear corrections are done by linearly interpolating the correction
% curves.

% Background subtraction is handled one of two ways. If background levels
% are included, then these images are averaged and subtracted from the
% single and double labeled images. These are called bsff images. If
% background images are not included, then the background (most likely
% pixel in the picture) is subtracted from the image. This will not account
% for variations in illumination. These images are named bs.

%--------------------------------------------------------------------------

% Created 9/5/12 by Katheryn Rothenberg
% Updated 9/14/12 by Wes Maloney - Translated original body of all_cor,
%       all_cor_func, deunder, norm_fret, align_images subfunctions
% Updated 9/25/12 by Katheryn Rothenberg - Translated original body of main
%       function. Debugged and running except for final .tif writing.
% Updated 11/27/12 by Katheryn Rothenberg - Edited the method of saving
%       .tif images from imwrite to using the Tiff class and fixed the
%       check for linear or nonlinear correction
% Updated 11/29/12 by Katheryn Rothenberg - Removed the spline_intrp
%       subfunction and replaced with a call to the cubic spline function.
%       This resulted in final images with very little variation from the
%       test images while using the test bleed through files. Also, added
%       the sourcefolder and dest folder fields to the params structure.
% Updated 06/17/14 by Katheryn Rothenberg - cut out all non-linear
%       calculations and considerations of cross-talk and dealing with
%       background images.

%--------------------------------------------------------------------------

% INPUTS:

% inputs that are required to run the function are as follows:
% af - base name for the acceptor channel images of a double labeled sample
% df - base name for the donor channel images of a double labeled sample
% fr - base name to find the FRET channel images of a double labeled sample
% abtfn - specifies acceptor bleed-through into the FRET channel, either
% number or file
% dbtfn - specifies donor bleed-through into the FRET channel, either
% number or file

% params will be a structure containing a field for each of the following
% parameters:
% sourcefolder - specifies the folder containing all the images being
% analyzed
% destfolder - specifies the destination folder for the output files and
% images
% ocimg - if set, the program will write out the images
% bit - specifies the bit of the image
% imin - the minimum intensity to allow into the output images, if either
% acceptor or donor image is below this intensity, it will be zeroed in
% output images
% outname - manually specifies output name
% double_norm - normalize to the intensity of the corrected acceptor and
% donor images
% donor_norm - normalize to the intensity of the corrected donor images

%--------------------------------------------------------------------------

% OUTPUTS:

% If ocimg is set, program writes out background subtracted version of the
% data images as well as a corrected FRET, and a normalized FRET image.
% Corrected FRET images are labeled c. and the normalized images are
% normalized to the acceptor and entitled cna.

%--------------------------------------------------------------------------

param = varargin{end};
if ~isfield(param,'bit')
    param.bit = 16;
end

param.bin = 1.0;

af = file_search(aexp,param.sourcefolder);
df = file_search(dexp,param.sourcefolder);
fr = file_search(fexp,param.sourcefolder);

% Check images found
if isempty(af)
    warning(['No Acceptor Images found, searched using: ',aexp])
end
if isempty(df)
    warning(['No Donor Images, searched using: ',dexp])
end
if isempty(fr)
    warning(['No FRET Images, searched using: ', fexp])
end

% Check number of images
if length(af)~=length(df)
    warning('Number of Donor and Acceptor Images Not the Same')
end
if length(af)~=length(fr)
    warning('Number of FRET and Acceptor Images Not the Same')
end
if length(df)~=length(fr)
    warning('Number of Donor and FRET Images Not the Same')
end

abt = double(abtfn{1});
dbt = double(dbtfn{1});

for i = 1:length(df)
    % Get fluorescent images
    axam = double(imread(fullfile(param.sourcefolder,af{i})));
    dxdm = double(imread(fullfile(param.sourcefolder,df{i})));
    dxam = double(imread(fullfile(param.sourcefolder,fr{i})));
    
    [axam, dxdm, dxam] = deover(axam,dxdm,dxam,param.bit);
    
    if isfield(param,'imin')
        [axam, dxdm, dxam] = deunder(axam,dxdm,dxam,param.imin);
    end
    
    % Following supplemental material of Chen, Puhl et al BJ Lett 2006
    iaa = axam;
    idd = dxdm;
    fc = dxam-abt.*iaa-dbt.*idd;
    if ~isfield(param,'leave_neg')
        iaa(iaa<0) = 0;
        idd(idd<0) = 0;
        fc(fc<0) = 0;
    end
    
    nafc = norm_fret(fc,iaa);
    if isfield(param,'outname')
        nfrn = ['_' param.outname fr{i}];
        ndfn = ['_' param.outname df{i}];
        nafn = ['_' param.outname af{i}];
    else
        nfrn = ['_' fr{i}];
        ndfn = ['_' df{i}];
        nafn = ['_' af{i}];
    end
    
    target_folder = fileparts([pwd '/' param.destfolder '/c' nfrn]);
    if (not(exist(target_folder,'dir')))
        mkdir(target_folder);
    end
    imwrite2tif(fc,[],fullfile(pwd,param.destfolder,['c' nfrn]),'single')
    
    target_folder = fileparts([pwd '/' param.destfolder '/cna' nfrn]);
    if (not(exist(target_folder,'dir')))
        mkdir(target_folder);
    end
    imwrite2tif(nafc,[],fullfile(pwd,param.destfolder,['cna' nfrn]),'single')
    
    if param.donor_norm
        ndfc = norm_fret(fc,idd);
        imwrite2tif(ndfc,fullfile(pwd,param.destfolder,['cnd' nfrn]),'single')
    end
    if param.double_norm
        nandfc = norm_fret(fc,iaa,idd);
        imwrite2tif(nandfc,fullfile(pwd,param.destfolder,['cnand' nfrn]),'single')
    end
    
    if param.ocimg
        target_folder = fileparts([pwd '/' param.destfolder '/bsd' ndfn]);
        if (not(exist(target_folder,'dir')))
            mkdir(target_folder);
        end
        imwrite2tif(idd,[],fullfile(pwd,param.destfolder,['bsd' ndfn]),'single')
        
        target_folder = fileparts([pwd '/' param.destfolder '/bsa' nafn]);
        if (not(exist(target_folder,'dir')))
            mkdir(target_folder);
        end
        imwrite2tif(iaa,[],fullfile(pwd,param.destfolder,['bsa' nafn]),'single')
    end
    
end

end

%--------------------------------------------------------------------------
% SUBFUNCTIONS:

function ftemp = norm_fret(f,n,varargin)

n_params = nargin;
if n_params== 2
    ftemp=f;
    ntemp=n;
    wnz=find(n ~= 0);
    wz =find(n == 0);
    if ~isempty(wnz)
        ftemp(wnz)=ftemp(wnz)./ntemp(wnz);
    end
    if ~isempty(wz)
        ftemp(wz)=0;
    end
elseif n_params == 3
    ftemp=f;
    ntemp=n;
    ntemp2=n2;
    wnz=find(n ~= 0 & n2 ~= 0);
    ftemp(wnz)=ftemp(wnz)./(ntemp(wnz).*ntemp2(wnz));
    ftemp(n == 0)=0;
    ftemp(n2 == 0)=0;
end
end

function [a,b,c] = deunder(a,b,c,thres)
% Identify pixels less than a certain value. If any of the pixels in a
% single image is less than the threshold, the pixel will be set to zero in
% all images

w=find(a < thres(1) | b < thres(2) | c < thres(3));
if ~isempty(w)
    a(w)=0;
    b(w)=0;
    c(w)=0;
end
end
