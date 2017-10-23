%% Set up
clear;
close all;
clc;
prefix = '';

%% Get Information
folder = input('Type the name of the folder that contains your images, \n make sure it is added to the path, \n and name your files so they look like \n"exp_01_w1channel1.TIF" and "exp_01_w2channel2.TIF": ','s');
% channels = {'TVFRET' 'Venus' 'Teal' 'FWCy5' 'FWDAPI' 'FWFITC' 'FWTR' 'CRFRET'};
channels = {'TVFRET' 'Venus' 'Teal'};
ref_channel = channels{1};
num_channels = length(channels);

% Prealocate space for PreParams
for i = 1:length(channels)
    PreParams.(channels{i}).xshift = [];
    PreParams.(channels{i}).yshift = [];
    PreParams.(channels{i}).k = [];
    PreParams.(channels{i}).ex = [];
    PreParams.(channels{i}).dark = [];
    PreParams.(channels{i}).shade = [];
end

% Prealocate reference channel x, y, k, ex params
k_param(1) = 0.0000000001;
ex_param(1) = 1;
xshift(1) = 0;
yshift(1) = 0;
PreParams.(channels{1}).xshift = xshift(1);
PreParams.(channels{1}).yshift = yshift(1);
PreParams.(channels{1}).k = k_param(1);
PreParams.(channels{1}).ex = ex_param(1);


%% 3D Image Registration using dftregistration.m
if isempty(file_search('\w+_reg_params.txt',folder));
    for i = 2:num_channels
        fretname = file_search(['\w+beads\w+' channels{i} '\w+' ref_channel '.TIF'],folder); % reference image
        othername = file_search(['\w+beads\w+' channels{i} '.TIF'],folder); % going to load stack of comparison images
        fretNumStacks = length(fretname);
        zdata_out = ones(1,7);
        for k = 1:fretNumStacks % number of image stacks
            fretinfo = imfinfo(fretname{k});
            otherinfo = imfinfo(othername{k});
            fretImage{i,k} = single(imread(fretname{k},11,'Info',fretinfo)); % FRET image 11 is "in focus"
            fretImageUsed{i,k} = fft2(fretImage{i,k});
            otherNumImages = length(otherinfo);
            for j = 1:otherNumImages % number of images in each stack
                otherImage{i,j,k} = single(imread(othername{k},j,'Info',otherinfo));
                otherImage1{i,j,k} = fft2(otherImage{i,j,k});
                [zdata{i,k}(j,1:4),~] = dftregistration(fretImageUsed{i,k},otherImage1{i,j,k},100);
                zdata{i,k}(j,5) = j; % Col5 = slice
            end
            % Calculate optimal z-plane
            zdata{i,k}(:,6) = k; % Col6 = stack
            x = zdata{i,k}(:,5); % slice
            y = zdata{i,k}(:,1); % RMSE
            sp = spline(x,y);
            [minval, minsite] = fnmin(sp);
            img_slices(k) = minsite;
            zdata{i,k}(:,7) = (minsite-11)*100; % Col7 = z offset (nm)
            row = round(minsite); % need to know which slice to use in subsequent registration calculations
            zdata_out = vertcat(zdata_out,zdata{i,k}(row,:)); % compile data for reference
            
            % Plot to make sure spline worked
            hold on
            xx = 0:.25:21;
            yy = spline(x,y,xx);
            plot(x,y,'o',xx,yy);
            scatter(minsite,minval);
            saveas(gcf, fullfile(pwd,folder,['Zopt_' channels{i} '_' num2str(k)]), 'jpeg')
            hold off
            close
        end
        zdata_out1 = double(zdata_out);
        zdata_out1(1,:) = [];
        save(fullfile(pwd,folder,[folder '_' ref_channel '_' channels{i} '_reg_params.txt']),'zdata_out1','-ascii');
        
        % Average optimal z-planes, x-shifts, and y-shifts calculated above
        img_slice(i) = round(mean(img_slices));
        xshift(i) = mean(zdata_out1(:,4));
        yshift(i) = mean(zdata_out1(:,3));
        PreParams.(channels{i}).xshift = xshift(i);
        PreParams.(channels{i}).yshift = yshift(i);
    end
else % Load txt file with x and y shifts previously calculated
    for i = 2:num_channels
        XYdataName = file_search(['\w+' channels{i} '\w+reg_params.txt'],folder);
        XYdata = load(XYdataName{1});
        xshift(i) = mean(XYdata(:,4));
        yshift(i) = mean(XYdata(:,3));
        PreParams.(channels{i}).xshift = xshift(i);
        PreParams.(channels{i}).yshift = yshift(i);
    end
end


%% Extract in focus images and register them
rehash
if isempty(file_search('reg_slice\w+.TIF',folder));
    for i = 2:num_channels
        fretname = file_search(['\w+beads\w+' channels{i} '\w+' ref_channel '.TIF'],folder); % reference image
        othername = file_search(['\w+beads\w+' channels{i} '.TIF'],folder); % going to load stack of comparison images
        fretNumStacks = length(fretname);
        slice = img_slice(i);
        for k = 1:fretNumStacks % number of image stacks
            %Load images
            fretinfo = imfinfo(fretname{k});
            otherinfo = imfinfo(othername{k});
            fretImage{i,k} = single(imread(fretname{k},11,'Info',fretinfo)); % FRET image 11 is "in focus"
            otherImage{i,k} = single(imread(othername{k},slice,'Info',otherinfo)); %Extract appropriate other images
            % Register and crop to eliminate edge pixels
            sz = size(otherImage{i,k});
            [y,x] = ndgrid(1:sz(1),1:sz(2));
            otherImage{i,k} = interp2(x,y,otherImage{i,k},x-xshift(i),y-yshift(i)); % register the other image to the FRET channel
            crop = [round(0.0246*sz(1)) round(0.0246*sz(1)) round(0.9509*sz(1)) round(0.9509*sz(1))];
            fretImage{i,k} = imcrop(fretImage{i,k},crop); % Crop FRET image
            otherImage{i,k} = imcrop(otherImage{i,k},crop); % Crop other image
            %Write to tif
            imwrite2tif(fretImage{i,k},[],fullfile(pwd,folder,['reg_slice_11_' fretname{k}]),'uint16');
            imwrite2tif(otherImage{i,k},[],fullfile(pwd,folder,['reg_slice_' num2str(slice) '_' othername{k}]),'uint16');
        end
    end
end

%% Particle Tracking to calculate radial distortion correction (k and ex params)
if isempty(file_search('\w+rad_params.txt',folder));
    rehash
    for i = 2:num_channels
        nameI1 = file_search(['reg_slice\w+' channels{i} '\w+' ref_channel '.TIF'],folder);
        nameI2 = file_search(['reg_slice\w+' channels{i} '.TIF'],folder);
        k_out = ones(1,2);
        for k = 1:length(nameI1)
            I1 = imread(nameI1{k});
            I2 = imread(nameI2{k});
            lambda = 1; %length scale of noise to be filtered out; typically 1 pixel
            w = 9; % should be a little greater than the radius of the largest features
            f1 = feature2D(I1,lambda,w);
            f2 = feature2D(I2,lambda,w);
            f1(:,6) = 1;
            f2(:,6) = 2;
            f1(:,7) = 1;
            f2(:,7) = 2;
            out = vertcat(f1,f2);
            
            % Particle Tracker
            [lub] = trackmem(out,5,2,2,0);
            x1 = lub(1:2:end,1);
            y1 = lub(1:2:end,2);
            x2 = lub(2:2:end,1);
            y2 = lub(2:2:end,2);
            
            % Make sure particle tracker isn't favoring particles at the center or
            % edges of pixels
            figure
            subplot(2,1,1)
            hist(mod(x1,1));
            subplot(2,1,2)
            hist(mod(x2,1));
            saveas(gcf, fullfile(pwd,folder,['mod_' channels{i} '_' num2str(k)]), 'jpeg')
            close
            
            % Visualize points connected by lines
            hold on
            for j = 1:length(x1)
                plot([x1(j) x2(j)], [y1(j) y2(j)], 'k-')
            end
            saveas(gcf, fullfile(pwd,folder,['quiver_' channels{i} '_' num2str(k)]), 'jpeg')
            hold off
            close
            
            % Shift x1, y1, x2, y2 to center, convert to polar coordinates and
            % normalize to R
            center = round(2039/2);
            x1 = x1-center;
            x2 = x2-center;
            y1 = y1-center;
            y2 = y2-center;
            [theta1,r1] = cart2pol(x1,y1);
            [theta2,r2] = cart2pol(x2,y2);
            R = sqrt(2)*center;
            r1 = r1./R;
            r2 = r2./R;
            w = find(r1>0.8);
            r1(w)=[];
            r2(w)=[];
            
            %Scatter plot to fit curves to where
            % y = s-x = TVFRETchannel location - FWCy5 channel location (r1-r2)
            % x = r = location of distorted FWCy5 channel (r2)
            scatter(r2,(r1-r2));
            axis([0 0.8 -.002 0.002]);
            saveas(gcf, fullfile(pwd,folder,['points_' channels{i} '_' num2str(k)]), 'jpeg')
            close
            
            % Use non linear tool to find parameters for curve fit
            x_in = r2;
            y_in = r1-r2;
            F_type = 5;
            F = @(k,r)k(1).*r.^k(2); %Equation used to radially distort images
            k0 = [1 1];
            [kp,resnorm,~,exitflag,output] = lsqcurvefit(F,k0,x_in,y_in); % use different value for k
            hold on
            scatter(x_in,y_in);
            scatter(x_in,F(kp,x_in));
            axis([0 0.8 -.002 0.002]);
            saveas(gcf, fullfile(pwd,folder,['radial_' channels{i} '_' num2str(k)]), 'jpeg')
            hold off
            close
            kp(:,1) = -kp(:,1); %Must switch sign of k term to correct radial distortion
            k_out = vertcat(k_out,kp);
        end
        k_out(1,:)=[];
        k_param(i) = mean(k_out(:,1));
        ex_param(i) = mean(k_out(:,2));
        PreParams.(channels{i}).k = k_param(i);
        PreParams.(channels{i}).ex = ex_param(i);
        save(fullfile(folder,[folder '_' ref_channel '_' channels{i} '_rad_params.txt']),'-ascii','k_out');
    end
else % Load txt file with k and ex params previously calculated
    for i = 2:num_channels
        RADdataName = file_search(['\w+' channels{i} '\w+rad_params.txt'],folder);
        RADdata = load(RADdataName{1});
        k_param(i) = mean(RADdata(:,1));
        ex_param(i) = mean(RADdata(:,2));
        PreParams.(channels{i}).k = k_param(i);
        PreParams.(channels{i}).ex = ex_param(i);
    end
end

%% Extract in focus images and radially correct them (just to output images and make sure they are preprocessed correctly
rehash
if isempty(file_search('rad_reg_slice\w+.TIF',folder));
    for i = 2:num_channels
        fretname = file_search(['reg_slice\w+beads\w+' channels{i} '\w+' ref_channel '.TIF'],folder); % reference image
        othername = file_search(['reg_slice\w+beads\w+' channels{i} '.TIF'],folder); % going to load stack of comparison images
        fretNumImgs = length(fretname);
        for k = 1:fretNumImgs % number of image stacks
            %Load images
            fretImage{i,k} = single(imread(fretname{k}));
            otherImage{i,k} = single(imread(othername{k}));
            otherImage{i,k} = lensdistort(otherImage{i,k},k_param(i),ex_param(i),'ftype',5,'bordertype','fit');
            %Write to 32bit .TIF
            imwrite2tif(fretImage{i,k},[],fullfile(pwd,folder,['rad_' fretname{k}]),'uint16');
            imwrite2tif(otherImage{i,k},[],fullfile(pwd,folder,['rad_' othername{k}]),'uint16');
        end
    end
end

%% Load dark images and calculate dark
rehash
if isempty(file_search('dark_mean\w+.TIF',folder));
    for i = 1:num_channels
        darkName = file_search(['dark_\w+' channels{i} '.TIF'],folder);
        a = imread(darkName{1});
        [r,c] = size(a);
        dark = zeros(r,c);
        for k = 1:length(darkName)
            dark = dark + single(imread(darkName{k}));
        end
        dark_mean{i} = dark./length(darkName);
        PreParams.(channels{i}).dark = dark_mean{i};
        imwrite2tif(dark_mean{i},[],fullfile(pwd,folder,['dark_mean_' channels{i} '.TIF']),'uint16');
    end
else
    for i = 1:num_channels
        dark_mean_Name = file_search(['dark_mean\w+' channels{i} '.TIF'],folder);
        dark_mean{i} = single(imread(dark_mean_Name{1}));
        PreParams.(channels{i}).dark = dark_mean{i};        
    end
end

%% Load shade images and correct them like the other images thusfar
rehash
if isempty(file_search('pre_shade\w+.TIF',folder));
    for i = 1:num_channels
        shadeImgNames = file_search(['shade\w+' channels{i} '.TIF'],folder);
        % Correct them and save them out
        for k = 1:length(shadeImgNames)
            img = single(imread(shadeImgNames{k}));
            img = img - dark_mean{i};
            sz = size(img);
            [y,x] = ndgrid(1:sz(1),1:sz(2));
            img = interp2(x,y,img,x-xshift(i),y-yshift(i));
            img = lensdistort(img,k_param(i),ex_param(i),'ftype',5,'bordertype','fit');
            sz = size(img);
            crop = [round(0.0246*sz(1)) round(0.0246*sz(1)) round(0.9509*sz(1)) round(0.9509*sz(1))];
            img = imcrop(img,crop);
            imwrite2tif(img,[],fullfile(folder,['pre_' shadeImgNames{k}]),'single');
        end
    end
end

%% Generate shade correct bnorm with corrected images
rehash
if isempty(file_search('bnorm\w+.TIF',folder))
    for i = 1:num_channels
        shade_correct_gen(['pre_shade\w+' channels{i} '.TIF'],folder)
        rehash
        bnorm_name = file_search(['bnorm\w+' channels{i}],folder);
        bnorm = single(imread(bnorm_name{1}));
        PreParams.(channels{i}).shade = bnorm;
    end
else
    for i = 1:num_channels
        bnorm_name = file_search(['bnorm\w+' channels{i}],folder);
        bnorm = single(imread(bnorm_name{1}));
        PreParams.(channels{i}).shade = bnorm;
    end
end

%% Save out PreParams into folder
save(fullfile(pwd,folder, 'PreParams.mat'),'-struct','PreParams');