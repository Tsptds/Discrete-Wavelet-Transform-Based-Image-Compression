clear all
clc
close all
% Test & Default / Alternate values are commented out

% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<SETTINGS>>>>>>>>>>>>>>>>>>>>>>>>>

% Enable/Disable rescaling - 2=wcodemat 1=rescale 0=disable
global rescaleWhenPlotting;rescaleWhenPlotting=2

% 1 to 12 inclusive, Approximations Excepted(12+4)
MaxEnergyBandsCount=12;

% NUMBER OF FILTER TAPS - 6, 8, 10
[h0, h1, g0, g1]=setTaps(10);

% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

% Iteration Lists
global It_1 It_2 It_3 It_4; % 4 Photos Each (16)

% Energy & SNR List
global En_list; 




% Load image
% image = imread('lena.jpg');
[file,path] = uigetfile({'*.jpg;*.jpeg'});
if isequal(file,0)
   disp('User selected Cancel');
%    clear all
   return
else
   disp(['User selected ', fullfile(path,file)]);
   image = imread(fullfile(path,file));
end



if size(image, 3) == 3
    image = rgb2gray(image);
end

imagefirst=image; % To Check Later, image will be iterated


%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<---------- DECOMPOSE 1-----------
It_1=decompose(image,h0,h1,It_1);

showDecomp(It_1);
sgtitle("Subbands of Iteration 1")
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<---------- DECOMPOSE 2-----------
image=It_1{1};
It_2=decompose(image,h0,h1,It_2);

showDecomp(It_2);
sgtitle("Subbands of Iteration 2")
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<---------- DECOMPOSE 3-----------
image=It_2{1};
It_3=decompose(image,h0,h1,It_3);

showDecomp(It_3);
sgtitle("Subbands of Iteration 3")
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<---------- DECOMPOSE 4-----------
image=It_3{1};
It_4=decompose(image,h0,h1,It_4);

showDecomp(It_4);
sgtitle("Subbands of Iteration 4")

% rec = idwt2(It_1{1}, It_1{3}, It_1{2}, It_1{4}, g0_6tap, g1_6tap); % LH and HL are reversed for dwt2
%     figure
%     showImage(rec)
%     imshow(imagefirst-uint8(rec));

%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<------FIND ENERGY BANDS---------
showNRG(It_1);
showNRG(It_2);
showNRG(It_3);
showNRG(It_4);

maxBands=maximums(En_list,MaxEnergyBandsCount);
fprintf("\n Max "+MaxEnergyBandsCount+" Energy Bands (Except Approximations, 1 to 16, Each 4 being 1 Iteration): \n")
disp(sort(maxBands(:)',"descend"))

%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<------SELECT BANDS TO RECONSTRUCT---
%                           Defaults are commented out above the rec lines

%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<---------- RECOMPOSE 4-----------
% rec = idwt2(It_4{1}, It_4{3}, It_4{2}, It_4{4}, g0, g1);
rec = idwt2(It_4{1}, It_4{3}, It_4{2}, It_4{4}, g0, g1); % LH and HL are reversed for dwt2
%     figure
%     showImage(RSC(rec))
%     imshow(It_3{1}-rec);
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<---------- RECOMPOSE 3-----------
% rec = idwt2(rec(1:size(It_3{1},1),1:size(It_3{1},2)), It_3{3}, It_3{2}, It_3{4}, g0, g1);
rec = idwt2(rec(1:size(It_3{1},1),1:size(It_3{1},2)), It_3{3}, It_3{2}, It_3{4}, g0, g1); 
% LH and HL are reversed for dwt2
 
%     figure
%     showImage(RSC(rec))
%     imshow(It_2{1}-rec);
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<---------- RECOMPOSE 2-----------
% rec = idwt2(rec(1:size(It_2{1},1),1:size(It_2{1},2)), It_2{3}, It_2{2}, It_2{4}, g0, g1);
rec = idwt2(rec(1:size(It_2{1},1),1:size(It_2{1},2)), It_2{3}, It_2{2}, It_2{4}, g0, g1); 
% LH and HL are reversed for dwt2
 
%     figure
%     showImage(RSC(rec))
%     imshow(It_1{1}-rec);
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<---------- RECOMPOSE 1-----------
% rec = idwt2(rec(1:size(It_1{1},1),1:size(It_1{1},2)), It_1{3}, It_1{2}, It_1{4}, g0, g1);
rec = idwt2(rec(1:size(It_1{1},1),1:size(It_1{1},2)), It_1{3}, It_1{2}, It_1{4}, g0, g1); 
% LH and HL are reversed for dwt2



%     figure
%     showImage(RSC(rec))
%     showImage(double(imagefirst)-rec);

[recNRG, recSNR]=imNRGSNR(rec);
fprintf("\n Reconstructed image Peak SNR: "+recSNR);
fprintf("\n Relative PSNR (Original vs Reconstructed) "+psnr(imagefirst,uint8(rec))+"\n")
% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< OG & Reconstruct Comparison------------
figure
subplot(1,2,1)
imshow(imagefirst);
title("Original Image")
subplot(1,2,2)
showImage(uint8(rec))
title("Reconstructed Image")
sgtitle("Comparison of Original and Reconstructed Images")
% showImage(double(imagefirst)-rec);

%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<-----FUNCTIONS-------------------


% OBSOLETE % RECOMPOSE
% function rec=recompose(it,fg1,fg2)
%     LL=[]; LH=[]; HL=[]; HH =[];
%     if isempty(it{1}) ==0
%         LL=it{1};
%     end
%     if isempty(it{2}) ==0
%         LH=it{2};
%     end
%     if isempty(it{3}) ==0
%         HL=it{3};
%     end
%     if isempty(it{4}) ==0
%         HH=it{4};
%     end
%     rec = idwt2(LL, HL, LH, HH, fg1, fg2); % LH and HL are reversed for dwt2
% end

% DECOMPOSE
function it=decompose(OGimg,fh1,fh2,it)
    [LL, HL, LH, HH] = dwt2(OGimg, fh1, fh2); % LH and HL are reversed for dwt2
    it{1}=LL;
    it{2}=LH;
    it{3}=HL;
    it{4}=HH;
    
end


%Rescale to 0, 255 range
function rescaled=RSC(im)
     global rescaleWhenPlotting
        if rescaleWhenPlotting==2
            rescaled = wcodemat(im,255);
        elseif rescaleWhenPlotting==1
            rescaled = rescale(im,0,255);
        else
            rescaled=im;
        end
end

%Converts to int to display
function []=showImage(img)
imshow(uint8(img));
end

% Show Decomposition
function []=showDecomp(IT)
    figure
    for k=1:numel(IT)
        subplot(2,2,k)
        showImage(RSC(IT{k}))
    end
end

% Show Energies
function []=showNRG(it)
    global En_list;
    fprintf("Image 1 is Approximation for Each Iteration \n")
    for e=1:4
        [en, snr]=imNRGSNR(it{e});
        disp("Energy of image "+e+" "+en+" \\ Peak SNR "+snr);
        En_list(1,(end+1))=en;
        En_list(2,(end))=snr;
    end
    fprintf('\n')
end

% Find Energy & SNR of Image
function [energy, SNR]=imNRGSNR(img)
    [W H]=size(img,1,2);
    average=sum(img(:))/(W*H);
    im_energy=0;
    for i=1:W
        for k=1:H
               im_energy=im_energy+(img(i,k)-average).^2;
        end
    end
    energy=im_energy/(W*H);     % Mean squared error
    SNR=20*log10(255)-10*log(energy);   % Peak signal to noise ratio
end

% Find max energy bands
function idx=maximums(lst,cnt)
    maxList=[];
    maxs=lst(1,:);
    maxs=maxs([2,3,4,6,7,8,10,11,12,14,15,16]);
    maxs=maxk(maxs,cnt);
    for i=1:numel(maxs)
       maxList(end+1)=maxs(i);
    end
    idx=zeros(1,cnt);
    for k=1:numel(maxList)
        idx(k)=find(lst(1,:)==maxList(k));
    end
end

% Set filter coefficients
function [h0, h1, g0, g1]=setTaps(tap)
%----------------------------------------------- 6 TAP
    h0_6tap = [0.0352 -0.0854 -0.1350 0.4599 0.8069 0.3327];
    h1_6tap = [-0.3327 0.8069 -0.4599 -0.1350 0.0854 0.0352];
    g0_6tap = [0.3327    0.8069    0.4599   -0.1350   -0.0854    0.0352];
    g1_6tap = [0.0352    0.0854   -0.1350   -0.4599    0.8069   -0.3327];
%----------------------------------------------- 8 TAP
    h0_8tap = [-0.0106    0.0329    0.0308   -0.1870   -0.0280    0.6309    0.7148    0.2304];
    h1_8tap = [-0.2304    0.7148    -0.6309   -0.0280   0.1870    0.0308    -0.0329   -0.0106];
    g0_8tap = [0.2304    0.7148    0.6309   -0.0280   -0.1870    0.0308    0.0329   -0.0106];
    g1_8tap = [-0.0106   -0.0329    0.0308    0.1870   -0.0280   -0.6309    0.7148   -0.2304];
% ----------------------------------------------- 10 TAP
    h0_10tap = [0.0033   -0.0126   -0.0062    0.0776   -0.0322   -0.2423    0.1384    0.7243    0.6038    0.1601];
    h1_10tap = [-0.1601    0.6038    -0.7243    0.1384   0.2423   -0.0322    -0.0776   -0.0062   0.0126    0.0033];
    g0_10tap = [0.1601    0.6038    0.7243    0.1384   -0.2423   -0.0322    0.0776   -0.0062   -0.0126    0.0033];
    g1_10tap = [0.0033    0.0126   -0.0062   -0.0776   -0.0322    0.2423    0.1384   -0.7243    0.6038   -0.1601];
    
    
    switch tap
        case 6
            h0=h0_6tap;h1=h1_6tap;g0=g0_6tap;g1=g1_6tap;
        case 8
            h0=h0_8tap;h1=h1_8tap;g0=g0_8tap;g1=g1_8tap;
        case 10
            h0=h0_10tap;h1=h1_10tap;g0=g0_10tap;g1=g1_10tap;
        otherwise
            disp("INVALID Number OF TAPS")
            return
    end
    tap
end