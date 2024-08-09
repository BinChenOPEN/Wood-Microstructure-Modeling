clear all, close all, clc
folder_save = 'SaveBirch';
load("volume_up_sampling.mat");

volume_up_sampling = volume_up_sampling(:,:,1:100);
volume_up_sampling_binary = zeros(size(volume_up_sampling),"logical");
% the threshold of the image intensity to get the binary image. If you use
% a large threshould, you will have a smaller cell wall thickness.
volume_up_sampling_binary(find(volume_up_sampling>0.5)) = 1;

% if the cell wall thickness is large. We first remove thickness of initial_thickness_removing
initial_thickness_removing = 7;
se = strel('cube',initial_thickness_removing);
Middle_Lamela = imerode(volume_up_sampling_binary,se);

porosity = 1;

while porosity>0.085
    % for i = 1:10
    % Change the parameters here. Larger value lead to thicker minimum middle lamela.
    conv_kernel = ones(21,21,21);
    conv_kernel = conv_kernel/prod(size(conv_kernel));
    % the convoluaiton take a lot of time for large structure
    C = convn(double(Middle_Lamela),conv_kernel,'same');

    % threshold should be smaller than 0.5. O.4 maybe is good.
    threshold = 0.4;
    D = Middle_Lamela&(C<=threshold);
    se = strel('cube',3);
    J1 = imerode(Middle_Lamela,se);

    %
    Middle_Lamela = J1|D;
    porosity = sum(double(Middle_Lamela(:)))/prod(size(Middle_Lamela))
    save(fullfile(folder_save,"Middle_Lamela.mat"),"Middle_Lamela","-v7.3");
end


