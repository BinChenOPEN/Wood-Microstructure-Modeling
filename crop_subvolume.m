folder = 'SaveBirch/LocalDistVolume';
folder_upsampling = 'SaveBirch/UpSampling';
mkdir(folder_upsampling);

file_all = dir(fullfile(folder,"*.bmp"));

if length(file_all) == 0
    file_all = dir(fullfile(folder,"*.tif"));
end

%% Crop a subvolume
volume_size = [1,1,1]*200;
region_origin = [201,201,101];
volume_crop = zeros(volume_size);
scale = 5;
t = 0;
for i = region_origin(3):volume_size(3)+region_origin(3)-1
    t = t+1;
    img = imread(fullfile(file_all(i).folder,file_all(i).name));
    volume_crop(:,:,t) = double(img(region_origin(1):volume_size(1)+region_origin(1)-1,region_origin(2):volume_size(2)+region_origin(2)-1))/255;

end

x_vector_interp = 1:1/scale:volume_size(1);
y_vector_interp = 1:1/scale:volume_size(2);
z_vector_interp = 1:1/scale:volume_size(3);
if 1
    [x,y,z] = meshgrid(1:volume_size(1),1:volume_size(2),1:volume_size(3));
    [xq,yq,zq] = meshgrid(x_vector_interp,y_vector_interp,z_vector_interp);
    volume_up_sampling = interp3(x,y,z,volume_crop,xq,yq,zq);

else

    %% Upsampling along in x y plane using interpolation
    volume_up_sampling = zeros(length(x_vector_interp),length(x_vector_interp),volume_size(3));
    parfor t = 1:volume_size(3)
        [X,Y] = meshgrid(1:volume_size(1),1:volume_size(2));
        [Xq,Yq] = meshgrid(x_vector_interp,y_vector_interp);
        Vq = interp2(X,Y,squeeze(volume_crop(:,:,t)),Xq,Yq);
        volume_up_sampling(:,:,t) = Vq;
    end


    %% Upsampling along z direction using interpolation

    parfor t = 1:volume_size(2)
        [X,Z] = meshgrid(x_vector_interp,1:volume_size(3));
        [Xq,Zq] = meshgrid(x_vector_interp,z_vector_interp);
        Vq = interp2(X,Z,squeeze(volume_up_sampling(:,T,:)),Xq,Zq);

        FileName = fullfile(folder_upsampling,['volImgRef_Y_',sprintf('%05d',t),'.bmp']);
        imwrite(Vq,FileName)
    end
end
save("volume_crop.mat","volume_crop","-v7.3");
save("volume_up_sampling.mat","volume_up_sampling","-v7.3");
