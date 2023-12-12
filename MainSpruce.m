clear, close all, clc,warning('off','all');
%% set parameters
% Volume structure size;
sizeVolume    = [1500,1500,750];

% The volume structrue is slightly enlarged for image interpolation. D
extraSZ   = [200,400,150];
sizeImEnlarge = sizeVolume+extraSZ; % The simulated structure is a bit larger than the expected one.
% The extra region will be removed afterwards.
SaveFolder = 'SaveSpruce/'; % Folder to save the data
mkdir(SaveFolder);


cellR              = 37.75; % The grid distance for the nodes we generated. In Unit of voxels.
% It is roughly the radius of the fibers.
xVector            = 5:cellR:sizeImEnlarge(1)-5; % the grid node points
yVector            = 5:cellR:sizeImEnlarge(2)-5; % the grid node points
[xGrid, yGrid]     = ndgrid(xVector,yVector);
cellLength         = 4877; % Average fiber length
cellLengthVariance = 1219; % Standard deviation of fiber length
vesselLength         = 4877; % Standard deviation of vessel length
vesselLengthVariance = 1219; % Standard deviation of vessel length
sliceInterestSpace   = 100; % We generate one slice every 100 slice.
% Along L direction, the intensity is interpolated from these slices
PeriodParameter      = 1000; % This parameter is related to the period of the year ring size
% ray cell basic length
rayCellLength      = 149.4; % ray cell length along radial direction
% ray cell disturbance
rayCell_variance   = 38.5; % ray cell length deviation along radial direction
% cell end thickness;
cellEndThick       = 4; % End of cell wall thickness along L direction
cellWallThick      = 4; % Cell wall thickness
% vessel wall thickness is thicker than general cell wall
vesselThicker      = 0; % Assume vessel is thicker than ray cells. It could be any value
rayHeight          = 50.2; % the width of the ray cell
% four neighbor grids. Do not change here
neighborLocal      = [-1, 0, 1, 0;
    0, -1, 0, 1];
% the grid node number
rayCellNum         = 8.44; % ray cell count in a group
rayCellNumStd      = 4.39; % ray cell count in a group
saveVolumeAs3D     = 0; % 1 or 0. If the volume size is too large. This operation may out of memory.
numGridNodes       = length(xGrid(:));

% For data transfer
Params.sizeIm         = sizeVolume;
Params.cellR          = cellR;
Params.xGrid          = xGrid;
Params.yGrid          = yGrid;
% Params.cellLengthDist = cellLengthDist;
Params.cellLength     = cellLength;
Params.rayCellLength  = rayCellLength;
Params.rayCell_variance  = rayCell_variance;
Params.cellEndThick   = cellEndThick;
Params.neighborLocal  = neighborLocal;
Params.rayHeight       = rayHeight;
Params.sizeImEnlarge  = sizeImEnlarge;
Params.cellThick      = cellWallThick;
Params.gridSize       = [length(xVector),length(yVector)];



%% Generate the distortion map for early wood and late wood
thickerAll    = [];
compress_all  = [];
Temp = 0;
for i = 1:20
    % Cell wall thickness distribution is different in early wood and late
    % wood
    featureSize1 = round(PeriodParameter+PeriodParameter*rand(1)); % The size of a year ring is roughly controled by this
    thickAmptitude = (12+6*rand(1)); % Control the amptitude of the cell wall thickness.
    thicker1 = thickAmptitude/featureSize1^4*(round(0.4*featureSize1):featureSize1)'.^4; % They are added to the original one
    thicker1 = thicker1-thicker1(1);
    thickerAll = [thickerAll;thicker1];
    
    % The fiber size is different in early wood and late wood. Their size
    % can be controled by applying a compression to the whole field. The
    % compression is designed
    k         = 0.6/3*featureSize1; % Change this value to control the extent of compression
    compress  = -k/featureSize1^3*(round(0.4*featureSize1):featureSize1)'.^3; % Control the extent of fiber compression in different regions
    compress  = compress + Temp - compress(1);
    compress_all = [compress_all;compress];
    Temp      = compress_all(end);
end

indStart = 2*PeriodParameter; % This value is not important. It should be roungly this value
thickerAll_validSub   = thickerAll(indStart+1:indStart+sizeImEnlarge(1));
compress_all_validSub = compress_all(indStart+1:indStart+sizeImEnlarge(1));
compress_all_validSub = compress_all_validSub - compress_all_validSub(1)...
    -(compress_all_validSub(end) - compress_all_validSub(1))*(1:sizeImEnlarge(1))'/sizeImEnlarge(1);

figure(1),
plot(thickerAll_validSub(extraSZ(1)/2+1:sizeVolume(1)+extraSZ(1)/2),'linewidth',1);
hold on,
plot(compress_all_validSub(extraSZ(1)/2+1:sizeVolume(1)+extraSZ(1)/2),'linewidth',1);
axis tight,
set(gcf,'color','w');
legend('Cell wall thickness increment','Compression','location','northeast');
legend('boxoff');
xlabel('R coordinate (voxels)');
ylabel('Y (voxels)');
ylim([-40,30])
set(gca,'fontsize',12),
saveas(gcf,[SaveFolder,'ThickDistribution'],'epsc');




% Specify the location of grid nodes and the thickness (with disturbance)
% for some typical slice
t  = 0;
sliceInterest = [1:sliceInterestSpace:sizeImEnlarge(3),sizeImEnlarge(3)];
for iSlice = sliceInterest
    t         = t+1;
    % Random location for the nodes. 3 and 1.5 are used to control the randomness of the nodes
    xGrid_interp(:,t)     = xGrid(:)+rand(numGridNodes,1)*3-1.5;
    yGrid_interp(:,t)     = yGrid(:)+rand(numGridNodes,1)*3-1.5;
    
    thicker_inerp           = thickerAll_validSub(round(xGrid_interp(:,t)));
    Thickness_interp_ray(:,t) = (cellWallThick-0.5)*ones(numGridNodes,1)+rand(numGridNodes,1)*1;
    
    Thickness_interp_fiber(:,t) = Thickness_interp_ray(:,t)+thicker_inerp;
end

% Interp the grid nodes location and thickness for all slices
parfor i = 1:length(xGrid(:))
    xGrid_all(i,:)     = spline(sliceInterest,xGrid_interp(i,:),1:sizeImEnlarge(3));
    yGrid_all(i,:)     = spline(sliceInterest,yGrid_interp(i,:),1:sizeImEnlarge(3));
    Thickness_all_fiber(i,:) = spline(sliceInterest,Thickness_interp_fiber(i,:),1:sizeImEnlarge(3));
    Thickness_all_ray(i,:) = spline(sliceInterest,Thickness_interp_ray(i,:),1:sizeImEnlarge(3));
end
Params.xGrid_all     = xGrid_all;
Params.yGrid_all     = yGrid_all;
Params.Thickness_all = Thickness_all_fiber;
Params.Thickness_all_ray = Thickness_all_ray;

if 0
    h0 = figure;
    plot(xGrid_all(1,:),1:length(xGrid_all(1,:)),'k.');
    hold on,
    plot(xGrid_all(2,:),1:length(xGrid_all(2,:)),'k.');
    set(gcf,'color','w'),axis tight,xlabel('X (voxel)'),ylabel('Z (voxel)')
    saveas(gcf,[SaveFolder,'grid_distri_Z'],'epsc');
end

if 0
    h = figure;
    plot(xGrid_all(1,:),1:length(Thickness_all_fiber(1,:)),'k-');
    set(gcf,'color','w'),axis tight,xlabel('X (voxel)'),ylabel('Z (voxel)')
    saveas(gcf,[SaveFolder,'thick_distri_Z'],'epsc');
end

if 0
    h1 = figure;
    plot(xGrid_interp(:,[1]),yGrid_interp(:,[1]),'r.');
    set(gcf,'color','w'),axis equal,axis tight,xlabel('X (voxel)'),ylabel('Y (voxel)')
    saveas(gcf,[SaveFolder,'grid'],'epsc');
end

% location of ray cells along T direction.
raycell_linspace = 9:10:length(xVector)-8;
raycellXindAll = raycell_linspace+rand(1,length(raycell_linspace))*10-5;
raycellXindAll = floor(raycellXindAll./2)*2+1;

% Specify the location of the vessels. The vessels should only be given at
% some positions. Addtionally, it shouold be within the grid.
x_rand_1 = round(((rand(60,1)*(length(xVector)-16))+8)/2)*2-1;
y_rand_1 = round(((rand(60,1)*(length(yVector)-14))+7)/4)*4;
y_rand_2 = round(((rand(50,1)*(length(yVector)-14))+7)/2)*2;
x_rand_2 = round(((rand(50,1)*(length(xVector)-16))+8)/2)*2;

x_rand_all = [x_rand_1;x_rand_2];
y_rand_all = [y_rand_1;y_rand_2];
Vessel_all = [x_rand_all,y_rand_all];

% Create a large volume image with all voxels being 0 (black)
volImgRef = 255*ones(sizeImEnlarge,'uint8');



t = 0;
s = 0;
raycellWidth = [];
% raycellXindAll means the ray cell group index along y dirction
for i = 1:length(raycellXindAll)
    
    % generate the raycell groups along z direciton
    m = sizeImEnlarge(3)/rayCellNum/rayHeight+5;
    raycellWidth_temple{i} = 0;
    
    % the space from the first stack of the rays to the top slice
    %% change made here
    rayCellSpace   = round(16*rand(ceil(m),1)+6);
    for j = 1 : m
        s = s+1;
        % have two column of ray cells
        raycellXind     = [raycellXindAll(i),raycellXindAll(i)+1];
        
        % for each column, it should have some ray cell stacks
        rayCellNumGroup      = round(rayCellNum+randn(1)*rayCellNumStd);
        if rayCellNumGroup <5
            rayCellNumGroup = 5;
        end
        if rayCellNumGroup >25
            rayCellNumGroup = 25;
        end
        
        %% for every positon along x direction, it may have several groups
        
        if j == 1
            raycellWidth_temple{i}  = raycellWidth_temple{i}(end)+...
                ([0:rayCellNumGroup]+rayCellSpace(j)+round(-30*rand(1)))*rayHeight;
        else
            raycellWidth_temple{i}  = raycellWidth_temple{i}(end)+...
                ([0:rayCellNumGroup]+rayCellSpace(j))*rayHeight;
        end
        
        %% Generate the ray cells in this group
        %% change made here
        if raycellWidth_temple{i}(1)<=sizeImEnlarge(3)-150 && raycellWidth_temple{i}(end)>=150
            t = t+1;
            raycellXindAll_update(t) = raycellXindAll(i);
            raycellWidth{t}  = round(raycellWidth_temple{i});
            keepRayCell(t) = i;
            raycellXind_all{t} = raycellXind;
        end
    end
end
% Generate small cells.
t = 0;
skipFiberColumn = [raycellXindAll,raycellXindAll+1];

for i = 2:2:length(xVector)-2
    for j = 2:2:length(yVector)-2
        if isempty(find(skipFiberColumn == j))
            t = t+1;
            
            % To control the length of the fibers
            fiberEndLocAll    = [];
            Initial = round(rand(1)*cellLength);
            fiberEndLocAll(1) = Initial;
            for s = 1:ceil(max((sizeImEnlarge(3)+4*cellLength)/cellLength+3))
                temp = min(3*cellLength,max(100,cellLength+randn(1)*cellLengthVariance));
                fiberEndLocAll(s+1) = round(fiberEndLocAll(s)+temp);
            end
            % This is a manully given value. To increase the randomness
            fiberEndLoc = fiberEndLocAll-4*cellLength+1;
            % The end of the vessel should be inside the volume.
            fiberEndLoc = fiberEndLoc(fiberEndLoc>=4 & fiberEndLoc<=sizeImEnlarge(3)-4);
            
            
            
            fiberEnd  = [];
            if ~isempty(fiberEndLoc)
                for iThick = 0:cellEndThick-1
                    fiberEnd = [fiberEnd;fiberEndLoc+iThick];
                end
            end
            
            i1 = i;
            % The arrangement of the cells should be stager. So every four
            % nodes, they should deviate along x direction.
            if mod(j,4) == 0
                i1 = i+1;
            end
            % Convert subscript to index of this cell
            ind  = length(xVector)*(j-1)+i1;
            % Skip some fibres
            skipCellThick = 0;
            for iSlice = 1:sizeImEnlarge(3)
                % used to store the four neighbored points for elipse fitting
                pointCoord = [];
                % small fibers
                % find the four neighbored nodes for elipse fitting
                for k = 1:4
                    indxNeighbor    = length(xVector)*(j-1+neighborLocal(2,k))+i1+neighborLocal(1,k);
                    pointCoord(k,:) = [xGrid_all(indxNeighbor,iSlice),yGrid_all(indxNeighbor,iSlice)];
                end
                % make the elipse more elipse
                if skipCellThick == 0
                    pointCoord(2,2)     = pointCoord(2,2)-1.5;
                    pointCoord(4,2)     = pointCoord(4,2)+1.5;
                end
                
                C_elipse                = fitElipse(pointCoord); % Estimate the coeffecients of the elipse
                % Then we can estimate the diameter along two direction for the fiber
                % The rectangle region covering the elipse.
                regionCellindx = ceil(ceil(C_elipse(3))+[-max(floor(C_elipse(1:2))):max(floor(C_elipse(1:2)))]);
                regionCellindy = ceil(ceil(C_elipse(4))+[-max(floor(C_elipse(1:2))):max(floor(C_elipse(1:2)))]);
                regionCellindx = max(1,min(regionCellindx)):min(max(regionCellindx),sizeImEnlarge(1));
                regionCellindy = max(1,min(regionCellindy)):min(max(regionCellindy),sizeImEnlarge(2));
                [regionCellx,regionCelly] = ndgrid(regionCellindx,regionCellindy);
                
                % External contour
                inElipse1      = (regionCellx-C_elipse(3)).^4./C_elipse(1)^4 + (regionCelly-C_elipse(4)).^4./C_elipse(2)^4;
                % Internal contour
                inElipse2      = (regionCellx-C_elipse(3)).^4./(C_elipse(1)-Thickness_all_fiber(ind,iSlice)-skipCellThick).^4 + ...
                    (regionCelly-C_elipse(4)).^4./(C_elipse(2)-Thickness_all_fiber(ind,iSlice)-skipCellThick).^4;
                %                 regionCell     = zeros(length(regionCellindx),length(regionCellindy));
                regionCell     = ones(length(regionCellindx),length(regionCellindy),'uint8');
                if isempty(find(iSlice==fiberEnd))
                    % if this slice is not the end of this cell
                    % inside the lumen, it should be black
                    regionCell(find((inElipse2<=1))) = 0;
                else
                    % if this slice is the end of this cell,
                    % all points should be white
                end
                volImgRef(regionCellindx,regionCellindy,iSlice) = volImgRef(regionCellindx,regionCellindy,iSlice).*regionCell;
            end
            
        end
    end
end
volImgRef(volImgRef>255) = 255;
if 0
    h2 = figure;
    imshow(volImgRef(:,:,1))
    set(gcf,'color','w'),axis equal,axis tight,xlabel('X (voxel)'),ylabel('Y (voxel)')
    saveas(h2,[SaveFolder,'surf1'],'epsc');
end


%% Generate ray cells
Folder_Init = [SaveFolder,'/volImgInit/'];
mkdir(Folder_Init);
parfor z = 1:size(volImgRef,3)
    imgZ = volImgRef(:,:,z);
    zStr = num2str(z);
    numStr = '00000';
    numStr(end-length(zStr)+1:end) = zStr;
    FileNameAll{z} = ['volImgRef_',numStr,'.bmp'];
    imgName = fullfile(Folder_Init,FileNameAll{z});
    imwrite(imgZ,imgName);
end
volImgRef_final = volImgRef;
clear volImgRef;

for t = 1:length(raycellWidth)
    volImgRef_final  = rayCellGenerateSpruce(volImgRef_final,raycellXind_all{t},raycellWidth{t},Params);
end

if 0
    h4 = figure;
    imshow(squeeze(volImgRef_final(240,:,:)))
    set(gcf,'color','w'),axis equal,axis tight,
    saveas(h4,[SaveFolder,'slice1_raycell'],'epsc');
end


%% Add complicated deformation to the volume image
[x,y]  = ndgrid(1:sizeImEnlarge(1),1:sizeImEnlarge(2));
u      = zeros(sizeImEnlarge(1:2));
v      = zeros(sizeImEnlarge(1:2));
u1     = zeros(sizeImEnlarge(1:2));
v1     = zeros(sizeImEnlarge(1:2));
t      = 0;
for i = 2:2:length(xVector)-1
    for j = 2:2:length(yVector)-1
        t = t+1;
        i1 = i;
        pointCoord = [];
        
        if mod(j,4) == 0
            i1 = i+1;
        end
        ind  = length(xVector)*(j-1)+i1;
        
        % Small cells
        islarge    = 0;
        indxPt     = length(xVector)*(j-1)+i1;
        pointCoord = [xGrid_all(indxPt,1),yGrid_all(indxPt,1)];
        
        %----------the value here are super important-------------%
        k          = [0.1,0.08,1.5+rand(1),3+3*rand(1)];
        
        u_temp     = localDistort(x,y,pointCoord(1),pointCoord(2),k);
        
        %----------the value here are super important-------------%
        k          = [0.1,0.08,1.5+rand(1),2*(1+rand(1))];
        
        v_temp     = localDistort(y,x,pointCoord(1),pointCoord(2),k);
        
        u          = u + sign(randn(1))*u_temp;
        v          = v + v_temp;
        
        if rand(1)<1/100
            % Deform for a large regions
            islarge    = 1;
            indxPt     = length(xVector)*(j-1)+i1;
            pointCoord = [xGrid_all(indxPt,1),yGrid_all(indxPt,1)];
            
            %----------the value here are super important-------------%
            k          = [0.015,0.01,3*(1+0.5*rand(1)),1.5*(1+1*rand(1))];
            u_temp     = localDistort(x,y,pointCoord(1),pointCoord(2),k);
            u1          = u1 + sign(randn(1))*u_temp;
            
            %----------the value here are super important-------------%
            k          = [0.02,0.015,4+2*rand(1),2*(1+1*rand(1))];
            v_temp     = localDistort(y,x,pointCoord(1),pointCoord(2),k);
            v1          = v1 + sign(randn(1))*v_temp;
        end
    end
end


if 1
    h5 = figure;
    surf(u(end:-1:1,:)),colorbar,view([0,90])
    set(gcf,'color','w'),shading interp, axis equal,axis tight,xlabel('Y (voxel)'),ylabel('X (voxel)')
    saveas(h5,[SaveFolder,'def_local'],'epsc');
end
% image interpolation

Folder_ImgSeries = [SaveFolder,'/volImgBackBone/'];
mkdir(Folder_ImgSeries);
parfor z = 1:size(volImgRef_final,3)
    imgZ = volImgRef_final(:,:,z);
    imgName = fullfile(Folder_ImgSeries,FileNameAll{z});
    imwrite(imgZ,imgName);
end
clear volImgRef_final;
clear volImgRef_final;

Folder_ImgSeries_LocalDist = [SaveFolder,'/LocalDistVolume/'];
mkdir(Folder_ImgSeries_LocalDist)

%% Impose local distortion by image interpolation
parfor z = 1:sizeImEnlarge(3)
    
    imgName = fullfile(Folder_ImgSeries,FileNameAll{z});
    ImCover = double(imread(imgName));
    
    [x_grid,y_grid]      = ndgrid(1:sizeImEnlarge(1),1:sizeImEnlarge(2));
    v_all_raycell = rayCell_shrinking_spruce(...
        Params,yVector,sizeImEnlarge,raycellWidth,raycellXindAll,rayHeight,cellR-1,z);
    
    u_compress_Mat = repmat(compress_all_validSub,1,sizeImEnlarge(2));
    x_interp = x_grid+u+u_compress_Mat;
    y_interp = y_grid+v+v_all_raycell;
    
    Vq       = griddata(x_interp(:),y_interp(:),ImCover(:),x_grid(:),y_grid(:),'linear');
    imInterp = reshape(Vq,sizeImEnlarge(1:2));
    
    imInterp(isnan(imInterp)) = 0;
    imInterp(imInterp<0)      = 0;
    imInterp(imInterp>0)    = 255;
    ImgTemp           = uint8(imInterp);
    
    imgName = fullfile(Folder_ImgSeries_LocalDist,FileNameAll{z});
    imwrite(ImgTemp,imgName);
end
% toc

%% image interpolation
Folder_ImgSeries_GlobalDist = [SaveFolder,'/GlobalDistVolume/'];
mkdir(Folder_ImgSeries_GlobalDist)

for t = 1:length(sliceInterest)-1
    if t == length(sliceInterest)-1
        indxZ    = sliceInterest(t):sliceInterest(t+1);
    else
        indxZ    = sliceInterest(t):sliceInterest(t+1)-1;
    end
    volImgLocalDistSub = [];
    for z = indxZ
        imgName = fullfile(Folder_ImgSeries_LocalDist,FileNameAll{z});
        imgZ = single(imread(imgName));
        volImgLocalDistSub = cat(3,volImgLocalDistSub,imgZ);
    end
    
    [x_grid,y_grid,z_grid] = ndgrid(1:single(sizeImEnlarge(1)),1:single(sizeImEnlarge(2)),single(indxZ));
    v_allz   = (x_grid-sizeImEnlarge(1)*2/5).^3/1e7/2+(y_grid-sizeImEnlarge(2)/2).*(z_grid-sizeImEnlarge(3)/2)/(10^4)/5;
    u_allz   = (y_grid-sizeImEnlarge(1)/3).^2/1e4/3+(x_grid-sizeImEnlarge(2)/3).*(z_grid-sizeImEnlarge(3)/3)/(10^4)/4;
    % w_allz   = (z_grid-sizeImEnlarge(3)/2).*(x_grid-sizeImEnlarge(1)/2)/1e4+5*sin(pi*x_grid/(sizeImEnlarge(1)/2)+z_grid/(sizeImEnlarge(3)/2));
    
    x_interp = x_grid-u_allz-repmat(u1,1,1,size(u_allz,3));
    y_interp = y_grid-v_allz;
    z_interp = z_grid;
    
    if t == 1
        if 1
            h6 = figure;
            surf(v_allz(end:-1:1,:,1)),hc = colorbar; view([0,90]);
            set(gcf,'color','w'),shading interp, axis equal,axis tight,xlabel('Y (voxel)'),ylabel('X (voxel)')
            saveas(h6,[SaveFolder,'def_global'],'epsc');
        end
    end
    
    % interpolation
    Vq       = uint8(interpn(x_grid,y_grid,z_grid,volImgLocalDistSub,...
        x_interp(:),y_interp(:),z_interp(:),'linear'));
    VolImg_Temp   = reshape(Vq,[sizeImEnlarge(1),sizeImEnlarge(2),length(indxZ)]);
    k = 0;
    for z = indxZ
        k = k+1;
        imgZ = VolImg_Temp(:,:,k);
        imgName = fullfile(Folder_ImgSeries_GlobalDist,FileNameAll{z});
        imwrite(imgZ,imgName);
    end
end
figure,
imshow(VolImg_Temp(:,:,1))
set(gcf,'color','w');

%% Save the center part of the volume
Folder_FinalVolumeSlice = [SaveFolder,'/FinalVolumeSlice/'];
mkdir(Folder_FinalVolumeSlice)

for i = round(extraSZ(3)/2)+1:sizeVolume(3)+round(extraSZ(3)/2)
    imgName = fullfile(Folder_ImgSeries_GlobalDist,FileNameAll{i});
    img = imread(imgName);
    imgName_save = fullfile(Folder_FinalVolumeSlice,FileNameAll{i-round(extraSZ(3)/2)});
    imwrite(img(round(extraSZ(1)/2)+1:sizeVolume(1)+round(extraSZ(1)/2),...
        round(extraSZ(2)/2)+1:sizeVolume(2)+round(extraSZ(2)/2)),imgName_save);
end

%% Save the center part of the volume
if saveVolumeAs3D
    Folder_FinalVolume = [SaveFolder,'/FinalVolume3D/'];
    mkdir(Folder_FinalVolume)
    for i = 1:sizeVolume(3)
        imgName = fullfile(Folder_FinalVolumeSlice,FileNameAll{i});
        img = imread(imgName);
        volumeFinal(:,:,i) = img;
    end
    save(fullfile(Folder_FinalVolume,'FinalVolume.mat'),'volumeFinal');
end
%%
