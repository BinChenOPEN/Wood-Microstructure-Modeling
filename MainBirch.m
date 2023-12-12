clear, close all, clc,warning('off','all');
%% set parameters
% Volume structure size;
sizeVolume    = [1500,1500,750];

% The volume structrue is slightly enlarged for image interpolation
extraSZ   = [200,400,150];
sizeImEnlarge = sizeVolume+extraSZ; % The simulated structure is a bit larger than the expected one.
% The extra region will be removed afterwards.
SaveFolder = 'SaveBirch/'; % Folder to save the data
mkdir(SaveFolder);


% cellwall2radiusRatio  = 1.38/(6.01+1.38);
cellR              = 14.5; % The grid distance for the nodes we generated. In Unit of voxels.
                         % It is roughly the radius of the fibers.
xVector            = 5:cellR:sizeImEnlarge(1)-5; % the grid node points
yVector            = 5:cellR:sizeImEnlarge(2)-5; % the grid node points
[xGrid, yGrid]     = ndgrid(xVector,yVector);

cellLength         = 2341;  % Average fiber length
cellLengthVariance = 581;  % Standard deviation of fiber length
sliceInterestSpace = 100;  % We generate one slice every 100 slice.
                           % Along L direction, the intensity is interpolated from these slices
vesselLength       = 780; % Average vessel length
vesselLengthVariance = 195; % Standard deviation of vessel length

% ray cell basic length
rayCellLength      = 62; % ray cell length along radial direction
% ray cell disturbance
rayCell_variance   = 15; % ray cell length deviation along radial direction
cellEndThick       = 4;  % End of cell wall thickness along L direction
cellWallThick      = 4;  % Cell wall thickness
isExistVessel      = 1;  % Include vessels, 1: with vessel, 0: no vessel
isExistRayCell     = 1;  % Is there any ray cells, 1: exist ray cells; 0: no ray cell
% vessel wall thickness is thicker than general cell wall, 0 means the
% sames
vesselThicker      = 1;  % Assume vessel is thicker than ray cells
raywHeight         = 42; % the height of the ray cell along fiber direction
raySpace           = 20; % the space between ray cell along T direction. The distance is raySpace*cellR
% four neighbor grids. Do not change here
neighborLocal      = [-1,  0, 1, 0;
                       0, -1, 0, 1];
% the grid node number
rayCellNum         = 11.33; % ray cell count in a group of ray cells
rayCellNumStd      = 3.39;  % ray cell count in a group of ray cells
VesselCount        = 50;   % This number is used to control the vessel number and distribution.
% Please tune with this parameter
writeLocalDeformData  = 0; % 1 or 0.
writeGlobalDeformData = 0; % 1 or 0.
saveVolumeAs3D        = 1; % 1 or 0. If the volume size is too large. This operation may out of memory.
numGridNodes          = length(xGrid(:)); % grid node number

% Store all parameters into structure Params
Params.sizeIm         = sizeVolume;
Params.cellR          = cellR;
Params.xGrid          = xGrid;
Params.yGrid          = yGrid;
% Params.cellLengthDist = cellLengthDist;
Params.cellLength     = cellLength;
Params.rayCellLength  = rayCellLength;
Params.rayCell_variance = rayCell_variance;
Params.cellEndThick   = cellEndThick;
Params.neighborLocal  = neighborLocal;
Params.rayHeight      = raywHeight;
Params.sizeImEnlarge  = sizeImEnlarge;
Params.cellThick      = cellWallThick;
Params.gridSize       = [length(xVector),length(yVector)];

%%
% Specify the location of grid nodes and the thickness (with disturbance)
% for some typical slices. Generate for every 100 slice
t  = 0;
sliceInterest = [1:sliceInterestSpace:sizeImEnlarge(3),sizeImEnlarge(3)];
for iSlice = sliceInterest
    t         = t+1;
    % Random location for the nodes. 3 and 1.5 are used to control the randomness of the nodes
    xGrid_interp(:,t)     = xGrid(:)+rand(numGridNodes,1)*3-1.5;
    yGrid_interp(:,t)     = yGrid(:)+rand(numGridNodes,1)*3-1.5;
    % A little bit randomness is included in the cell wall thickness.
    Thickness_interp(:,t) = (cellWallThick-0.5)*ones(numGridNodes,1)+rand(numGridNodes,1)*1;
end

% Interp the grid nodes location and thickness for all slices using spline
% interpolation
parfor i = 1:length(xGrid(:))
    xGrid_all(i,:)     = spline(sliceInterest,xGrid_interp(i,:),1:sizeImEnlarge(3));
    yGrid_all(i,:)     = spline(sliceInterest,yGrid_interp(i,:),1:sizeImEnlarge(3));
    Thickness_all(i,:) = spline(sliceInterest,Thickness_interp(i,:),1:sizeImEnlarge(3));
end
Params.xGrid_all     = xGrid_all;
Params.yGrid_all     = yGrid_all;
Params.Thickness_all = Thickness_all;


%% change the following to if 1 to see and save some figures.
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
    plot(xGrid_all(1,:),1:length(Thickness_all(1,:)),'k-');
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
raycell_linspace = 20:raySpace:length(yVector)-10;
if isExistRayCell
    raycellXindAll = raycell_linspace+rand(1,length(raycell_linspace))*10-5;
    raycellXindAll = floor(raycellXindAll./2)*2+1;
end

% Specify the location of the vessels. The vessels should only be given at
% some positions. Addtionally, it shouold be within the grid.
x_rand_1 = round(((rand(VesselCount,1)*(length(xVector)-16))+8)/2)*2-1;
y_rand_1 = round(((rand(VesselCount,1)*(length(yVector)-14))+7)/4)*4;
y_rand_2 = round(((rand(VesselCount/2,1)*(length(yVector)-14))+7)/2)*2;
x_rand_2 = round(((rand(VesselCount/2,1)*(length(xVector)-16))+8)/2)*2;

x_rand_all = [x_rand_1;x_rand_2];
y_rand_all = [y_rand_1;y_rand_2];
Vessel_all = [x_rand_all,y_rand_all];

%% tune with location of the vessels. The vessels should not be too close
% Remove some vessel that too close to the other vessels
for i = 1:length(Vessel_all)
    Vessel_all_temp = Vessel_all(i,:);
    %     Vessel_all_temp(i,:) = [];
    Vessel_all_dist =  Vessel_all-Vessel_all_temp;
    if Vessel_all(i,1) == 0
    else
        mark0 = find(Vessel_all_dist(:,1)<=6 &...
            Vessel_all_dist(:,1)>=-6 &...
            abs(Vessel_all_dist(:,2))<=4);
        Vessel_all(mark0,:) = zeros(length(mark0),2);
    end
    Vessel_all(i,:) = Vessel_all_temp;
end

Vessel_all(find(Vessel_all(:,1)==0),:) = [];

% The vessels should not be too close to the ray cells regions
m = 0;
Vessel =[];
for i = 1:length(Vessel_all)
    t = 0;
    if isExistRayCell
        for j = 1:length(raycellXindAll)
            if Vessel_all(i,2)-raycellXindAll(j)>=-3 && Vessel_all(i,2)-raycellXindAll(j)<=4
                t = t+1;
            end
        end
    end
    if t == 0
        m = m+1;
        Vessel(m,:) = Vessel_all(i,:);
    end
end
Vessel_all = Vessel;


% The vessels are extended. Some vessels are extended to be
% double-vessel cluster, some are triple-vessel clusters.
Vessel_all_extend = [];
for i = 1:length(Vessel_all)
    
    Vessel_all_temp = Vessel_all(i,:);
    Vessel_all_dist = Vessel_all-Vessel_all_temp;
    
    mark0 = find(Vessel_all_dist(:,1)<=24 &...
        Vessel_all_dist(:,1)>=-8 & abs(Vessel_all_dist(:,2))<=8);
    
    mark1 = find(Vessel_all_dist(:,1)<=12 &...
        Vessel_all_dist(:,1)>=-6 & abs(Vessel_all_dist(:,2))<=6);
    
    if length(mark0) > 1
        Vessel_all_extend = [Vessel_all_extend;Vessel_all(i,:)];
        possibility = rand(1);
        if length(mark1) <=1
            if possibility<0.2
                temp = [Vessel_all(i,1)+6+sign(randn(1)),Vessel_all(i,2)+sign(randn(1))*2];
                Vessel_all_extend = [Vessel_all_extend;temp];
            else
                if possibility<0.5
                    temp = [Vessel_all(i,1)+6,Vessel_all(i,2)];
                    Vessel_all_extend = [Vessel_all_extend;temp];
                end
                
            end
        end
    else
        if Vessel_all(i,1)+12<length(xVector) && Vessel_all(i,2)+10<length(yVector)
            temp0 = [Vessel_all(i,1)+5+sign(randn(1)),Vessel_all(i,2)];
            possibility = rand(1);
            if possibility<0.3
                temp  = [temp0;temp0(1)+5,temp0(2)+2*sign(randn(1))];
            else
                temp  = [temp0;temp0(1)+5+sign(randn(1)),temp0(2)];
            end
            Vessel_all_extend = [Vessel_all_extend;Vessel_all(i,:);temp];
        else
            Vessel_all_extend = [Vessel_all_extend;Vessel_all(i,:)];
        end
        
    end
end
Vessel_all = Vessel_all_extend;

% check the extended vessel locations again. They should satisfy the former
% two conditions
m = 0;
Vessel =[];
if isExistVessel
    for i = 1:length(Vessel_all)
        t = 0;
        if isExistRayCell
            for j = 1:length(raycellXindAll)
                if Vessel_all(i,2)-raycellXindAll(j)>=-3 && Vessel_all(i,2)-raycellXindAll(j)<=4
                    t = t+1;
                end
            end
        end
        if t == 0
            if Vessel_all(i,2) <= length(yVector)-3 && ...
                    Vessel_all(i,1) <= length(xVector)-3
                m = m+1;
                Vessel(m,:) = Vessel_all(i,:);
            end
        end
    end
end
%
%% The fibers in the vessels need to be removed.
indxSkipAll  = [];
indxVesselCen = [];
indxSkip     = zeros(6,1);
for i = 1:size(Vessel,1)
    % skip all the cells inside the vessel
    indxSkip(1) = sub2ind([length(xVector),length(yVector)],Vessel(i,1)-1,Vessel(i,2)-2);
    indxSkip(2) = sub2ind([length(xVector),length(yVector)],Vessel(i,1)+1,Vessel(i,2)-2);
    indxSkip(3) = sub2ind([length(xVector),length(yVector)],Vessel(i,1)-2,Vessel(i,2));
    indxSkip(4) = sub2ind([length(xVector),length(yVector)],Vessel(i,1)+2,Vessel(i,2));
    indxSkip(5) = sub2ind([length(xVector),length(yVector)],Vessel(i,1)-1,Vessel(i,2)+2);
    indxSkip(6) = sub2ind([length(xVector),length(yVector)],Vessel(i,1)+1,Vessel(i,2)+2);
    %     indxSkip(7) = sub2ind([length(xVector),length(yVector)],Vessel(i,1),Vessel(i,2));
    
    % The six points used to fit the vessel
    indxVessel{i}(:,1) = sub2ind([length(xVector),length(yVector)],Vessel(i,1)-3,Vessel(i,2)-1);
    indxVessel{i}(:,2) = sub2ind([length(xVector),length(yVector)],Vessel(i,1)-3,Vessel(i,2)+1);
    indxVessel{i}(:,3) = sub2ind([length(xVector),length(yVector)],Vessel(i,1),Vessel(i,2)-3);
    indxVessel{i}(:,4) = sub2ind([length(xVector),length(yVector)],Vessel(i,1),Vessel(i,2)+3);
    indxVessel{i}(:,5) = sub2ind([length(xVector),length(yVector)],Vessel(i,1)+3,Vessel(i,2)-1);
    indxVessel{i}(:,6) = sub2ind([length(xVector),length(yVector)],Vessel(i,1)+3,Vessel(i,2)+1);
    
    % These grid nodes are used for vessel fitting (6 nodes on the surface
    % of the vessel), so they are skipped.
    indxSkipAll  = [indxSkipAll;indxSkip];
    % The central points on the vessels
    indxVesselCen = [indxVesselCen; sub2ind([length(xVector),length(yVector)],Vessel(i,1),Vessel(i,2))];
end
%

% Create a large volume image with all voxels being 1 (White)
volImgRef = 255*ones(sizeImEnlarge,'uint8');


%% This section is used to control the ray cell distribution in the volume
t = 0;
s = 0;
raycellWidth = [];
keepRayCell = [];
raycellXindAll_update = [];
% raycellXindAll means the ray cell group index along y dirction
if isExistRayCell
    for i = 1:length(raycellXindAll)
        
        % generate the raycell groups along z direciton. Number of ray cells
        % clusters
        m = sizeImEnlarge(3)/rayCellNum/raywHeight+6;
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
                    ([0:rayCellNumGroup]+rayCellSpace(j)+round(-30*rand(1)))*raywHeight;
            else
                raycellWidth_temple{i}  = raycellWidth_temple{i}(end)+...
                    ([0:rayCellNumGroup]+rayCellSpace(j))*raywHeight;
            end
            
            %% Generate the ray cells in this group
            if raycellWidth_temple{i}(1)<=sizeImEnlarge(3)-150 && raycellWidth_temple{i}(end)>=150
                t = t+1;
                raycellXindAll_update(t) = raycellXindAll(i);
                raycellWidth{t}  = round(raycellWidth_temple{i});
                keepRayCell(t) = i;
                raycellXind_all{t} = raycellXind;
            end
        end
        
    end
end
% Generate small fibers.
t = 0;
SkipFibre = [];
if isExistRayCell && ~isempty(keepRayCell)
    skipFiberColumn = [raycellXindAll(squeeze(keepRayCell)),raycellXindAll(squeeze(keepRayCell))+1];
else
    skipFiberColumn = [];
end
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
                if ~isempty(find(indxSkipAll==ind))
                    % if this node exist on the surface of vessels, skip it.
                    continue
                else
                    % small fibers
                    % find the four neighbored nodes for elipse fitting
                    for k = 1:4
                        indxNeighbor    = length(xVector)*(j-1+neighborLocal(2,k))+i1+neighborLocal(1,k);
                        pointCoord(k,:) = [xGrid_all(indxNeighbor,iSlice),yGrid_all(indxNeighbor,iSlice)];
                    end
                    % make the elipse more elipse
                    if skipCellThick == 0
                        pointCoord(2,2)     = pointCoord(2,2)-2;
                        pointCoord(4,2)     = pointCoord(4,2)+2;
                    end
                    
                    C_elipse                = fitElipse(pointCoord); % Estimate the coeffecients of the elipse.
                    % Then we can estimate the diameter along two direction for the fiber
                    % The rectangle region covering the elipse.
                    regionCellindx = ceil(ceil(C_elipse(3))+[-max(floor(C_elipse(1:2))):max(floor(C_elipse(1:2)))]);
                    regionCellindy = ceil(ceil(C_elipse(4))+[-max(floor(C_elipse(1:2))):max(floor(C_elipse(1:2)))]);
                    regionCellindx = max(1,min(regionCellindx)):min(max(regionCellindx),sizeImEnlarge(1));
                    regionCellindy = max(1,min(regionCellindy)):min(max(regionCellindy),sizeImEnlarge(2));
                    [regionCellx,regionCelly] = ndgrid(regionCellindx,regionCellindy);
                    
                    % External contour of the elipse
                    inElipse1      = (regionCellx-C_elipse(3)).^2./C_elipse(1)^2 + (regionCelly-C_elipse(4)).^2./C_elipse(2)^2;
                    % Internal contour of the elipse
                    inElipse2      = (regionCellx-C_elipse(3)).^2./(C_elipse(1)-Thickness_all(ind,iSlice)-skipCellThick).^2 + ...
                        (regionCelly-C_elipse(4)).^2./(C_elipse(2)-Thickness_all(ind,iSlice)-skipCellThick).^2;
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
end
volImgRef(volImgRef>255) = 255;
if 0
    h2 = figure;
    imshow(volImgRef(:,:,1))
    set(gcf,'color','w'),axis equal,axis tight,xlabel('X (voxel)'),ylabel('Y (voxel)')
    saveas(h2,[SaveFolder,'surf1'],'epsc');
end
%% Generate large vessels.

t = 0;
for i = 2:2:length(xVector)-1
    for j = 2:2:length(yVector)-1
        t = t+1;
        
        % This is used to control the length of the cells.
        
        vesselEndLocAll    = [];
        Initial = round(rand(1)*vesselLength);
        vesselEndLocAll(1) = Initial;
        for s = 1:ceil(max((sizeImEnlarge(3)+5*vesselLength)/vesselLength+3))
            temp = min(3*vesselLength,max(100,vesselLength+randn(1)*vesselLengthVariance));
            vesselEndLocAll(s+1) = round(vesselEndLocAll(s)+temp);
        end
        
        vesselEndLoc = vesselEndLocAll-4*vesselLength+1;
        
        % The end of vessel should be inside the volume
        vesselEndLoc = vesselEndLoc(vesselEndLoc>=4 & vesselEndLoc<=sizeImEnlarge(3)-4);
        
        
        vesselEnd  = [];
        for iThick = 0:cellEndThick-1
            vesselEnd = [vesselEnd;vesselEndLoc+iThick];
        end
        
        for iSlice = 1:sizeImEnlarge(3)
            i1 = i;
            % used to store the six neighbored points for elipse fitting
            pointCoord = [];
            
            if mod(j,4) == 0
                i1 = i+1;
            end
            
            % the location index of the vessel
            ind  = length(xVector)*(j-1)+i1;
            
            iVessel = find(ind==indxVesselCen);
            if ~isempty(iVessel)
                % large vessels. Fit six points to get the boundary of the
                % vessel
                pointCoord = [xGrid_all([indxVessel{iVessel}]',iSlice),yGrid_all([indxVessel{iVessel}]',iSlice)];
                C_elipse   = fitElipse(pointCoord);
                
                % Fit the surface for the elipse
                for t1 = round(C_elipse(3))+[-max(floor(C_elipse(1:2))):max(floor(C_elipse(1:2)))]
                    for t2 = round(C_elipse(4))+[-max(floor(C_elipse(1:2))):max(floor(C_elipse(1:2)))]
                        % Make sure the point is in the image
                        if t1>0 && t1<sizeImEnlarge(1) && t2>0 && t2<sizeImEnlarge(2)
                            
                            inElipse1 = (t1-C_elipse(3))^2/C_elipse(1)^2 + (t2-C_elipse(4))^2/C_elipse(2)^2;
                            inElipse2 = (t1-C_elipse(3))^2/(C_elipse(1)-Thickness_all(ind,iSlice)-vesselThicker)^2 + ...
                                (t2-C_elipse(4))^2/(C_elipse(2)-Thickness_all(ind,iSlice)-vesselThicker)^2;
                            
                            if inElipse1<=1 && inElipse2>=1
                                % on the vessel wall, all points should be
                                % white
                                volImgRef(t1,t2,iSlice) = 255;
                            else
                                % inside the vessel, all points should be
                                % black
                                
                                if inElipse2 < 1
                                    if isempty(find(iSlice==vesselEnd))
                                        volImgRef(t1,t2,iSlice) = 0;
                                    else
                                        volImgRef(t1,t2,iSlice) = 255;
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

if 0
    h3 = figure;
    imshow(volImgRef(:,:,1))
    set(gcf,'color','w'),axis equal,axis tight,
    saveas(h3,[SaveFolder,'slice1_inital'],'epsc');
end


% Folder_Init = [SaveFolder,'/volImgInit/'];
% mkdir(Folder_Init);
parfor z = 1:size(volImgRef,3)
    imgZ = volImgRef(:,:,z);
    zStr = num2str(z);
    numStr = '00000';
    numStr(end-length(zStr)+1:end) = zStr;
    FileNameAll{z} = ['volImgRef_',numStr,'.tif'];
    FileUAll{z} = ['u_volImgRef_',numStr,'.csv'];
    FileVAll{z} = ['v_volImgRef_',numStr,'.csv'];
    FileVRayAll{z} = ['v_ray_volImgRef_',numStr,'.csv'];
    FileWAll{z} = ['w_volImgRef_',numStr,'.csv'];
    
    if 0
        imgName = fullfile(Folder_Init,FileNameAll{z});
        imwrite(imgZ,imgName);
    end
end
volImgRef_final = volImgRef;
clear volImgRef;
t = 0;
%% Generate ray cells. Ray cells are inserted into the volume by calling function rayCellGenerate
% raycellXindAll means the ray cell group index along y dirction
if isExistRayCell
    for t = 1:length(raycellWidth)
        volImgRef_final    = rayCellGenerate(volImgRef_final,raycellXind_all{t},raycellWidth{t},Params);
    end
end

if 1
    h4 = figure;
    imshow(squeeze(volImgRef_final(240,:,:)))
    set(gcf,'color','w'),axis equal,axis tight,
    saveas(h4,[SaveFolder,'slice1_raycell'],'epsc');
end
% save the volume image

%% Add complicated deformation to the volume image
% The deformation fields are generated seperatedly. Then, they are summed
% together.
% Here u, v are initialized to be zero. Then they are summed.
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
        if ~isempty(find(indxSkipAll==ind))
            continue
        else
            iVessel = find(ind==indxVesselCen);
            if  isempty(iVessel)
                % Small fibers
                islarge    = 0;
                indxPt     = length(xVector)*(j-1)+i1;
                pointCoord = [xGrid_all(indxPt,1),yGrid_all(indxPt,1)];
                
                %----------the value here are super important-------------%
                k          = [0.08,0.06,2+rand(1),1*(1+rand(1))];
                u_temp     = localDistort(x,y,pointCoord(1),pointCoord(2),k);
                %----------the value here are super important-------------%
                k          = [0.08,0.06,2+rand(1),1*(1+rand(1))];
                v_temp     = localDistort(y,x,pointCoord(2),pointCoord(1),k);
                signTemp   = sign(randn(1));
                u          = u - signTemp*u_temp;
                v          = v - signTemp*v_temp;
                
            else
                
                if 1
                    % Large vessels
                    islarge    = 1;
                    indxPt     = length(xVector)*(j-1)+i1;
                    pointCoord = [xGrid_all(indxPt,1),yGrid_all(indxPt,1)];
                    %----------the value here are super important-------------%
                    k          = [0.04,0.03,1*(1+rand(1)),2+1*rand(1)];
                    u_temp     = localDistort(x,y,pointCoord(1),pointCoord(2),k);
                    v_temp     = localDistort(y,x,pointCoord(2),pointCoord(1),k);
                    
                    signDisp   = sign(randn(1));
                    u          = u - signDisp* u_temp;
                    %                 v          = v - signDisp* v_temp;
                end
            end
        end
        if rand(1)<1/100 % here just means the possibility.
            % If it is smaller, that means less distortion will be added
            % Deform for a large regions, even larger than a vessel
            islarge    = 1;
            indxPt     = length(xVector)*(j-1)+i1;
            pointCoord = [xGrid_all(indxPt,1),yGrid_all(indxPt,1)];
            
            if randn(1)>0
                %----------the value here are super important-------------%
                k          = [0.01,0.008,1.5*(1+rand(1)),1*(1+1*rand(1))];
                u_temp     = localDistort(x,y,pointCoord(1),pointCoord(2),k);
                
                u1          = u1 + sign(randn(1))*u_temp;
            else
                %----------the value here are super important-------------%
                k          = [0.01,0.008,1.5+1.5*rand(1),1*(1+1*rand(1))];
                v_temp     = localDistort(y,x,pointCoord(2),pointCoord(1),k);
                v1          = v1 + sign(randn(1))*v_temp;
            end
            
        end
    end
end

if 0
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


% Create te folder to save the results

Folder_ImgSeries_LocalDist = [SaveFolder,'/LocalDistVolume/'];
mkdir(Folder_ImgSeries_LocalDist)
if writeLocalDeformData
    Folder_dips_LocalDistU = [SaveFolder,'/LocalDistVolumeDispU/'];
    mkdir(Folder_dips_LocalDistU)
    Folder_dips_LocalDistV = [SaveFolder,'/LocalDistVolumeDispV/'];
    mkdir(Folder_dips_LocalDistV)
    Folder_dips_LocalDistVRay = [SaveFolder,'/LocalDistVolumeDispVRay/'];
    mkdir(Folder_dips_LocalDistVRay)
end
%% Impose local distortion by image interpolation
parfor (z = 1:sizeImEnlarge(3),4)
    imgName = fullfile(Folder_ImgSeries,FileNameAll{z});
    ImCover = double(imread(imgName));
    [x_grid,y_grid]      = ndgrid(1:sizeImEnlarge(1),1:sizeImEnlarge(2));
    if isExistRayCell
        v_all_raycell = rayCell_shrinking(...
            Params,yVector,sizeImEnlarge,raycellWidth,raycellXindAll_update,v,raywHeight,cellR-cellWallThick/2,z);
        x_interp = x_grid+u;
        y_interp = y_grid+v+v_all_raycell;
    else
        
        x_interp = x_grid+u;
        y_interp = y_grid+v;
    end
    F       = scatteredInterpolant(x_interp(:),y_interp(:),ImCover(:),'linear');
    Vq = F(x_grid(:),y_grid(:));
    
    imInterp = reshape(Vq,sizeImEnlarge(1:2));
    imInterp(isnan(imInterp)) = 0;
    imInterp(imInterp<127)    = 0;
    imInterp(imInterp>=127)   = 255;
    ImgTemp                   = uint8(imInterp);
    
    if writeLocalDeformData
        u_FileName = fullfile(Folder_dips_LocalDistU,FileUAll{z});
        writematrix(u,u_FileName)
        v_FileName = fullfile(Folder_dips_LocalDistV,FileVAll{z});
        writematrix(v,v_FileName)
        
        v_ray_FileName = fullfile(Folder_dips_LocalDistVRay,FileVRayAll{z});
        writematrix(v_all_raycell,v_ray_FileName)
    end
    imgName = fullfile(Folder_ImgSeries_LocalDist,FileNameAll{z});
    imwrite(ImgTemp,imgName);
end

%% image interpolation

Folder_ImgSeries_GlobalDist = [SaveFolder,'/GlobalDistVolume/'];
mkdir(Folder_ImgSeries_GlobalDist)
if writeGlobalDeformData
    Folder_dips_GlobalDistU = [SaveFolder,'/GlobalDistVolumeDispU/'];
    mkdir(Folder_dips_GlobalDistU)
    Folder_dips_GlobalDistV = [SaveFolder,'/GlobalDistVolumeDispV/'];
    mkdir(Folder_dips_GlobalDistV)
end

%% Deform several slices simultaneous using global distortion
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
    
    %         v_allz   = (x_grid-sizeImEnlarge(1)/3).^2/1e4/8+(y_grid-sizeImEnlarge(2)/2).*(z_grid-sizeImEnlarge(3)/2)/(10^4)/9;
    %         u_allz   = (y_grid-sizeImEnlarge(1)/3).^2/1e4/9+(x_grid-sizeImEnlarge(2)/3).*(z_grid-sizeImEnlarge(3)/3)/(10^4)/7;
    v_allz   = (x_grid-sizeImEnlarge(1)/3).^2/1e4/8+(y_grid-sizeImEnlarge(2)/2).*(x_grid-sizeImEnlarge(1)/3)/1e4/9;
    u_allz   = (y_grid-sizeImEnlarge(1)/3).^2/1e4/9+(x_grid-sizeImEnlarge(2)/2).*(y_grid-sizeImEnlarge(1)/3)/1e4/8;
    
    %         w_allz   = (z_grid-sizeImEnlarge(3)/2).*(x_grid-sizeImEnlarge(1)/2)/1e4+5*sin(pi*x_grid/(sizeImEnlarge(1)/2)+z_grid/(sizeImEnlarge(3)/2));
    
    u_allz_temp = -u_allz-repmat(u1,1,1,size(u_allz,3));
    v_allz_temp = -v_allz-repmat(v1,1,1,size(u_allz,3));
    
    %     u_allz_temp = -repmat(u1,1,1,size(x_grid,3));
    %     v_allz_temp = -repmat(v1,1,1,size(x_grid,3));
    
    x_interp = x_grid+u_allz_temp;
    y_interp = y_grid+v_allz_temp;
    z_interp = z_grid;
    
    if t == 1
        if 1
            h6 = figure;
            surf(u_allz_temp(end:-1:1,:,1)),hc = colorbar; view([0,90]);
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
        
        if writeGlobalDeformData
            disp_u = u_allz_temp(:,:,k);
            dispUName = fullfile(Folder_dips_GlobalDistU,FileUAll{z});
            writematrix(disp_u,dispUName)
            
            disp_v = v_allz_temp(:,:,k);
            dispVName = fullfile(Folder_dips_GlobalDistV,FileVAll{z});
            writematrix(disp_v,dispVName)
        end
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
