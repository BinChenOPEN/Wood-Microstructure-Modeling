function [v_all] = rayCell_shrinking_spruce_20220322(...
    Params,yVector,sizeImEnlarge,raycellWidth,raycellXindAll,raywidth,raysize,indxZ)

xGrid_all = Params.xGrid_all;
yGrid_all = Params.yGrid_all;
Thickness_all = Params.Thickness_all;

xLeft = xGrid_all(:,indxZ);
yLeft = yGrid_all(:,indxZ);
Thickness_slice = Thickness_all(:,indxZ);
x_node_grid = reshape(xLeft,Params.gridSize);
y_node_grid = reshape(yLeft,Params.gridSize);
thick_node_grid = reshape(Thickness_slice,Params.gridSize);

[xgrid,ygrid]  = ndgrid(1:sizeImEnlarge(1),1:sizeImEnlarge(2));


for i = 1:length(raycellXindAll)
    v_all_temp  = zeros(sizeImEnlarge(1),sizeImEnlarge(2),length(indxZ));
%     ycenter = (yVector(raycellXindAll(i))+yVector(raycellXindAll(i)+1))/2;
    y_node_grid1 = interp1(x_node_grid(:,raycellXindAll(i)),...
        y_node_grid(:,raycellXindAll(i)),1:sizeImEnlarge(1),'spline')';
    y_node_grid2 = interp1(x_node_grid(:,raycellXindAll(i)+2),...
        y_node_grid(:,raycellXindAll(i)+2),1:sizeImEnlarge(1),'spline')';
    R            = (y_node_grid2-y_node_grid1)./2-2;
    ycenter      = (y_node_grid2+y_node_grid1)./2;
    
    v1 = zeros(sizeImEnlarge(1:2));
    for j = 1:sizeImEnlarge(1)
        v1(j,1:ycenter(j)-R(j)) = R(j);
        v1(j,ycenter(j)+R(j):end) = -R(j);
        v1(j,ycenter(j)-R(j) :ycenter(j)+R(j)) = ...
            -(ygrid(j,ycenter(j)-R(j):ycenter(j)+R(j))-ycenter(j));
    end
    
    
    
    minIndx1 = max(1,round(min(raycellWidth{i})+1*raysize));
    maxIndx1 = min(round(max(raycellWidth{i})+raywidth-1*raysize),sizeImEnlarge(3));
    if maxIndx1-minIndx1>10
    indx1 = minIndx1:maxIndx1;
    else
        indx1 = [];
    end
    
    minIndx0 = max(1,round(min(raycellWidth{i})-1*raysize));
    maxIndx0 = min(round(max(raycellWidth{i})+raywidth+1*raysize),sizeImEnlarge(3));
    if maxIndx0-minIndx0>10
    indx0 = minIndx0:maxIndx0;
    end
    indx2 = 1:sizeImEnlarge(3);
    indx2(indx0) = [];
    
    minIndx3 = max(1,round(min(raycellWidth{i})-3*raysize));
    maxIndx3 = min(round(min(raycellWidth{i})+1*raysize),sizeImEnlarge(3));
    if maxIndx3 - minIndx3 >10
    indx3 = minIndx3: maxIndx3;
    else
        indx3 = [];
    end
    
    minIndx4 = max(1,round(max(raycellWidth{i})+raywidth-1*raysize));
    maxIndx4 = min(round(max(raycellWidth{i})+raywidth+3*raysize),sizeImEnlarge(3));
    if maxIndx4 - minIndx4 >10
    indx4 = minIndx4: maxIndx4;
    else
        indx4 = [];
    end
    
    coeff1 = zeros(1,sizeImEnlarge(3));
    coeff2 = zeros(1,sizeImEnlarge(3));
    
    coeff1(indx1) = 0;
    coeff1(indx2) = 1;
    
    coeff2(indx1) = 1;
    coeff2(indx2) = 0;
    
    if length(indx3) ~= 0
        coeff1(indx3) = [length(indx3)-1:-1:0]./round(4*raysize);
        coeff2(indx3) = [0:length(indx3)-1]./round(4*raysize);
    end
    if length(indx4) ~= 0
        coeff1(indx4) = [0:length(indx4)-1]./round(4*raysize);
        coeff2(indx4) = [length(indx4)-1:-1:0]./round(4*raysize);

    end
%     coeff2(indx3) = [0:round(4*raysize)]./round(4*raysize);
%     coeff2(indx4) = [round(4*raysize):-1:0]./round(4*raysize);
    

    v2 = zeros(sizeImEnlarge(1:2));
    temp2_1 = v2;
    temp2_2 = v2;
    for j  = 1:sizeImEnlarge(1) 
        temp2_1(j,:) = R(j)*(1./(1+exp(0.01.*(ygrid(j,:)-ycenter(j)+R(j))))-0.5);
        temp2_2(j,:) = R(j)*(1./(1+exp(0.01.*(ygrid(j,:)-ycenter(j)-R(j))))-0.5);
        v2(j,1:ycenter(j)-R(j))   = temp2_1(j,1:ycenter(j)-R(j)) ;
        v2(j,ycenter(j)+R(j):end) = temp2_2(j,ycenter(j)+R(j):end);
        
        v2(j,1:ycenter(j)-R(j))            = v2(j,1:ycenter(j)-R(j)) + R(j)/2;
        v2(j,ycenter(j)+R(j):end)          = v2(j,ycenter(j)+R(j):end)-R(j)/2;
        v2(j,ycenter(j)-R(j) :ycenter(j)+R(j)) = ...
            -(ygrid(j,ycenter(j)-R(j):ycenter(j)+R(j))-ycenter(j))/2;
    end
    
        
    for j = 1:length(indxZ)
        v_all_temp(:,:,j) = v_all_temp(:,:,j)+coeff1(indxZ(j)).*v1+coeff2(indxZ(j)).*v2;
    end
    v_all_structure{i} = v_all_temp;
end
v_all  = zeros(sizeImEnlarge(1),sizeImEnlarge(2),length(indxZ));
for i = 1:length(raycellXindAll)
    v_all = v_all + v_all_structure{i};
end