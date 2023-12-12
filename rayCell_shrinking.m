function [v_all] = rayCell_shrinking(...
    Params,yVector,sizeImEnlarge,raycellWidth,raycellXindAll,v,raywidth,raysize,indxZ)
% ThiS function is used to compress the ray cell regions
xGrid_all = Params.xGrid_all;
yGrid_all = Params.yGrid_all;
Thickness_all = Params.Thickness_all;
cellThick      = Params.cellThick;

xLeft = xGrid_all(:,indxZ);
yLeft = yGrid_all(:,indxZ);
Thickness_slice = Thickness_all(:,indxZ);
x_node_grid = reshape(xLeft,Params.gridSize);
y_node_grid = reshape(yLeft,Params.gridSize);
thick_node_grid = reshape(Thickness_slice,Params.gridSize);

[xgrid,ygrid]  = ndgrid(1:sizeImEnlarge(1),1:sizeImEnlarge(2));

raycellXindUnique = unique(raycellXindAll);
% raycellXindUniqueNum = hist(raycellXindAll,unique(raycellXindAll));



k = 0;
for i = 1:length(raycellXindUnique)
    coeff1 = ones(1,sizeImEnlarge(3));
    coeff2 = zeros(1,sizeImEnlarge(3));
    v1 = zeros(sizeImEnlarge(1:2));
    v2 = zeros(sizeImEnlarge(1:2));
%     for t = 1:raycellXindUnique(i)
        k = k+1;
        v_all_temp  = zeros(sizeImEnlarge(1),sizeImEnlarge(2),length(indxZ));
        
        %     ycenter = (yVector(raycellXindAll(i))+yVector(raycellXindAll(i)+1))/2;
        y_node_grid1 = interp1(x_node_grid(:,raycellXindAll(k)),...
            y_node_grid(:,raycellXindAll(k)),1:sizeImEnlarge(1),'spline')';
        y_node_grid2 = interp1(x_node_grid(:,raycellXindAll(k)+2),...
            y_node_grid(:,raycellXindAll(k)+2),1:sizeImEnlarge(1),'spline')';
        
        
        indv1 = sub2ind(size(v),[1:sizeImEnlarge(1)]',round(y_node_grid1));
        indv2 = sub2ind(size(v),[1:sizeImEnlarge(1)]',round(y_node_grid2));
        
        v_node_grid1 = v(indv1);
        v_node_grid2 = v(indv2);
        
        %     v_node_grid1 = interp1(x_node_grid(:,raycellXindAll(i)),...
        %         v_local1,1:sizeImEnlarge(1),'spline')';
        %     v_node_grid2 = interp1(x_node_grid(:,raycellXindAll(i)+2),...
        %         v_local2,1:sizeImEnlarge(1),'spline')';
        %
        R            = (y_node_grid2+v_node_grid2-y_node_grid1-v_node_grid1)./2-cellThick/2+0.5;
        ycenter      = (y_node_grid2+y_node_grid1)./2;
        
        
        for j = 1:sizeImEnlarge(1)
            indx = max(cellThick/2,round(ycenter(j)-R(j))) :min(ycenter(j)+R(j),sizeImEnlarge(2)-cellThick/2+1);
            if length(indx)>2
                v1(j,indx) = ...
                    -(indx-ycenter(j));
                
                v1(j,1:indx(1)-1) = R(j);
                v1(j,indx(end)+1:end) = -R(j);
            end
        end
        
        
        minIndx1 = max(1,round(min(raycellWidth{k})+1*raysize));
        maxIndx1 = min(round(max(raycellWidth{k})+raywidth-1*raysize),sizeImEnlarge(3));
        if maxIndx1-minIndx1>10
            indx1 = minIndx1:maxIndx1;
        else
            indx1 = [];
        end
        
%         minIndx0 = max(1,round(min(raycellWidth{k})-1*raysize));
%         maxIndx0 = min(round(max(raycellWidth{k})+raywidth+1*raysize),sizeImEnlarge(3));
%         if maxIndx0-minIndx0>10
%             indx0 = minIndx0:maxIndx0;
%         end
%         indx2 = 1:sizeImEnlarge(3);
%         indx2(indx0) = [];
        
        minIndx3 = max(1,round(min(raycellWidth{k})-3*raysize));
        maxIndx3 = min(round(min(raycellWidth{k})+1*raysize),sizeImEnlarge(3));
        if maxIndx3 - minIndx3 >10
            indx3 = minIndx3: maxIndx3;
        else
            indx3 = [];
        end
        
        minIndx4 = max(1,round(max(raycellWidth{k})+raywidth-1*raysize));
        maxIndx4 = min(round(max(raycellWidth{k})+raywidth+3*raysize),sizeImEnlarge(3));
        if maxIndx4 - minIndx4 >10
            indx4 = minIndx4: maxIndx4;
        else
            indx4 = [];
        end
        

        
        coeff1(indx1) = 0;        
        coeff2(indx1) = 1;
        
        if length(indx3) >5
            coeff1(indx3) = [length(indx3)-1:-1:0]./round(4*raysize);
            coeff2(indx3) = [0:length(indx3)-1]./round(4*raysize);
        end
        if length(indx4) >5
            coeff1(indx4) = [0:length(indx4)-1]./round(4*raysize);
            coeff2(indx4) = [length(indx4)-1:-1:0]./round(4*raysize);
        end
        
        
        temp2_1 = v2;
        temp2_2 = v2;
        for j  = 1:sizeImEnlarge(1)
            temp2_1(j,:) = 2*R(j)*(1./(1+exp(0.02.*(ygrid(j,:)-ycenter(j)+R(j))))-0.5);
            temp2_2(j,:) = 2*R(j)*(1./(1+exp(0.02.*(ygrid(j,:)-ycenter(j)-R(j))))-0.5);
            v2(j,1:ycenter(j)-R(j)) = temp2_1(j,1:ycenter(j)-R(j)) ;
            v2(j,ycenter(j)+R(j):end) = temp2_2(j,ycenter(j)+R(j):end);
        end
%     end
    
    for j = 1:length(indxZ)
        v_all_temp(:,:,j) = v_all_temp(:,:,j)+coeff1(indxZ(j)).*v1+coeff2(indxZ(j)).*v2;
    end
    u_all_structure{i} = v_all_temp;
end
v_all  = zeros(sizeImEnlarge(1),sizeImEnlarge(2),length(indxZ));
for i = 1:length(raycellXindUnique)
    v_all = v_all + u_all_structure{i};
end