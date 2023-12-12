function volImgRef_final = rayCellGenerate(volImgRef_final,raycellXind,raycellWidth,Params)
% function to generate ray cells
rayCellLength  = Params.rayCellLength;
rayCell_variance  = Params.rayCell_variance;
rayHeight       = Params.rayHeight;
xGrid_all      = Params.xGrid_all;
yGrid_all      = Params.yGrid_all;
Thickness_all  = Params.Thickness_all;
sizeImEnlarge  = Params.sizeImEnlarge;
xGrid          = Params.xGrid;
cellEndThick   = Params.cellEndThick;
cellThick      = Params.cellThick;

% raycellWidth(raycellWidth<-rayHeight | raycellWidth>sizeImEnlarge(3)+rayHeight) = [];
for i_raycolumn = raycellXind
    % between first column and second column need a deviation
    raycolumn_rand = 1/2*rayHeight;
    s = 1;
    vesselEndLoc = [];
    vesselEndLoc(1) = -rayCellLength*rand(1);
    while vesselEndLoc(s)<sizeImEnlarge(1)+rayCellLength
        s = s+1;
        temp = rayCellLength + rayCell_variance*randn(1);
        if temp<rayCellLength/3
            temp=2*rayCellLength;
        else
            if temp>2*rayCellLength
                temp=2*rayCellLength;
            end
        end
        vesselEndLoc(s) = vesselEndLoc(s-1)+temp;
    end
    vesselEndLoc = round(vesselEndLoc);
%     vesselEndLoc = [1:rayCellLength:sizeImEnlarge(1)-rayCellLength/3]';
%     vesselEndLoc = vesselEndLoc+round(rayCell_variance*rand(length(vesselEndLoc),1));
    vesselEndLoc(vesselEndLoc>sizeImEnlarge(1)+rayCellLength/2) = [];
%     vesselEndLoc1 = vesselEndLoc;
    % for each column of ray cells, 
%     endRandom    = round(rayCell_variance*1/4+rayCell_variance/2*rand(length(raycellWidth)));
    for i_rayrow = 1:length(vesselEndLoc)-1
        m2 = 0;
        for j_slice = round(raycellWidth)
            m2 = m2+1;
            
            if mod(i_raycolumn,2) == 0
                k = round(j_slice+round(raycolumn_rand));
            else
                k = j_slice;
            end
            t              = round((max(1,k)+min(k+rayHeight,sizeImEnlarge(3)))/2);
            if t >= 1 && t<sizeImEnlarge(3)-10
            xGrid_t        = reshape(xGrid_all(:,t),size(xGrid));
            yGrid_t        = reshape(yGrid_all(:,t),size(xGrid));
            thickGrid_t    = reshape(Thickness_all(:,t),size(xGrid));

            yInterp1_C     = spline(xGrid_t(:,i_raycolumn),yGrid_t(:,i_raycolumn),1:sizeImEnlarge(1))-1.5; 
            thickInterp_C  = spline(xGrid_t(:,i_raycolumn),thickGrid_t(:,i_raycolumn),1:sizeImEnlarge(1));

            yInterp2_C     = spline(xGrid_t(:,i_raycolumn+1),yGrid_t(:,i_raycolumn+1),1:sizeImEnlarge(1))+1.5; 

            cell_center    = [(1:length(yInterp2_C))',round(yInterp2_C(:)+yInterp1_C(:))/2,...
                    repmat((max(1,k)+min(k+rayHeight,sizeImEnlarge(3)))/2,length(yInterp2_C),1)];
            cell_r         = [(yInterp2_C(:)-yInterp1_C(:))/2,repmat((min(k+rayHeight,sizeImEnlarge(3))-max(1,k))/2,length(yInterp2_C),1)];
            
            vesselEndLoc_column = vesselEndLoc+round(mod(m2,2)*rayCellLength/2);
            if vesselEndLoc_column(i_rayrow+1)<=sizeImEnlarge(1)-1 && vesselEndLoc_column(i_rayrow)>= 2

                

                cell_neighPt   = [vesselEndLoc_column(i_rayrow),vesselEndLoc_column(i_rayrow+1);
                                   round(yInterp1_C(vesselEndLoc_column(i_rayrow))),...
                                   round(yInterp2_C(vesselEndLoc_column(i_rayrow)));...
                                   max(1,k),min(k+rayHeight,sizeImEnlarge(3))];
                
                if i_rayrow == 1
                    i_valid = vesselEndLoc_column(i_rayrow)+cellEndThick:vesselEndLoc_column(i_rayrow+1)-cellEndThick/2;
                    
                else
                    if i_rayrow == length(vesselEndLoc)-1
                        i_valid = vesselEndLoc_column(i_rayrow)+cellEndThick/2-1:vesselEndLoc_column(i_rayrow+1)-cellEndThick;
                    else
                        i_valid = vesselEndLoc_column(i_rayrow)+cellEndThick/2-1:vesselEndLoc_column(i_rayrow+1)-cellEndThick/2;
                    end
                end
                

                for i = vesselEndLoc_column(i_rayrow)+1:vesselEndLoc_column(i_rayrow+1)
                    
                    if j_slice == min(raycellWidth)
                    volImgRef_final(i,...
                                yInterp1_C(i):yInterp2_C(i),...
                                cell_center(i,3):cell_neighPt(3,2)) = 255; 
                    else
                        if j_slice == max(raycellWidth)
                            volImgRef_final(i,...
                                        yInterp1_C(i):yInterp2_C(i),...
                                        cell_neighPt(3,1):cell_center(i,3)) = 255; 
                        else
                            volImgRef_final(i,...
                                        yInterp1_C(i):yInterp2_C(i),...
                                        cell_neighPt(3,1):cell_neighPt(3,2)) = 255; 
                        end
                    end
                    
                    if find(round(i_valid) == i)
                        for j = cell_neighPt(2,1):cell_neighPt(2,2)
                            for s = cell_neighPt(3,1):cell_neighPt(3,2)
                                if (j-cell_center(i,2))^2/(cell_r(i,1)-thickInterp_C(i))^2 ...
                                            + (s-cell_center(i,3))^2/(cell_r(i,2)-thickInterp_C(i))^2 <1
                                            volImgRef_final(i,j,s) = 0;
                                else
                                    if ((j-cell_center(i,2))^2/(cell_r(i,1)-thickInterp_C(i)*2/3)^2 ...
                                            + (s-cell_center(i,3))^2/(cell_r(i,2)-thickInterp_C(i))^2) >=1 ...
                                        && ((j-cell_center(i,2))^2/(cell_r(i,1))^2 ...
                                            + (s-cell_center(i,3))^2/(cell_r(i,2))^2) <=1
                                        volImgRef_final(i,j,s) = 255;
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
end

