function fisher_n = single_neuron_FIM(curve2_n)
%
% calculate fisher information for single neurons 2-D tuning curve

% construct the fisher information matrix (equation 29)
[FX,FY] = gradient(curve2_n);

for dim1 = 1:2
    for dim2 = 1:2
        
        % crop response matrix and calculate approx partial first derivative along specified dimension
        if dim1 == 1
            d_dim1 = FY;
        elseif dim1 == 2
            d_dim1 = FX;
        end
        
        if dim2 == 1
            d_dim2 = FY;
        elseif dim2 == 2
            d_dim2 = FX;
        end
        
        fisher_n(:,:,dim1,dim2) = d_dim1.*d_dim2.*(1./curve2_n);
        
        % remove NaNs
        fisher_n(isnan(fisher_n)) = 0;
        
    end
end

% if(0)
%     for dim1 = 1:2
%         for dim2 = 1:2
%             
%             % crop response matrix and calculate approx partial first derivative along specified dimension
%             if dim1 == 1
%                 d_dim1 = diff(curve2_n(:,1:end-1),1,1);
%             elseif dim1 == 2
%                 d_dim1 = diff(curve2_n(1:end-1,:),1,2);
%             end
%             
%             if dim2 == 1
%                 d_dim2 = diff(curve2_n(:,1:end-1),1,1);
%             elseif dim2 == 2
%                 d_dim2 = diff(curve2_n(1:end-1,:),1,2);
%             end
%             
%             fisher_n(:,:,dim1,dim2) = d_dim1.*d_dim2.*(1./curve2_n(1:end-1,1:end-1));
%             
%             % remove NaNs
%             fisher_n(isnan(fisher_n)) = 0;
%             
%         end
%     end
% end


