function [xb, yb, rb, max_H] = find_circle(XY, xc_min, xc_max, yc_min, yc_max, r_min, r_max)

% Parameter space for centers and radius:
step_sizeX=(xc_max-xc_min)/50;
xc = xc_min:step_sizeX:xc_max;
step_sizeY=(yc_max-yc_min)/50;
yc = yc_min:step_sizeY:yc_max;
step_size=(r_max-r_min)/50;
rc = r_min:step_size:r_max;
H = zeros(length(xc), length(yc), length(rc));

% Optimal parameters
for i=1:length(xc)
    for j=1:length(yc)
        for k=1:length(rc)
            cur_sph = (XY(:,1) - xc(i)).^2+(XY(:,2) - yc(j)).^2 - rc(k)^2;
            H(i,j,k) = sum(abs(cur_sph)<step_size/2);
        end
    end
end


[max_H,idx_max_H] = max(H(:));
[I1,I2,I3] = ind2sub(size(H),idx_max_H);

xb = xc(I1);
yb = yc(I2);
rb = rc(I3);



disp("");

end