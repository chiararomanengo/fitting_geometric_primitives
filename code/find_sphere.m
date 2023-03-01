function rb = find_sphere(r,XY)

% Parameter space for centers and radius:
step_size=r/50;
rc = 0.95*r:(1.05*r-0.95*r)/50:1.05*r;
H = zeros(length(rc),1);

% Optimal parameters
for i=1:length(rc)
    cur_sph = (XY(:,1)).^2+(XY(:,2)).^2+(XY(:,3)).^2-rc(i)^2;
    H(i) = sum(abs(cur_sph)<step_size/2);
end


[max_H,idx_max_H] = max(H);

rb = rc(idx_max_H);


end