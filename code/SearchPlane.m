function [output, mfe] = SearchPlane(xyz)

try
    
[coord, maxCoord]=PlanesHT0_search(xyz);

mAus=1;
for j=1:size(coord,1)
    rho=coord(j,1);
    theta=coord(j,2);
    phi=coord(j,3);
    if abs(cos(phi))>10^(-3)
        x_piano = linspace(min(xyz(:,1)),max(xyz(:,1)),800);
        y_piano=linspace(min(xyz(:,2)),max(xyz(:,2)),800);
        piano = zeros(length(x_piano)*length(y_piano),3);
        for i=1:length(x_piano)
            for k=1:length(y_piano)
                piano( (i-1)*length(x_piano)+k,1) = x_piano(i);
                piano( (i-1)*length(x_piano)+k,2) = y_piano(k);
                piano( (i-1)*length(x_piano)+k,3) = (rho-x_piano(i)*cos(theta)*sin(phi)-y_piano(k)*sin(phi)*sin(theta))/cos(phi);
            end
        end
    else
        if abs(cos(theta))>10^(-3)
            y_piano=linspace(min(xyz(:,2)),max(xyz(:,2)),800);
            z_piano=linspace(min(xyz(:,3)),max(xyz(:,3)),800);
            piano = zeros(length(y_piano)*length(z_piano),3);
            for i=1:length(y_piano)
                for k=1:length(z_piano)
                    piano( (i-1)*length(y_piano)+k,1) = (rho-y_piano(i)*sin(theta)*sin(phi))/(sin(phi)*cos(theta));
                    piano( (i-1)*length(y_piano)+k,2) = y_piano(i);
                    piano( (i-1)*length(y_piano)+k,3) = z_piano(k);
                end
            end
        else
            x_piano=linspace(min(xyz(:,1)),max(xyz(:,1)),800);
            z_piano=linspace(min(xyz(:,3)),max(xyz(:,3)),800);
            piano = zeros(length(x_piano)*length(z_piano),3);
            for i=1:length(x_piano)
                for k=1:length(z_piano)
                    piano( (i-1)*length(x_piano)+k,1) = x_piano(i);
                    piano( (i-1)*length(x_piano)+k,2) = (rho)/(sin(phi)*sin(theta));
                    piano( (i-1)*length(x_piano)+k,3) = z_piano(k);
                end
            end
        end
    end   
    
    [I,dist] = knnsearch(piano,xyz);
    m=MFE(xyz,dist);
    
    if m<mAus
        mAus=m;
        mfe=m;
        output=[cos(theta)*sin(phi) sin(phi)*sin(theta) cos(phi) -rho];
    end   
end
catch
    mfe=NaN;
    output=[];   
    return
end


end