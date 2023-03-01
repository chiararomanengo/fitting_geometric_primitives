function [coeffImpl, mfe]=cone(A,xyzAus,t)

    D=0;
    n=0;
    
    PtMin =0.85*A; %
    PtMax =1.15*A; %
   
    Eps = (1/20)*(PtMax(1)-PtMin(1));

    [Coord, Max,H]=extract_superficie('Cono', xyzAus, PtMin, PtMax, Eps,D,n,t);

    if Max>0
        mfe=10;
        for j=1:numel(Coord)
            a=Coord(j);
            
            theta=2*pi;
            N=500;
            h2=2*max([max(xyzAus(:,3)) abs(min(xyzAus(:,3)))]);
            cono = zeros(N^2,3);
            for i=1:size( cono,1)
                u = rand(1)*theta;
                v = -h2+rand(1)*(h2*2);
                cono(i,1)=sin(a)*cos(u)*v;
                cono(i,2)=sin(a)*sin(u)*v;
                cono(i,3)=cos(a)*v;
            end
            
            [I,dist] = knnsearch(cono,xyzAus);
            
            m=MFE(xyzAus,dist);
            if m<mfe
                mfe=m;
                aFin=a;
            end
            coeffImpl=aFin;
        end      
    else
        coeffImpl=[];
        mfe=NaN;
    end
end