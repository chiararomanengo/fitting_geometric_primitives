function [coeffImpl, mfe, Max]=cylinderZ(A,B,Xi)

    D=0;
    n=0;
    t=0;
    
    PtMin = [0.9*A, 0.9*B]; %
    PtMax = [1.1*A, 1.1*B];
    Eps = [(1/10)*(PtMax(1)-PtMin(1)) (1/10)*(PtMax(2)-PtMin(2))];

    [Coord, Max,H]=extract_superficie('cilindroZ', Xi(:,1:2), PtMin, PtMax, Eps,D,n,t);

    if Max>0 
        mfe=10;
        for j=1:numel(Coord)/2
            a=Coord(j);
            b= Coord(numel(Coord)/2+j);
            if abs(a-A)<abs(A-b)
                a=Coord(j);
            else
                a=Coord(numel(Coord)/2+j);
            end
            b=a;
            
            theta=0:pi/500:2*pi;
            z=min(Xi(:,3)):(max(Xi(:,3))-min(Xi(:,3)))/(numel(theta)-1):max(Xi(:,3));
            if numel(z)==0
                coeffImpl=[];
                mfe=NaN;
            else
               cilindro1=zeros(numel(theta)*numel(theta),3);
                theta=repmat(theta,1,numel(theta));
                z=repmat(z,numel(z),1);
                z=reshape(z,[1,numel(theta)]);
                cilindro1(:,3)=z;
                cilindro1(:,1)=a*cos(theta);
                cilindro1(:,2)=b*sin(theta);
                
                [I,dist] = knnsearch(cilindro1,Xi);
                
                m=MFE(Xi,dist);
                if m<mfe
                    mfe=m;
                    aFin=a;
                    bFin=b;
                end
                coeffImpl=[aFin bFin];
            end
        end        
    else
        coeffImpl=[];
        mfe=NaN;
    end
end