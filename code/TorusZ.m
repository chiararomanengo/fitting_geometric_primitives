function [coeffImpl, mfe,Max]=TorusZ(A,B,xyzAus)

    D=0;
    n=0;
    t=0;

    PtMin=[0.99*A 0.99*B];
    PtMax=[1.01*A 1.01*B];

    Eps = (1/10)*[(PtMax(1)-PtMin(1)) (PtMax(2)-PtMin(2))];

    [Coord, Max,H]=extract_superficie('TorusZ', xyzAus, PtMin, PtMax, Eps,D,n,t);

    if Max>0 
        mfe=10;
        for s=1:numel(Coord)/2
            a=Coord(s);
            b=Coord(numel(Coord)/2+s);
            
            theta=0:pi/700:2*pi;
            curva=zeros(numel(theta),3);
            curva(:,1)=cos(theta);
            curva(:,2)=sin(theta);
            curva(:,3)=zeros(numel(theta),1);
            t=0:2*pi/(numel(theta)-1):2*pi;
            
            torus=zeros(numel(theta)*numel(t),3);
            for j=1:numel(theta)
                for k=1:numel(t)
                    torus(k+(j-1)*numel(t),1)=curva(j,1)*(a+b*cos(t(k)));
                    torus(k+(j-1)*numel(t),3)=curva(j,3)+b*sin(t(k));
                    torus(k+(j-1)*numel(t),2)=curva(j,2)*(a+b*cos(t(k)));
                end
            end


            [I,dist] = knnsearch(torus,xyzAus);

            m=MFE(xyzAus,dist);
            if m<mfe
                mfe=m;
                aFin=a;
                bFin=b;
            end
            coeffImpl=[aFin bFin];
        end
    else
        coeffImpl=[];
        mfe=NaN;
    end

end