function [output, mfe] = SearchTorusNoise(X)

trX=mean(X(:,1));
trY=mean(X(:,2));
trZ=mean(X(:,3));

X(:,1)=X(:,1)-trX;
X(:,2)=X(:,2)-trY;
X(:,3)=X(:,3)-trZ;

maxValue=max(vecnorm(X'));
coord=PlanesHT0_search(X);
indAus=[];
mAus=1;
for j=1:size(coord,1)
    rho=coord(j,1);
    theta=coord(j,2);
    phi=coord(j,3);
    if (j>1 && (abs(phi-coord(j-1,3))>10^(-3) || abs(phi-coord(j-1,3)-pi)>10^(-3))) || j==1
        if abs(cos(phi))>10^(-3)
            x_piano = linspace(min(X(:,1)),max(X(:,1)),100);
            y_piano=linspace(min(X(:,2)),max(X(:,2)),100);
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
                y_piano=linspace(min(X(:,2)),max(X(:,2)),100);
                z_piano=linspace(min(X(:,3)),max(X(:,3)),100);
                piano = zeros(length(y_piano)*length(z_piano),3);
                for i=1:length(y_piano)
                    for k=1:length(z_piano)
                        piano( (i-1)*length(y_piano)+k,1) = (rho-y_piano(i)*sin(theta)*sin(phi))/(sin(phi)*cos(theta));
                        piano( (i-1)*length(y_piano)+k,2) = y_piano(i);
                        piano( (i-1)*length(y_piano)+k,3) = z_piano(k);
                    end
                end
            else
                x_piano=linspace(min(X(:,1)),max(X(:,1)),100);
                z_piano=linspace(min(X(:,3)),max(X(:,3)),100);
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
        
        [I,dist] = knnsearch(piano,X);
        cost=0.03;
        ind=find(dist<=cost*maxValue);
        while numel(ind)==0
            cost=cost+0.01;
            ind=find(dist<cost*maxValue);
        end
        
        [I,dist] = knnsearch(piano,X(ind,:));
        m=MFE(X(ind,:),dist);
        if m<mAus && numel(ind)>numel(indAus)           
            indAus=ind;
            mAus=m;
            vPlaneTot=[cos(theta)*sin(phi) sin(phi)*sin(theta) cos(phi)];
            XAus=X(ind,:);
        end
    end   
end

[U,~,~] = pca(XAus);
n_bf = cross(U(:,1),U(:,2));

R = rot_mat([0;0;1],n_bf);
XAus = XAus*R;
zAus=mean(XAus(:,3));   
xy=XAus(:,1:2);

try
    mBB = minBoundingBox(xy');
catch
    ind=[];
    mfe=NaN;
    output=[];
    return
end

if (norm(mBB(:,1)-mBB(:,2)) < norm(mBB(:,1)-mBB(:,4)) )
    m = (mBB(2,4)-mBB(2,1))/(mBB(1,4)-mBB(1,1));
    alpha = atan(m);
else
    m = (mBB(2,2)-mBB(2,1))/(mBB(1,2)-mBB(1,1));
    alpha = atan(m);
end
R1 = [cos(-alpha) sin(-alpha); -sin(-alpha) cos(-alpha)];
xy=xy*R1;

Aus=xy(xy(:,1)==max(xy(:,1)),:);
Aus2=xy(xy(:,1)==min(xy(:,1)),:);
TrY1=(Aus(1,2)+Aus2(1,2))/2;
TrX1=(max(xy(:,1))+min(xy(:,1)))/2;
xy(:,1)=xy(:,1)-TrX1;
xy(:,2)=xy(:,2)-TrY1;

if (abs(min(xy(:,2)))<=0.7*max(xy(:,2)) || abs(max(xy(:,2)))<=0.7*abs(min(xy(:,2)))) % caso mezzo toro
    corda=max(xy(:,1))-min(xy(:,1));
    f=(max(xy(:,2))-min(xy(:,2)));
    if f < 0.1*corda
        r=1.2*((corda^2/(8*f))+f/2);
    else
        r=((corda^2/(8*f))+f/2);
    end
    if abs(min(xy(:,2)))<=0.7*max(xy(:,2))
        TrY11=r-max(xy(:,2));
        xy(:,2)=xy(:,2)+TrY11;
        TrY11=-TrY11;
    else
        TrY11=r-abs(min(xy(:,2)));
        xy(:,2)=xy(:,2)-TrY11;
    end
    TrX11=0;
    A=r;
else
    if abs(max(xy(:,1))+min(xy(:,1)))/2<0.1*max(xy(:,1)) && abs(max(xy(:,2))+min(xy(:,2)))/2<0.1*max(xy(:,1)) % caso toro inter0
        TrX11=(max(xy(:,1))+min(xy(:,1)))/2;
        TrY11=(max(xy(:,2))+min(xy(:,2)))/2;
        xy(:,2)=xy(:,2)-TrY11;
        xy(:,1)=xy(:,1)-TrX11;
        A=max(xy(:,1));
    else
        if (numel(find((xy(:,1).^2+xy(:,2).^2)<(max(xy(:,1))/10)^2))==0) 
            TrX11=0;
            TrY11=0;
            A=max(xy(:,1));
        else 
            aus=xy(xy(:,2)==max(xy(:,2)),:);
            corda=max(xy(:,1))-min(xy(:,1));
            f=(max(xy(:,2))-min(xy(:,2)))/2;
            r=(corda^2/(8*f))+f/2;
            if abs(aus(1))<0.01
                TrY11=r+max(xy(:,2));
                xy(:,2)=xy(:,2)+TrY11;
                TrY11=-TrY11;
                A=max(xy(:,2));
            else
                TrY11=r-abs(min(xy(:,2)));
                xy(:,2)=xy(:,2)-TrY11;
                TrX11=0;
                A=min(abs(xy(:,2)));
            end
            AUS=0;
        end
    end
end

aus=max(max(xy(:,1)),max(xy(:,2)));
[xb, yb, rb] = find_circle(xy, -aus, aus, -aus, aus, 0.5*A,1.5*A);
xyz=[xy zeros(size(xy,1),1)];
a=rb;
theta=0:pi/200:2*pi;
circ=ones(numel(theta),3);
circ(:,3)=circ(:,3)*zAus;
circ(:,1)=xb+a*cos(theta);
circ(:,2)=yb+a*sin(theta); 
[I,dist] = knnsearch(circ,xyz);
mAus=MFE(xyz,dist);

points=circ;
TrY=TrY1+TrY11;
TrX=TrX1+TrX11;
points(:,1)=points(:,1)+TrX;
points(:,2)=points(:,2)+TrY;
points(:,1:2)=points(:,1:2)*R1';
points=points*R';

XAus=points;

pc=pointCloud(X);
step=0.3;
PP=pcdownsample(pc,'gridAverage',maxValue*step);
PP=PP.Location;

step=step-0.1;
while size(PP,1)/size(X,1)<0.05 && step>=0.1
    PP=pcdownsample(pc,'gridAverage',maxValue*step);
    PP=PP.Location;
    step=step-0.1;
end

ind1=[];
v=[];
b=[];
P=[];
Pn=[];
for t=1:size(PP,1)
    [I,dist] = knnsearch(PP(t,:),X);
    cost=0.05;
    ind=find(dist<=cost*maxValue);
    while numel(ind)<10
        cost=cost+0.02;
        ind=find(dist<cost*maxValue);
    end
    XXX=X(ind,:);

    [coord, maxCoord]=PlanesHT0(XXX);
    
    for j=1:size(coord,1)
        rho=coord(j,1);
        theta=coord(j,2);
        phi=coord(j,3);
        if abs(cos(phi))>10^(-3)
            x_piano = linspace(min(X(:,1)),max(X(:,1)),100);
            y_piano=linspace(min(X(:,2)),max(X(:,2)),100);
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
                y_piano=linspace(min(X(:,2)),max(X(:,2)),100);
                z_piano=linspace(min(X(:,3)),max(X(:,3)),100);
                piano = zeros(length(y_piano)*length(z_piano),3);
                for i=1:length(y_piano)
                    for k=1:length(z_piano)
                        piano( (i-1)*length(y_piano)+k,1) = (rho-y_piano(i)*sin(theta)*sin(phi))/(sin(phi)*cos(theta));
                        piano( (i-1)*length(y_piano)+k,2) = y_piano(i);
                        piano( (i-1)*length(y_piano)+k,3) = z_piano(k);
                    end
                end
            else
                x_piano=linspace(min(X(:,1)),max(X(:,1)),100);
                z_piano=linspace(min(X(:,3)),max(X(:,3)),100);
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
        
        [I,dist] = knnsearch(piano,X);
        cost=0.01;
        ind=find(dist<=cost*maxValue);
        while numel(ind)==0
            cost=cost+0.01;
            ind=find(dist<cost*maxValue);
        end
        
        ind1=[ind1; ind];
        ind1=unique(ind1);
        v=[v; cos(theta)*sin(phi) sin(phi)*sin(theta) cos(phi)];
        b=[b; rho];
        P=[P; PP(t,:)];
        indexAus=randi(size(piano,1));
        Pn=[Pn; piano(indexAus,:)];  
    end
end

vFin=[];
PFin=[];
Ppiano=[];
for s=1:size(v,1)
    aus=zeros(size(v,1),1);
    for j=1:size(v,1)
        if P(j,1)==P(s,1)
            aus(j)=1;
        end
    end
    if numel(find(aus))>1
        vAus=mean(v(find(aus),:));
        contr=1;
        for k=1:size(vFin,1)
            if vFin(k,1)==vAus(1)
                contr=0;
            end
        end
        if contr
            vFin=[vFin; vAus];
            PFin=[PFin; P(s,:)];
            Ppiano=[Ppiano; Pn(s,:)];
        end
    else
        vAus=v(find(aus),:);
        contr=1;
        for k=1:size(vFin,1)
            if vFin(k,1)==vAus(1)
                contr=0;
            end
        end
        if contr
            vFin=[vFin; vAus];
            PFin=[PFin; P(s,:)];
            Ppiano=[Ppiano; Pn(s,:)];
        end
    end
end

vVert=[];
coeffV=[];
for i=1:size(vFin,1)
    n=vFin(i,:);
    [I,dist] = knnsearch(XAus,PFin(i,:));
    point=XAus(I,:);
    t=point-PFin(i,:);
    if norm(t)~=0 && norm(n)~=0
        v1=cross(n,t);
        vVert=[vVert; v1];
        coeffV=[coeffV; -v1(1)*PFin(i,1)-v1(2)*PFin(i,2)-v1(3)*PFin(i,3)];
    end
end

C=[];
soglia=0.005;
GoodPoint=[];
GoodNorm=[];
radius1=[];
while size(C,1)<10
    for t=1:size(vVert,1)
        rho=coeffV(t);
        a=vVert(t,1);
        b=vVert(t,2);
        c=vVert(t,3);  
        if norm([a b c])>10^(-6)
            if abs(c)>10^(-3)
                x_piano = linspace(min(X(:,1)),max(X(:,1)),100);
                y_piano=linspace(min(X(:,2)),max(X(:,2)),100);
                piano = zeros(length(x_piano)*length(y_piano),3);
                for i=1:length(x_piano)
                    for k=1:length(y_piano)
                        piano( (i-1)*length(x_piano)+k,1) = x_piano(i);
                        piano( (i-1)*length(x_piano)+k,2) = y_piano(k);
                        piano( (i-1)*length(x_piano)+k,3) = (-rho-x_piano(i)*a-y_piano(k)*b)/c;
                    end
                end
            else
                if abs(a)>10^(-3)
                    y_piano=linspace(min(X(:,2)),max(X(:,2)),100);
                    z_piano=linspace(min(X(:,3)),max(X(:,3)),100);
                    piano = zeros(length(y_piano)*length(z_piano),3);
                    for i=1:length(y_piano)
                        for k=1:length(z_piano)
                            piano( (i-1)*length(y_piano)+k,1) = (-rho-y_piano(i)*b)/(a);
                            piano( (i-1)*length(y_piano)+k,2) = y_piano(i);
                            piano( (i-1)*length(y_piano)+k,3) = z_piano(k);
                        end
                    end
                else
                    x_piano=linspace(min(X(:,1)),max(X(:,1)),100);
                    z_piano=linspace(min(X(:,3)),max(X(:,3)),100);
                    piano = zeros(length(x_piano)*length(z_piano),3);
                    for i=1:length(x_piano)
                        for k=1:length(z_piano)
                            piano( (i-1)*length(x_piano)+k,1) = x_piano(i);
                            piano( (i-1)*length(x_piano)+k,2) = (-rho)/(b);
                            piano( (i-1)*length(x_piano)+k,3) = z_piano(k);
                        end
                    end
                end
            end
            
            [I,dist] = knnsearch(piano,X);
            cost=0.02;
            ind=find(dist<=cost*maxValue);
            while numel(ind)<15
                cost=cost+0.02;
                ind=find(dist<cost*maxValue);
            end
            XXX=X(ind,:);
            
            parN=round((size(XXX,1)/10)*3);
            minPts=round((size(XXX,1)/10));
            
            [IDX2,DDD] = knnsearch(XXX,XXX,'k', parN);
            epsilon = mean(DDD(:, parN));
            IDX2=DBSCAN(XXX, epsilon, minPts);
            numCl=max(IDX2);
            
            if numCl<3
                for i=1:numCl-1
                    ind1=find(IDX2==i);
                    for j=1:numCl
                        ind2=find(IDX2==j);
                        if numel(ind1)>numel(ind2)
                            indMax=ind1;
                        else
                            indMax=ind2;
                        end
                    end
                end
                
                if numCl==1
                    indMax=find(IDX2==1);
                end
                
                xyz=XXX(indMax,:);
                
                vPlane=vVert(t,:)/norm(vVert(t,:));
                
                [U,~,~] = pca(xyz);
                n_bf = cross(U(:,1),U(:,2));
                
                R = rot_mat([0;0;1],n_bf);
                xyz = xyz*R;
                
                zAus=mean(xyz(:,3));
                xy=xyz(:,1:2);
                try
                    mBB = minBoundingBox(xy');
                catch
                    ind=[];
                    mfe=NaN;
                    output=[];
                end
                
                if (norm(mBB(:,1)-mBB(:,2)) < norm(mBB(:,1)-mBB(:,4)) )
                    m = (mBB(2,4)-mBB(2,1))/(mBB(1,4)-mBB(1,1));
                    alpha = atan(m);
                else
                    m = (mBB(2,2)-mBB(2,1))/(mBB(1,2)-mBB(1,1));
                    alpha = atan(m);
                end
                R1 = [cos(-alpha) sin(-alpha); -sin(-alpha) cos(-alpha)];
                xy=xy*R1;
                
                Aus=xy(xy(:,1)==max(xy(:,1)),:);
                Aus2=xy(xy(:,1)==min(xy(:,1)),:);
                TrY1=(Aus(1,2)+Aus2(1,2))/2;
                TrX1=(max(xy(:,1))+min(xy(:,1)))/2;
                xy(:,1)=xy(:,1)-TrX1;
                xy(:,2)=xy(:,2)-TrY1;
                if (abs(min(xy(:,2)))<=0.5*max(xy(:,2)) || abs(max(xy(:,2)))<=0.5*abs(min(xy(:,2)))) % caso mezzo toro
                    corda=max(xy(:,1))-min(xy(:,1));
                    f=(max(xy(:,2))-min(xy(:,2)));
                    if f < 0.1*corda
                        r=1.2*((corda^2/(8*f))+f/2);
                    else
                        r=((corda^2/(8*f))+f/2);
                    end
                    if abs(min(xy(:,2)))<=0.5*max(xy(:,2))
                        TrY11=r-max(xy(:,2));
                        xy(:,2)=xy(:,2)+TrY11;
                        TrY11=-TrY11;
                    else
                        TrY11=r-abs(min(xy(:,2)));
                        xy(:,2)=xy(:,2)-TrY11;
                    end
                    TrX11=0;
                    A=r;
                else
                    if abs(max(xy(:,1))+min(xy(:,1)))/2<0.1*max(xy(:,1)) && abs(max(xy(:,2))+min(xy(:,2)))/2<0.1*max(xy(:,1)) % caso toro inter0
                        TrX11=(max(xy(:,1))+min(xy(:,1)))/2;
                        TrY11=(max(xy(:,2))+min(xy(:,2)))/2;
                        xy(:,2)=xy(:,2)-TrY11;
                        xy(:,1)=xy(:,1)-TrX11;
                        A=max(xy(:,1));
                    else
                        if (numel(find((xy(:,1).^2+xy(:,2).^2)<(max(xy(:,1))/10)^2))==0) 
                            TrX11=0;
                            TrY11=0;                          
                            A=max(xy(:,1));
                        else
                            aus=xy(xy(:,2)==max(xy(:,2)),:);
                            corda=max(xy(:,1))-min(xy(:,1));
                            f=(max(xy(:,2))-min(xy(:,2)))/2;
                            r=(corda^2/(8*f))+f/2;
                            if abs(aus(1))<0.01
                                TrY11=r+max(xy(:,2));
                                xy(:,2)=xy(:,2)+TrY11;
                                TrY11=-TrY11;
                                A=max(xy(:,2));
                            else
                                TrY11=r-abs(min(xy(:,2)));
                                xy(:,2)=xy(:,2)-TrY11;
                                TrX11=0;
                                A=min(abs(xy(:,2)));
                            end
                            AUS=0;
                        end
                    end
                end
                
                aus=max(max(xy(:,1)),max(xy(:,2)));
                [xb, yb, rb] = find_circle(xy, -aus, aus, -aus, aus, 0.5*A,1.5*A);
                xyz=[xy zeros(size(xy,1),1)];
                a=rb;
                theta=0:pi/400:2*pi;
                circ=zeros(numel(theta),3);
                circ(:,1)=xb+a*cos(theta);
                circ(:,2)=yb+a*sin(theta); 
                [I,dist] = knnsearch(circ,xyz);
                mAus=MFE(xyz,dist);
                
                if mAus<soglia 
                    TrY=TrY1+TrY11+yb;
                    TrX=TrX1+TrX11+xb;
                    P=[TrX TrY zAus];
                    P(1:2)=P(1:2)*R1';
                    P=P*R';
                    C=[C; P];
                    GoodPoint=[GoodPoint; PFin(t,:)];
                    GoodNorm=[GoodNorm; vFin(t,:)];
                    radius1=[radius1; rb];
                end
            end
        end
    end
    soglia=soglia+0.005;
end

parN=round((size(C,1)/10)*3);
minPts=round((size(C,1)/5)*1);

[IDX2,DDD] = knnsearch(C,C,'k', parN);
epsilon = mean(DDD(:, parN));
IDX2=DBSCAN(C, epsilon, minPts);

if numel(find(IDX2==1))>5
    aus=C(find(IDX2==1),:);
    sumC=0;
    for s=1:size(aus,1)-1
        for k=s+1:size(aus,1)
            sumC=sumC+sum(aus(s,:)-aus(k,:));
        end
    end
    if sumC>0
        C=C(find(IDX2==1),:);
        radius1=radius1(find(IDX2==1),:);
    end
end

[coord, maxCoord]=PlanesHT0_search_center(C);


mAus=1;
for j=1:size(coord,1)
    rho=coord(j,1);
    theta=coord(j,2);
    phi=coord(j,3);
    if (j>1 && (abs(phi-coord(j-1,3))>10^(-3) || abs(phi-coord(j-1,3)-pi)>10^(-3))) || j==1
        if abs(cos(phi))>10^(-3)
            x_piano = linspace(min(X(:,1)),max(X(:,1)),400);
            y_piano=linspace(min(X(:,2)),max(X(:,2)),400);
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
                y_piano=linspace(min(X(:,2)),max(X(:,2)),400);
                z_piano=linspace(min(X(:,3)),max(X(:,3)),400);
                piano = zeros(length(y_piano)*length(z_piano),3);
                for i=1:length(y_piano)
                    for k=1:length(z_piano)
                        piano( (i-1)*length(y_piano)+k,1) = (rho-y_piano(i)*sin(theta)*sin(phi))/(sin(phi)*cos(theta));
                        piano( (i-1)*length(y_piano)+k,2) = y_piano(i);
                        piano( (i-1)*length(y_piano)+k,3) = z_piano(k);
                    end
                end
            else
                x_piano=linspace(min(X(:,1)),max(X(:,1)),400);
                z_piano=linspace(min(X(:,3)),max(X(:,3)),400);
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
        [I,dist] = knnsearch(piano,C);
        m=MFE(C,dist);
        if m<mAus
            mAus=m;
            ax=[cos(theta)*sin(phi) sin(phi)*sin(theta) cos(phi)];
        end
    end
end

if (norm(ax-vPlaneTot)<1.1 || norm(ax+vPlaneTot)<1.1)
    ax=vPlaneTot;
end

vPlane=ax;

z=[0 0 1];
if norm(z-vPlane)>2*10^(-1) && norm(z+vPlane)>2*10^(-1)
    R = rot_mat(z',vPlane');
    xyz=C*R;
    n=ax;
else
    R = eye(3);
    xyz=C*R;
    n=z;
end

zAus=mean(xyz(:,3));
xy=xyz(:,1:2);
try
    mBB = minBoundingBox(xy');
catch
    ind=[];
    mfe=NaN;
    output=[];
    return
end

if (norm(mBB(:,1)-mBB(:,2)) < norm(mBB(:,1)-mBB(:,4)) )
    m = (mBB(2,4)-mBB(2,1))/(mBB(1,4)-mBB(1,1));
    alpha = atan(m);
else
    m = (mBB(2,2)-mBB(2,1))/(mBB(1,2)-mBB(1,1));
    alpha = atan(m);
end
R1 = [cos(-alpha) sin(-alpha); -sin(-alpha) cos(-alpha)];
xy=xy*R1;

Aus=xy(xy(:,1)==max(xy(:,1)),:);
Aus2=xy(xy(:,1)==min(xy(:,1)),:);
TrY1=(Aus(1,2)+Aus2(1,2))/2;
TrX1=(max(xy(:,1))+min(xy(:,1)))/2;
xy(:,1)=xy(:,1)-TrX1;
xy(:,2)=xy(:,2)-TrY1;

if (abs(min(xy(:,2)))<=0.7*max(xy(:,2)) || abs(max(xy(:,2)))<=0.7*abs(min(xy(:,2)))) 
    corda=(max(xy(:,1))-min(xy(:,1)));
    f=(max(xy(:,2))-min(xy(:,2)));
    if f < 0.3*corda 
        r1=mean(radius1);
        r=1.2*((corda^2/(8*f))+f/2);
        if r/r1>10
            r=corda;
        end
    else
        r=((corda^2/(8*f))+f/2);
    end
    if abs(min(xy(:,2)))<=0.7*max(xy(:,2))
        TrY11=r-max(xy(:,2));
        xy(:,2)=xy(:,2)+TrY11;
        TrY11=-TrY11;
    else
        TrY11=r-abs(min(xy(:,2)));
        xy(:,2)=xy(:,2)-TrY11;
    end
    TrX11=0;
    A=r;
else
    if abs(max(xy(:,1))+min(xy(:,1)))/2<0.1*max(xy(:,1)) && abs(max(xy(:,2))+min(xy(:,2)))/2<0.1*max(xy(:,1)) % caso toro inter0
        TrX11=(max(xy(:,1))+min(xy(:,1)))/2;
        TrY11=(max(xy(:,2))+min(xy(:,2)))/2;
        xy(:,2)=xy(:,2)-TrY11;
        xy(:,1)=xy(:,1)-TrX11;
        A=max(xy(:,1));
    else
        if (numel(find((xy(:,1).^2+xy(:,2).^2)<(max(xy(:,1))/10)^2))==0) 
            TrX11=0;
            TrY11=0;
            A=max(xy(:,1));
        else 
            aus=xy(xy(:,2)==max(xy(:,2)),:);
            corda=max(xy(:,1))-min(xy(:,1));
            f=(max(xy(:,2))-min(xy(:,2)))/2;
            r=(corda^2/(8*f))+f/2;
            if abs(aus(1))<0.01
                TrY11=r+max(xy(:,2));
                xy(:,2)=xy(:,2)+TrY11;
                TrY11=-TrY11;
                A=max(xy(:,2));
            else
                TrY11=r-abs(min(xy(:,2)));
                xy(:,2)=xy(:,2)-TrY11;
                TrX11=0;
                A=min(abs(xy(:,2)));
            end
            AUS=0;
        end
    end
end

aus=max(max(xy(:,1)),max(xy(:,2)));
[xb, yb, rb,max_H] = find_circle(xy, -aus, aus, -aus, aus, 0.9*A,1.1*A);

mAus=1;
xyz=[xy zeros(size(xy,1),1)];
if max_H<5
    [coord, maxCoord, Ps]= extract_curve('circle', xy, 0.9*A, 1.1*A, A/50);
    if maxCoord>max_H
        xb=0;
        yb=0;
        for k=1:numel(coord)
            a1=coord(k);
            theta=0:pi/200:2*pi;
            circ=zeros(numel(theta),3);
            circ(:,1)=a1*cos(theta);
            circ(:,2)=a1*sin(theta); 
            [I,dist] = knnsearch(circ,xyz);
            m=MFE(xyz,dist);
            if m<mAus
                a=coord(k);
                mAus=m;
            end
        end
    else
        a=rb;
    end
else
    a=rb;
end

TrY=TrY1+TrY11+yb;
TrX=TrX1+TrX11+xb;
V=[TrX TrY zAus];
V(1:2)=V(1:2)/R1;
V=V/R;

r2=a;
r1=mean(radius1);

TrVx=V(1,1);
TrVy=V(1,2);
TrVz=V(1,3);

xyz=X;
xyz(:,1)=X(:,1)-TrVx;
xyz(:,2)=X(:,2)-TrVy;
xyz(:,3)=X(:,3)-TrVz;

z=[0 0 1];
if norm(z-n)>10^(-2) && norm(z+n)>10^(-2)
    R = rot_mat(z',n');
    xyz=xyz*R;
else
    R = eye(3);
    xyz=xyz*R;
end

[radius2, mfe2]=TorusZ(r2,r1,xyz);

V(1)=V(1)+trX;
V(2)=V(2)+trY;
V(3)=V(3)+trZ;

output=[radius2(1) radius2(2) n V];
mfe=mfe2;
end

