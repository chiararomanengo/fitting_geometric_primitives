function [output, mfe]=SearchCylinder(X)

try
    
trX=mean(X(:,1));
trY=mean(X(:,2));
trZ=mean(X(:,3));

X(:,1)=X(:,1)-trX;
X(:,2)=X(:,2)-trY;
X(:,3)=X(:,3)-trZ;

maxValue=max(vecnorm(X'));
pc=pointCloud(X);
step=0.2;
PP=pcdownsample(pc,'gridAverage',maxValue*step);
PP=PP.Location;

ind1=[];
v=[];
b=[];
P=[];
Pn=[];
mfe=[];

for t=1:size(PP,1)
    cost=0.02;
    [ind,dist]=rangesearch(X,PP(t,:),cost*maxValue);
    ind=ind{1,1};
    while numel(ind)<10
        cost=cost+0.05;
        [ind,dist]=rangesearch(X,PP(t,:),cost*maxValue);
        ind=ind{1,1};
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
        
        [I,dist] = knnsearch(piano,X(ind,:));
        m=MFE(X(ind,:),dist);
        
        ind1=[ind1; ind];
        ind1=unique(ind1);
        v=[v; cos(theta)*sin(phi) sin(phi)*sin(theta) cos(phi)];
        b=[b; rho];
        P=[P; PP(t,:)];
        indexAus=randi(size(piano,1));
        Pn=[Pn; piano(indexAus,:)];
        mfe=[mfe; m];
    end
end

vFin=[];
PFin=[];
Ppiano=[];
bFin=[];
mfeFin=[];
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
            bFin=[bFin; b(s)];
            PFin=[PFin; P(s,:)];
            Ppiano=[Ppiano; Pn(s,:)];
            mfeFin=[mfeFin; mfe(s)];
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
            bFin=[bFin; b(s)];
            PFin=[PFin; P(s,:)];
            Ppiano=[Ppiano; Pn(s,:)];
            mfeFin=[mfeFin; mfe(s)];
        end
    end
end

vAus=vFin;
bAus=bFin;
PAus=PFin;

soglia=0.01;
ind=find(mfeFin>0.01);
vAus(ind,:)=[];
PAus(ind,:)=[];
bAus(ind)=[];

k=0;
while numel(bAus)<3
    soglia=soglia+0.001;
    ind=find(mfeFin>soglia);
    vAus=vFin;
    bAus=bFin;
    PAus=PFin;
    vAus(ind,:)=[];
    PAus(ind,:)=[];
    bAus(ind)=[];
    k=k+1; 
    if k>numel(mfeFin)
        output=[];
        mfe=NaN;
       return 
    end
end

ax=[];
for j=1:size(vAus,1)-1
    v1=vAus(j,:);
    for k=j+1:size(vAus,1)
        v2=vAus(k,:);
        if norm(v1-v2)>10^(-2) && norm(v1+v2)>10^(-2)
            coeff=cross(v1,v2); 
            for s=1:3
                if abs(coeff(s))<10^(-10)
                    coeff(s)=0;
                end
            end
            if coeff(1)<0
                coeff=-coeff;
            else
                if coeff(1)==0
                    if coeff(2)<0
                        coeff(2:3)=-coeff(2:3);
                    else
                        if coeff(2)==0
                            if coeff(3)<0
                                coeff(3)=-coeff(3);
                            end
                        end
                    end
                end
            end
            ax=[ax; coeff];
        end
    end
end

if size(ax,1)>1
    ax=mean(ax);
end
n=ax/norm(ax);

for i=1:numel(n)
    if abs(n(i)-1)<10^(-2) || abs(n(i)+1)<10^(-2)
        n(i)=1;
        if i==1
            n(2)=0;
            n(3)=0;
        else
            if i==2
                n(1)=0;
                n(3)=0;
            else
                n(1)=0;
                n(2)=0;
            end
        end
    end
end

M=mean(X);
X= bsxfun(@minus, X, M);
z=[0 0 1];
if norm(z-n)>10^(-2) && norm(z+n)>10^(-2)
    R = rot_mat(z',n');
    xyz=X*R;
else
    R = eye(3);
    xyz=X*R;
end

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


if (abs(min(xy(:,2)))<=0.1*max(xy(:,2)) || abs(max(xy(:,2)))<=0.1*abs(min(xy(:,2)))) 
    corda=max(xy(:,1))-min(xy(:,1));
    f=(max(xy(:,2))-min(xy(:,2)));
    r=(corda^2/(8*f))+f/2;
    if abs(min(xy(:,2)))<=0.1*max(xy(:,2))
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
[xb, yb, rb, max_H] = find_circle(xy, -aus/2, aus/2, -aus/2, aus/2, 0.8*A,1.2*A);

mAus=1;
xyzCirc=[xy zeros(size(xy,1),1)];
if max_H<2
    [coord, maxCoord, Ps]= extract_curve('circle', xy, 0.8*A, 1.2*A, A/100);
    if maxCoord>max_H
        xb=0;
        yb=0;
        for k=1:numel(coord)
            a1=coord(k);
            theta=0:pi/200:2*pi;
            circ=zeros(numel(theta),3);
            circ(:,1)=a1*cos(theta);
            circ(:,2)=a1*sin(theta);
            [I,dist] = knnsearch(circ,xyzCirc);
            m=MFE(xyzCirc,dist);
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


xy(:,1)=xy(:,1)-xb;
xy(:,2)=xy(:,2)-yb;
xyz(:,1:2)=xy;
[radius1, mfe1, max1]=cylinderZ(a,a,xyz);

TrY=TrY1+TrY11+yb;
TrX=TrX1+TrX11+xb;
P=[TrX TrY 0];
P(1:2)=P(1:2)/R1;
P=P/R;
P=P+M;
P(1)=P(1)+trX;
P(2)=P(2)+trY;
P(3)=P(3)+trZ;

a=radius1(1);
output=[a n P];
mfe=mfe1;

catch
    mfe=NaN;
    output=[];   
    return
end


end

