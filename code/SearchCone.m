function [output, mfe] = SearchCone(X)

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
        %     aus=1;
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

[V, maxCoord]=searchVertex(vAus,bAus,2*maxValue);

if size(V,1)>1
   V=mean(V); 
end

ax=[];
V=V';
for i=1:size(vAus,1)-1
    n=vAus(i,:);
    t=V'-PAus(i,:);
    v1=cross(n,t);   
    for j=i+1:size(vAus,1)
        n=vAus(j,:);
        t=V'-PAus(j,:);
        v2=cross(n,t);
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

v=ax;
indices=zeros(size(v,1),1);
for i=1:size(v,1)-1
    v1=v(i,:);
    for k=i+1:size(v,1)
        v2=v(k,:);
        if norm(v1-v2)<5*10^(-1)
           indices(i)= indices(i)+1;
           indices(k)= indices(k)+1;
        end
    end
    
end

ax=ax(indices==max(indices),:);

if size(ax,1)>1
    ax=mean(ax);
end
n=ax/norm(ax);

for i=1:numel(n)
   if abs(n(i)-1)<3*10^(-2) || abs(n(i)+1)<3*10^(-2)
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

TrVx=V(1);
TrVy=V(2);
TrVz=V(3);

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

rho=sqrt(xyz(:,1).^2+xyz(:,2).^2+xyz(:,3).^2);
[M,I]=max(rho);

aus=xyz(I,:);
A=sqrt(aus(1)^2+aus(2)^2);
if aus(3)>0 
    C=aus(3);
else
    C=abs(min(aus(3)));
end

theta=0:pi/50:2*pi;
t=-1:1/(numel(theta)-1):1;

v=[0 0 1];
u=[C/A 0 1];
A=acos(sum(v.*u)/(norm(v)*norm(u)));
A=pi/2-A;

[radius1, mfe1]=cone(A,xyz,t);

mAus=1;
if numel(radius1)>1
    for j=1:numel(radius1)
        alpha=radius1(j);
        theta=2*pi;
        N=80;
        h2=2*max(xyz(:,3));
        X = zeros(N^2,3);
        for i=1:size(X,1)
            u = rand(1)*theta;
            v = -h2+rand(1)*(h2*2);
            X(i,1)=sin(alpha)*cos(u)*v;
            X(i,2)=sin(alpha)*sin(u)*v;
            X(i,3)=cos(alpha)*v;
        end
        [I,dist] = knnsearch(X,xyz);
        
        m1=MFE(xyz,dist);
        if m1<mAus
            alphaOpt=alpha;
            mAus=m1;
        end
    end    
else
    alphaOpt=radius1;
end
mfe=min([mAus mfe1]);

V(1)=V(1)+trX;
V(2)=V(2)+trY;
V(3)=V(3)+trZ;

output=[radius1(1) alphaOpt n V'];

catch
    mfe=NaN;
    output=[];   
    return
end


end