function [coord, maxCoord]=PlanesHT0(X)

rho= max(sqrt(X(:,1).^2+X(:,2).^2+X(:,3).^2));
a=0:rho/50:rho;
b=0:3:357;
b=(pi*b)/180;
c=0:2:180;
c=(pi*c)/180;

soglia=rho;
soglia(find(soglia==0))=10^(-2);

Na=numel(a);
Nb=numel(b);
Nc=numel(c);

H=zeros(Na,Nb,Nc);

for j=1:Na
    for k=1:Nb
        for h=1:Nc
            H(j,k,h)=numel(find(abs(a(j)-X(:,1)*cos(b(k))*sin(c(h))-X(:,2)*sin(c(h))*sin(b(k))-X(:,3)*cos(c(h)))<soglia/200));
        end
    end
end

maxCoord=max(max(max(H)));

coord=[];
for j=1:Na
    for k=1:Nb
        for h=1:Nc
            if  H(j,k,h)==maxCoord
                coord=[coord; a(j) b(k) c(h)];
            end            
        end
    end
end


end