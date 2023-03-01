function [coord, maxCoord]=searchVertex(v,rho,maxValue)

a=-maxValue:maxValue/150:maxValue;
b=-maxValue:maxValue/150:maxValue;
c=-maxValue:maxValue/150:maxValue;

Na=numel(a);
Nb=numel(b);
Nc=numel(c);

H=zeros(Na,Nb,Nc);

soglia=rho;
soglia(find(soglia==0))=10^(-2);

for j=1:Na
    for k=1:Nb
        for h=1:Nc
            H(j,k,h)=numel(find(abs(rho-v(:,1)*a(j)-v(:,2)*b(k)-v(:,3)*c(h))<soglia/20));
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