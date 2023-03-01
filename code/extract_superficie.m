function  [coord, maxCoord, H]= extract_superficie(curveId, xy, PtMin, PtMax, Eps,D,n,s)
a=PtMin(1):Eps(1):PtMax(1);
Na=numel(a);
if numel(Eps)>1
    b=PtMin(2):Eps(2):PtMax(2);
    Nb=numel(b);
    H=zeros(Na,Nb);
    [A,B]=HT(curveId,xy,b,s);
    
    for i=1:size(A,1)
        for j=1:Na
            for k=1:Nb
                aus=1;
                for t=1:size(A,2)
                    if (a(j)-Eps(1)/2<=A(i,t) && A(i,t)<(a(j)+Eps(1)/2)) && (b(k)-Eps(2)/2<=B(i,t) && B(i,t)<(b(k)+Eps(2)/2) && aus)
                        H(j,k)=H(j,k)+1;
                        aus=0;
                    end
                end
            end
        end
    end
    [row,col]=find(H == max(max(H)));
    aCoord=a(row);
    bCoord=b(col);
    coord=[aCoord, bCoord];
    maxCoord=max(max(H));

else
    H=zeros(Na,1);
    b=0;
    [A,B]=HT(curveId,xy,b,s);
    for i=1:size(A,1)
        for j=1:Na
            aus=1;
            for t=1:size(A,2)
                if (a(j)<=A(i,t) && A(i,t)<(a(j)+Eps(1)) && aus)
                    H(j)=H(j)+1;
                    aus=0;
                end
            end
        end
    end
    [row]=find(H == max(H));
    coord=a(row);
    maxCoord=max(H);
 
end