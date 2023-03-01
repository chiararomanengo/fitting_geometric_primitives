function  [coord, maxCoord, Ps]= extract_curve(curveId, xy, PtMin, PtMax, Eps)
    a=PtMin(1):Eps(1):PtMax(1);
    Na=numel(a);
    if numel(Eps)>1
        b=PtMin(2):Eps(2):PtMax(2);
        Nb=numel(b);
        H=zeros(Na,Nb);
        [A,B]=HT(curveId,xy,a);
        for i=1:size(A,1)
            for j=1:Na
                for k=1:Nb
                    aus=1;
                    for t=1:size(A,2)
                        if (a(j)<=A(i,t) && A(i,t)<(a(j)+Eps(1))) && (b(k)<=B(i,t) && B(i,t)<(b(k)+Eps(2)) && aus)
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
   
         p=(60/100)*maxCoord;
        Ps=zeros(size(xy));
        for i=1:size(A,1)
            for j=1:Na
                for k=1:Nb
                    for t=1:size(A,2)
                        if (a(j)<=A(i,t) && A(i,t)<(a(j)+Eps(1))) && (b(k)<=B(i,t) && B(i,t)<(b(k)+Eps(2)))
                            if H(j,k)>p
                                Ps(i,:)=xy(i,:);
                            end
                        end
                    end   
                end     
            end  
        end
    else 
        H=zeros(Na,1);
        [A,B]=HT(curveId,xy,a);
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
   
        p=(60/100)*maxCoord;
        Ps=zeros(size(xy));
        for i=1:size(A,1)
            for j=1:Na
                for t=1:size(A,2)
                    if (a(j)<=A(i,t) && A(i,t)<(a(j)+Eps(1)))
                        if H(j)>p
                           Ps(i,:)=xy(i,:);
                        end
                    end
                end   
            end     
        end  
    end
   
    Ps=Ps(find(sum(Ps')),:);
   
end