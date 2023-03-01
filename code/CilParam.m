function [A,B]= CilParam(xyz)

    t=(2*pi)/250:(2*pi)/250:2*pi;
    A=zeros(size(xyz,1),numel(t));
    B=zeros(size(xyz,1),numel(t));
    for j=1:numel(t)
         mat=[cos(t(j)) 0; 0 sin(t(j))];
         mat=pinv(mat);
        for i=1:size(xyz,1)
               if (abs(cos(t(j)))>0.00001 && abs(sin(t(j)))>0.00001)                  
                    sol=mat*xyz(i,:)';
                    A(i,j)=sol(1);
                    B(i,j)=sol(2);
               end
        end     
    end
    
    col=find(sum(A));
    A=A(:,col);
    B=B(:,col);
    col=find(sum(B));
    A=A(:,col);
    B=B(:,col);

end