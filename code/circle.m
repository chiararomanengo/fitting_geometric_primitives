function [A,B]=ellisse(xyz)

    t=0:(2*pi)/1000:2*pi;
    A=zeros(size(xyz,1),numel(t));
    B=zeros(size(xyz,1),numel(t));
    for i=1:size(xyz,1)
        for j=1:numel(t)
              if (abs(sin(t(j)))>0.00000000001 && abs(cos(t(j)))>0.00000000001)
                    mat=[cos(t(j)) 0; 0 sin(t(j))];
                    sol=pinv(mat)*xyz(i,:)';
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