function [A,B]= ConeAlpha(xyz,v)

    u=0.1:(2*pi-0.1)/(numel(v)-1):2*pi;
    A=zeros(size(xyz,1),numel(u)*numel(v));

for j=1:length(u)
    for k=1:length(v)
        mat=[cos(u(j))*v(k); v(k)*sin(u(j)); v(k)];
        mat=pinv(mat);
        for i=1:size(xyz,1)
            sol=mat*xyz(i,:)';
            A(i,(j-1)*length(u)+k) = sol(1);
        end
    end
end

col=find(sum(A));
A=A(:,col);
B=0;

end