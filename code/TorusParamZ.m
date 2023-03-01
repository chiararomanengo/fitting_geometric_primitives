function [A,B]= TorusParamZ(xyz)

    v=0:(2*pi)/50:2*pi;
    u=0:(2*pi)/(numel(v)-1):2*pi;
    A=zeros(size(xyz,1),numel(u)*numel(v));
    B=zeros(size(xyz,1),numel(u)*numel(v));

for j=1:length(u)
    for k=1:length(v)
        mat=[cos(v(k)) cos(v(k))*cos(u(j)); sin(v(k)) sin(v(k))*cos(u(j)); 0 sin(u(j))];
        mat=pinv(mat);
        for i=1:size(xyz,1)
            sol=mat*xyz(i,:)';
            A(i,(j-1)*length(u)+k) = sol(1);
            B(i,(j-1)*length(u)+k) =sol(2);
        end
    end
end

end