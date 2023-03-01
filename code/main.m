clc; clear; close all;

%% Select a point cloud

t=4;

nameFile = strcat('../test/pointCloud/pointCloud', int2str(t) , '.txt');

xyz=load(nameFile);

figure
hold on
scatter3(xyz(:,1),xyz(:,2),xyz(:,3),'.k');
axis equal

%%

mfe=ones(5,1);
output=cell(5,1);

patchType=cell(1,1);
outputFin=cell(1,1);
mfeFin=ones(1,1);
 
i=1;
[output{1,i}, mfe(1,i)]=SearchPlane(xyz);
if mfe(1,i)>10^(-3) || isnan(mfe(1,i))
    try
        [output{2,i},mfe(2,i)]=SearchCylinder(xyz);
    catch
        coeff{2,i}=[];
        coeffParam{2,i}=[];
        mfe(2,i)=NaN;
    end
    if mfe(2,i)>10^(-3) || isnan(mfe(2,i))
        [output{3,i},mfe(3,i)]=SearchSphere(xyz);
        if mfe(3,i)>10^(-3) || isnan(mfe(3,i))
            try
                [output{4,i}, mfe(4,i)]=SearchCone(xyz);
            catch

                mfe(4,i)=NaN;
            end
            if mfe(4,i)>10^(-3) || isnan(mfe(4,i))
                try
                    [output{5,i},mfe(5,i)]=SearchTorus(xyz);
                catch
                    mfe(5,i)=NaN;
                    output{5,i}=[];
                end
            end
        end
    end
end

[M,I] = min(mfe);

if I(i)==5 && abs((mfe(5,i)-mfe(3,i)))<5*10^(-3)
    I(i)=3;
end

if I(i)==2 && abs((mfe(4,i)-mfe(2,i)))<5*10^(-3)
    I(i)=4;
end

switch I(i)
    case 1
        patchType{1,i}='Plane';
        outputFin{1,i}=[1; output{1,i}'];
        mfeFin(i)=mfe(1,i);
    case 2
        patchType{1,i}='Cylinder';
        outputFin{1,i}=[2; output{2,i}'];
        mfeFin(i)=mfe(2,i);
    case 3
        patchType{1,i}='Sphere';
        outputFin{1,i}=[3; output{3,i}'];
        mfeFin(i)=mfe(3,i);
    case 4
        patchType{1,i}='Cone';
        aus=output{4,i};
        aus(1)=[];
        outputFin{1,i}=[4; aus'];
        mfeFin(i)=mfe(4,i);
    case 5
        patchType{1,i}='Torus';
        outputFin{1,i}=[5; output{5,i}'];
        mfeFin(i)=mfe(5,i);
end


%%

nameFile = strcat('./testHT/pointCloud', int2str(t) , '_prediction.txt');

descriptors=outputFin{1,i};
writematrix(descriptors,nameFile);

