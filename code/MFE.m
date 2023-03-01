function mfe=MFE(xyz,dist)

    base=max(xyz(:,1))-min(xyz(:,1));
    h=max(xyz(:,2))-min(xyz(:,2));
    diag=sqrt(base^2+h^2);
    h=max(xyz(:,3))-min(xyz(:,3));
    Fin=sqrt(diag^2+h^2);
    
    mfe=mean(dist)/Fin;

end