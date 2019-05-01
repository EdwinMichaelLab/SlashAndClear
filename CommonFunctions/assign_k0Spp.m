function k0Spp=assign_k0Spp(bCulex,k0_Cx,k0_An)
k0Spp=zeros(size(bCulex));
if length(bCulex) == 1 && bCulex == 1
    k0Spp(1)=k0_Cx;
elseif length(bCulex) == 1 && bCulex == 0
    k0Spp(1)=k0_An;
elseif length(bCulex) > 1
    for jj=1:length(bCulex)
        if bCulex(jj)==1
            k0Spp(jj)=k0_Cx;
        else
            k0Spp(jj)=k0_An;
        end
    end
end
end