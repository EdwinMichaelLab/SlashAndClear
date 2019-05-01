function mat=wormMatingProb(w,k)
mat = 1 - (1 + w./(2*k)).^(-(1+k));
end