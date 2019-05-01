% R_L, R_U, k1, k2, n, A

param_min = [50; 200; 0.5; 2; 0.75; 0.05];
param_max = [250; 400; 3; 10; 0.95; 0.4];

% sample parameters
n = 500; 
P = param_min + (param_max-param_min).*lhsdesign(n,length(param_min),'criterion','correlation')';

% sampled monthly rainfall
out = [1  213.3726546   17.58811 336.3574
2  140.8576131   17.91813 284.1882
3  156.7036108   73.10406 295.5932
4  249.6492613  123.14563 361.6807
5  301.4888766  183.92816 426.9220
6  351.0909768  261.30948 495.4988
7  360.6242472  203.05669 467.1044
8  238.7224165   85.65108 335.6017
9   77.5569640  -78.55215 166.5825
10  -0.7517445 -150.36941 115.8439
11  52.7450363  -88.79477 146.0304
12 189.8713339   54.89617 312.6493];

minR = out(:,3);
minR(minR<0) = 0;
maxR = out(:,4);

m = 12;
R = zeros(m,n);
for t = 1:m
    
    x = mod(t,12);
    if x == 0
        x = 12;
    end
    
    R(t,:) = minR(x)+(maxR(x) - minR(x)).*lhsdesign(n,1,'criterion','correlation')';
    
end
 

save('MBRFunParams.mat','P','R');