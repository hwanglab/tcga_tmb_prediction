%----usage-----
% compute multi-scale LBP features
% reference: Ojala, Timo, Matti Pietik‰inen, and Topi M‰enp‰‰. 
% "Multiresolution gray-scale and rotation invariant texture classification with local binary patterns." 
% IEEE Transactions on Pattern Analysis & Machine Intelligence 7 (2002): 971-987.

function FF=xu_LBP(Rimg,mapping)

rr=[2,4,6,8]; % this is typically multiscale LBP parameter setting
pp=mapping.samples;

FF=[];
for t=1:length(rr)
    L=LBP(Rimg,rr(t),pp,mapping);
    FF=[FF,L];
end