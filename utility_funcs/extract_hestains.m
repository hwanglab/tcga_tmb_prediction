% extract H,E vectors for image color normalization
% for paper reference, see?Macenko, Marc, et al. "A method for normalizing histology slides for quantitative analysis." 
% 2009 IEEE International Symposium on Biomedical Imaging: From Nano to Macro. IEEE, 2009.
% author: Hongming Xu, Ph.D. Cleveland Clinic
% feel free to use it, but should give original author credits

function HE=extract_hestains(RGB,bwTissue)

Io=240; 
beta=0.15;
alpha=1;

Irgb=reshape(RGB,[],3);
I=double(Irgb(bwTissue,:));
% calculate optical density
OD = -log((I+1)/Io);
% remove transparent pixels
ODhat = OD(~any(OD < beta, 2), :);
% calculate eigenvectors
[V, ~] = eig(cov(ODhat));

% project on the plane spanned by the eigenvectors corresponding to the two
% largest eigenvalues
That = ODhat*V(:,2:3);

% find the min and max vectors and project back to OD space
phi = atan2(That(:,2), That(:,1));

minPhi = prctile(phi, alpha);
maxPhi = prctile(phi, 100-alpha);

vMin = V(:,2:3)*[cos(minPhi); sin(minPhi)];
vMax = V(:,2:3)*[cos(maxPhi); sin(maxPhi)];

% a heuristic to make the vector corresponding to hematoxylin first and the
% one corresponding to eosin second
if vMin(1) > vMax(1)
    HE = [vMin vMax];
else
    HE = [vMax vMin];
end