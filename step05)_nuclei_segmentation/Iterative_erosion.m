function t_e=Iterative_erosion(bw,SE)
t_e=bw;
thr=0.33;
while (sum(t_e(:))/sum(bw(:))>thr)
    t_e=imerode(t_e,SE);
end

