function bwGT=wsi_read_GT(fullFileName,rn,cn,tt)
s=xml2struct(fullFileName);
rr=s.Annotations.Annotation.Regions.Region;
bwGT=false(rn,cn);
if isstruct(rr)
    bw=false(rn,cn);
    v=rr.Vertices;
    m=cell2mat(v.Vertex);
    temp=[m.Attributes];
    X={temp.X};
    X=str2double(X);
    
    Y={temp.Y};
    Y=str2double(Y);
    
    X=round(X/tt);
    Y=round(Y/tt);
    
    bw=roipoly(bw,X,Y);
    bwGT=bwGT|bw;
else
    for n=1:length(rr)
        bw=false(rn,cn);
        v=rr{1,n}.Vertices;
        m=cell2mat(v.Vertex);
        temp=[m.Attributes];
        X={temp.X};
        X=str2double(X);
        
        Y={temp.Y};
        Y=str2double(Y);
        
        X=round(X/tt);
        Y=round(Y/tt);
        
        bw=roipoly(bw,X,Y);
        bwGT=bwGT|bw;
    end
end
