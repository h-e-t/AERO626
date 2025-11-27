function S = skewSym(v)
    a = size(v,1)==3;
    n = size(v,round(1+a));
    V = permute(v,[round(2-a),3,round(1+a)]);
    I = repmat(eye(3),[1,1,n]);
    S = cross(repmat(V,[1,3]),I);
end