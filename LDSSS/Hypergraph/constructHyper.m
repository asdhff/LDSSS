function [H] = constructHyper(Z)
v = size(Z, 2);
H = cell(1,v);
for i=1:v
      Weight{i} = constructW_PKN((abs(Z{i})+abs(Z{i}'))./2, 5);
      Diag_tmp = diag(sum(Weight{i}));
      H{i} = Diag_tmp - Weight{i};
end
end

