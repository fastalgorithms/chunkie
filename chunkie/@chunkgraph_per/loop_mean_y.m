function y = loop_mean_y(obj,edges)
%CURVE_MEAN_Y mean y-coordinate of the points making up an edge list.

s = 0;
cnt = 0;

for jj = 1:numel(edges)
    rr = obj.echnks(abs(edges(jj))).r(:,:);
    s = s + sum(rr(2,:));
    cnt = cnt + size(rr,2);
end

y = s/cnt;
end
