function odtime = getodtime(d1, od_cut, times, timePromoterStart1)
odtime = zeros(size(d1,2),1);
for i=1:size(d1, 2)
    idx1 = find(d1(:,i) > od_cut);
    idx2 = find(times >= timePromoterStart1);
    odtime(i)=times(min(idx1(idx1>=min(idx2))));    
end

