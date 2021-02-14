ranseq = randi([0 1],500,1);
re = ranseq*0.75+(1-ranseq)*0.25;

consumption = zeros(1,500);
asset = zeros(1,500);
b = 0;
for i=1:500
w = b + re(i);
    if w <= 2
        w = round(w,3);
        [~,id] = min(abs(x-w));
        b = g(id);
        asset(i) = b;
        consumption(i) = w-(b*0.97);
    else % the wealth exceeds the bound of v function,  stop the loop
        asset(i) = -1;
        consumption(i) = -1;
        b = 2;
    end
end
plot(consumption)
plot(asset)