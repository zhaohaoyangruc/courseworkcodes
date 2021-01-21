function w = reserwage(d,l,u,b,g)
% w: reservation wage to be reported
% d: discount factor
% l: lower bound of wage offer
% u: upper bound of wage offer
% b: unemployment benefit
% g: initial guess of w
fun = @(x,c) max(x,c);
w0 = g;
w1 = g+0.1; % just initialization
while abs(w1-w0)>10^(-6)
    w0 = w1;
    w1 = (1-d)*b+(d./(u-l))*integral(@(x) fun(x,w0),l,u);
end
w = w1;