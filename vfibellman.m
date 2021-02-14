function [v,g] = vfibellman(a,b,s,alphap,betap,q,x)
% a: the lower bound of the wealth 
% b: the upper bound of the wealth
% s: INTEGER, 10^(-s) is the step length
% q: bond price
% x: initial guess of bellman equation. row vector with length
%       1+(b-a)/s.
p = 10^(-1*s); % step length
w = a:p:b; % vector of possible wealth level
n = length(w); 
v = zeros(1,n); % Initialize bellman equation result vector
g = zeros(1,n); % Initialize policy function result vector
maxc = zeros(1,n); % Initialize maximizer result vector

er = 1; % Initialize error 
itr = 0; % Iteration counter;
while er>10^(-3)
    itr = itr + 1; 
    for i=1:n
        c = p:p:w(i); % possible consumption level
        % consumption level cannot be negative or 0 (log utility)
        % consumption level cannot be higher than wealth (no borrow)
        m = length(c);
        wi = w(i).*ones(1,m); % row vector of w(i)
        fw1 = (wi-c)/q+alphap; % the wealth tomorrow, alpha realized
        fw1 = round(fw1,s); % round to a level in w
        fw2 = (wi-c)/q+1-alphap; % the wealth tomorrow, 1-alpha realized
        fw2 = round(fw2,s); % round to a level in w
        
        y = zeros(1,m); % Initialize intermediate bellman value vector
        for j=1:m % calculate bellman expression for given consumption level
            if fw1(j) > b % future wealth level exceed the upper bound
                    id1 = -1; % revise to an indicator
            else 
                [~,id1] = min(abs(w-fw1(j))); % choose the wealth level closest
            end
            
            if fw2(j) > b % future wealth level exceed the upper bound
                    id2 = -1; % revise to an indicator
            else 
                [~,id2] = min(abs(w-fw2(j))); % choose the wealth level closest
            end
            
            if (id1 == -1) || (id2 == -1)
                y(j) = log(c(j))-1000; % Too small consumption cannot be optimal
            else
                y(j) = log(c(j)) + betap*0.5*x(id1) + betap*0.5*x(id2);
            end
        end
        % choose the consumption level that maximize bellman
        [v(i),id3] = max(y); % find the maximizer location
        maxc(i) = c(id3); % record the maximizer
        g(i) = (w(i)-maxc(i))/q; % policy function is given by budget constraint
    end
    diff = abs(v-x);
    er = max(diff); % update error level
    disp(er)
    disp(itr)
    x = v;
end