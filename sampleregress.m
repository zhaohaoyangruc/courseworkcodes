function b = sampleregress(n,r)
    % n: sample size
    % r: repeat times
    % b: the vector of estimated coefficients
    
    b = zeros(r,1); % initialize result vector
    
    for i=1:r
        % data generation
        x = -1 + 2*rand(n,1); % generate x
        u = lognrnd(0,1,n,1); % generate u
        y = x + u; % calculate y
        
        % regression
        mdl = fitlm(x,y,'Intercept',false); % regression y on x, without intercept
        b(i) = mdl.Coefficients.Estimate;
    end
    
    b = sqrt(n).*(b-1); % CLT expression
end


