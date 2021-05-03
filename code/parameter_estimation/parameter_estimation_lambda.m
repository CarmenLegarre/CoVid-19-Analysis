global gamma Lambda mu p2 rho sigma pa pe ps n

% pre-difined parameters
p = load("parameters_CB.mat");
gamma = p.gamma; Lambda = p.Lambda; mu = p.mu; rho = p.rho;
sigma = p.sigma; pa = p.pa; pe = p.pe; ps = p.ps; 
p2 = 0.16; n = 0;

lambda_v = zeros(length(t_cb), 1);
p1_v = zeros(length(t_cb), 1);
k_v = zeros(length(t_cb), 1);

for i = 1: length(t_cb(:, 1))
    % data range to estimate the parameters
    Q1 = Q(t_cb(i, 1):t_cb(i, 2));
    tp1 = tp(t_cb(i, 1):t_cb(i, 2));
    
    % estimate the parameters of the ARX model
    [a(i, :), a_error, a_variance] = param_est(Q1, 3);
    [b(i, :), b_error, b_variance] = param_est(tp1, 2);
 
    % estimate the parameters lambda, k and p1 from ARX model parameters
    [lambda, k, p1] = solution_equation(a(i, :), b(i, :));

    % save data
    lambda_v(i) = lambda;
    k_v(i) = k;
    p1_v(i) = p1;
end



function a = variance(param_values1, param_variances1, param_values2, param_variances2)
    param_variances = zeros(4);
    param_variances(1:3, 1:3) = param_variances1(:, :);
    param_variances(4, 4) = param_variances2(1, 1);
    
    param_values=[param_values1';param_values2'];
    
    J = zeros(3, 4);
    h = 1e-5;
    
    for i=1:length(J(:, 1))
        
        [lambda1, k1, p11] = solution_equation4(param_values1, param_values2);
        
        for j = 1:length(J(1, :))
            
            param_values1(j) = param_values(j)+h;
            [lambda2, k2, p12] = solution_equation4(param_values(1:3), param_values(4));
            param_values(j) = param_values(j)-h;
            
            if i == 1
                J(i, j) = (lambda2-lambda1)/h;
            elseif i == 2
                J(i, j) = (k2-k1)/h;
            elseif i == 3
                J(i, j) = (p12-p11)/h;
            end
            
        end
        
    end
    J = J(1:3, 1:3);
    param_variances = param_variances(1:3, 1:3);
    a = J*param_variances*J';
end

function [a, param_error, param_variance] = param_est(data, order)  
    phi = zeros(length(data(order:end-1)),order);
    y = zeros(length(data(order+1:end)), 1);

    for i = 1:order
        phi(:, i) = -data(order-i+1:end-i);
    end
    
    y(:) = data(order+1:end);
    a = inv(phi'*phi)*phi'*y;
    
    % One step ahead prediction
    est_output = zeros(1, length(data));
    for i = 1:order
        est_output(i) = data(i);
    end
    
    for i = order+1:length(data)
        for j = 1:order
            est_output(i) = est_output(i) - data(i-j)*a(j);
        end
    end
    
    variance = param_est_error(data, est_output', order);
    param_error = zeros(order, 1);
    param_variance = zeros(order, 1);
    
    for i = 1:order
        param_variance(i) = variance(i, i);
        param_error(i) = sqrt(variance(i, i))/(sqrt(length(data)))/length(data);
    end
    param_variance = variance;
end

function a = param_est_error(output, estimated_output, order)
    phi = zeros(length(output(order:end-1)), order);
    for i = 1:order
        phi(:, i) = -output(order-i+1:end-i);
    end
    R = (phi'*phi);
    v = var(output-estimated_output);
    a = v*inv(R);
    
end


function [lambda_v, k_v, p1_v] = solution_equation(a, b)
    global gamma p2 sigma pe ps 
    
    ranges = 0:0.001:2;
    for i = 1:length(ranges)
        
        [lambda, k, p1] = solution_equation1(a, ranges(i));
        c2 = sigma*k*ps*(p2+(1-p2)*(1-p1*lambda))+(lambda*(1-p2)+gamma)*(k*pe-sigma);
        c1 = -(sigma+(lambda*(1-p2)+gamma)-k*pe);
        f1 = 2+c1+b(1)*(lambda*(1-p2)*(1+ps*p1*k));
        f2 = c2-c1-1+b(2)*(lambda*(1-p2)*(1+ps*p1*k))-1;
        f = sigma+(lambda*(1-p2)+gamma)-k*pe-2-b(1)*(lambda*(1-p2)*...
            (1+ps*p1*k));
        f = 2*ps*k*lambda*(1-p2)*(b(1)*f1+b(2)*f2);
        if i == 1
            x1 = ranges(i);
            y1 = f;
        elseif i > 2
            x2 = ranges(i);
            y2 = f;
        end
        
        if i > 2
            if y1*y2<0
                break
            end
        end
            
    end
    
    condition = y2*y1 < 0;
    
    if ~condition
        lambda_v = nan;
        k_v = nan;
        p1_v = nan;
    end
    
    i = 0;
    while condition
        i = i+1;
        x = x1 + (x2-x1)/2;
        [lambda, k, p1] = solution_equation1(a, x);
        c2 = sigma*k*ps*(p2+(1-p2)*(1-p1*lambda))+(lambda*(1-p2)+gamma)*(k*pe-sigma);
        c1 = -(sigma+(lambda*(1-p2)+gamma)-k*pe);
        f1 = 2+c1+b(1)*(lambda*(1-p2)*(1+ps*p1*k));
        f2 = c2-c1-1+b(2)*(lambda*(1-p2)*(1+ps*p1*k))-1;
        
        f = sigma+(lambda*(1-p2)+gamma)-k*pe-2-b(1)*(lambda*(1-p2)*...
            (1+ps*p1*k));
        f = 2*ps*k*lambda*(1-p2)*(b(1)*f1+b(2)*f2);
        if y1*f > 0
            y1 = f;
            x1 = p1;
        elseif y1*f < 0
            y2 = f;
            x2 = p1;
        end
        
        
        if abs(f)<1e-20
            lambda_v = lambda;
            k_v = k;
            p1_v = p1;
            condition = false;
        elseif abs(f)<1e-6 && (i>100 || i>99)
            lambda_v = lambda;
            k_v = k;
            p1_v = p1;
            condition = false;
            
        elseif (x<0 && x>1) || i>100
            abs(f)
            lambda_v = nan;
            k_v = nan;
            p1_v = nan;
            condition = false;
        end
        
    end

end

function [lambda_v, k_v, p1_v] = solution_equation1(a, ranges_p1)
    global gamma p2 sigma pe ps alpha delta n
    
    ranges_lambda = 0:0.1:1.5;
    
    for i = 1:length(ranges_p1)
        
        p1 = ranges_p1(i);
        
        for j = 1:length(ranges_lambda)
            
            lambda = ranges_lambda(j);
            k = (lambda*(1-p2) + (gamma+alpha+delta+sigma) -a(1) - 3)/pe;
            
            f1 = lambda*(1-p2)*(alpha+delta+sigma)+k*...
                (-pe*(gamma+alpha+delta)-sigma*ps)+k*lambda*(1-p2)*(sigma*ps*p1*(n+1)-pe)-...
                a(2)-2*a(1)-3+(alpha+delta)*gamma+sigma*(gamma+alpha+delta);
            f2 = -lambda*(1-p2)*sigma+k*(pe*gamma+sigma*ps+lambda*(1-p2)*...
                (pe-sigma*ps*p1*(n+1)))+(1+a(2)+a(3)+a(1))/(alpha+delta)-sigma*gamma;
            
            a1 = -pe*(gamma+alpha+delta)-sigma*ps+lambda*(1-p2)*(-pe+sigma*ps*p1*(n+1));
            a2 = pe*gamma+sigma*ps+lambda*(1-p2)*(pe-sigma*ps*p1*(n+1));
            
            f = 2*a1*f1+2*a2*f2;
            if j == 1
                x1 = ranges_lambda(j);
                y1 = f;
            elseif j > 2
                x2 = ranges_lambda(j);
                y2 = f;
            end
            
            if j > 2
                if y1*y2<0
                    break
                end
            end
        
        end
            
    
    condition = y2*y1 < 0;
    
    if ~condition
        lambda_v = nan;
        k_v = nan;
        p1_v = nan;
    end
    
    ii = 0;
    while condition
        ii = ii+1;
        x = x1 + (x2-x1)/2;
        lambda = x;
        
        k = (lambda*(1-p2) + (gamma+alpha+delta+sigma) -a(1) - 3)/pe;
            
        f1 = lambda*(1-p2)*(alpha+delta+sigma)+k*...
            (-pe*(gamma+alpha+delta)-sigma*ps)+k*lambda*(1-p2)*(sigma*ps*p1*(n+1)-pe)-...
            a(2)-2*a(1)-3+(alpha+delta)*gamma+sigma*(gamma+alpha+delta);
        f2 = -lambda*(1-p2)*sigma+k*(pe*gamma+sigma*ps+lambda*(1-p2)*...
            (pe-sigma*ps*p1*(n+1)))+(1+a(2)+a(3)+a(1))/(alpha+delta)-sigma*gamma;
        
        a1 = -pe*(gamma+alpha+delta)-sigma*ps+lambda*(1-p2)*(-pe+sigma*ps*p1*(n+1));
        a2 = pe*gamma+sigma*ps+lambda*(1-p2)*(pe-sigma*ps*p1*(n+1));
        
        f = 2*a1*f1+2*a2*f2;
        if y1*f > 0
            y1 = f;
            x1 = lambda;
        elseif y1*f < 0
            y2 = f;
            x2 = lambda;
        end
        
        if abs(f)<1e-14
            lambda_v(i) = lambda;
            k_v(i) = k;
            p1_v(i) = p1;
            condition = false;
            
        elseif (x<0 && x>1.5) || ii>100
            lambda_v = nan;
            k_v = nan;
            p1_v = nan;
            condition = false;
        end
        
    end
    
    end

end

function j = jacobian(a, b, lambda, k, p1)
global gamma p2 sigma pe ps alpha delta
f1 = lambda*(1-p2)*(alpha+delta+sigma)+k*...
            (-pe*(gamma+alpha+delta)-sigma*ps)+k*lambda*(1-p2)*(sigma*ps*p1-pe)-...
            a(2)-2*a(1)-3+(alpha+delta)*gamma+sigma*(gamma+alpha+delta);
        
f2 = -lambda*(1-p2)*sigma+k*(pe*gamma+sigma*ps+lambda*(1-p2)*...
            (pe-sigma*ps*p1))+(1+a(2)+a(3)+a(1))/(alpha+delta)-sigma*gamma;
       
j = [-1 0 0 0; -4*f1+2*f2/(alpha+delta) -2*f1+f2/(alpha+delta) f2/(alpha+delta) 0;
    0 0 0 -lambda*(1-p2)*1+ps*p1*lambda*k];

end

