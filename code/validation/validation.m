close all
global alpha gamma Lambda mu p2 rho sigma pa pe ps delta lambda k p1
% pre-difined parameters
p = load("parameters_CB.mat");
gamma = p.gamma; Lambda = p.Lambda; mu = p.mu; rho = p.rho;
sigma = p.sigma; pa = p.pa; pe = p.pe; ps = p.ps; 
p2 = 0.16;
for j = 1:3
    
x_real = date_cb(t_cb(1, 1)+tau2+ts:t_cb(end, 2)+tau2+ts);
x_sim = []; y_sig1 = []; y_sig2 = []; y_sim = [];

for i = 1:length(t_cb(:, 1))

    Q1 = Q(t_cb(i, 1):t_cb(i, 2));
    tp1 = tp(t_cb(i, 1):t_cb(i, 2));
    st1 = st(t_cb(i, 1):t_cb(i, 2));
    
    lambda = lambda_v(i);
    k = k_v(i);
    p1 = p1_v(i);
    if isnan(lambda)
        continue
    end
    
    step = 1;
    steps = 1:step:3*step;
    x0(1) = (Q1(3)-2*Q1(2)+Q1(1))/step^2; x0(1) = Q3(t_cb(i, 1));
    x0(2) = (Q1(2)-Q1(1))/step; x0(2) = Q2(t_cb(i, 1));
    x0(3) = Q1(1);
    
    st0(3) = st1(1)^2; 
    st0(2) = st1(2)^2+st1(1)^2; st0(2) = st2(t_cb(i, 1))^2;
    st0(1) = st1(3)^2+4*st1(2)^2+st1(1)^2; st0(1) = st3(t_cb(i, 1)) ^2;
    st0 = sqrt(st0);

    % prediction with the parameters
    ysigma = solution_sigma(st0, t_cb(i, 1), t_cb(i, 2));
    y2 = solution(x0, t_cb(i, 1), t_cb(i, 2));
    %t2 = t_cb(i, 1):t_cb(i,2);
    t2 = date_cb(t_cb(i, 1)+tau2+ts:t_cb(i, 2)+1+tau2+ts);
    
    
    if j == 1
        x_real = date_cb(t_cb(1, 1)+ts:t_cb(end, 2)+ts);
        t2 = date_cb(t_cb(i, 1)+ts:t_cb(i, 2)+1+ts);
        ysigma1 = y2+ysigma; 
        ysigma2 = y2-ysigma; 
        y_real = Q(1:length(x_real));
        y_sim = [y_sim; y2';nan];
        
    elseif j == 2
        ysigma1 = y2*alpha+y2*alpha.*sqrt((ysigma./y2).^2+(error_alpha/alpha)^2); 
        ysigma2 = y2*alpha-y2*alpha.*sqrt((ysigma./y2).^2+(error_alpha/alpha)^2); 
        y_real = dp(1:length(x_real));
        y_sim = [y_sim; y2'*alpha;nan];
    else
        ysigma1 = y2*delta+y2*delta.*sqrt((ysigma./y2).^2+(error_delta/delta)^2); 
        ysigma2 = y2*delta-y2*delta.*sqrt((ysigma./y2).^2+(error_delta/delta)^2); 
        y_real = r(1:length(x_real));
        y_sim = [y_sim; y2'*delta;nan];
    end
    x_sim = [x_sim; t2];
    ysigma1(ysigma1<0) = 0; ysigma1(isnan(ysigma1)) = 0;
    ysigma2(ysigma2<0) = 0; ysigma2(isnan(ysigma2)) = 0;
    y_sig1 = [y_sig1; ysigma1';nan]; y_sig2 = [y_sig2; ysigma2';nan];

end

figure_result(j, x_real, y_real, x_sim, y_sim, y_sig1, y_sig2)

end


function y = solution(x0, t0, tf)
    global alpha gamma mu p2 sigma pe ps delta lambda p1 k
    
    a1 = -(lambda*(1-p2)+gamma+delta+alpha+sigma+3*mu-k*pe);
    a2 = (-(alpha+delta+mu)*(lambda*(1-p2)+gamma+mu)-(sigma-k*pe+mu)*...
        (lambda*(1-p2)+gamma+alpha+delta+2*mu)+sigma*k*ps*(p2+(1-p2)*(1-p1*lambda)));
    a3 = (k*pe-sigma-mu)*(alpha+delta+mu)*(lambda*(1-p2)+gamma+mu)+...
        sigma*k*ps*(p2+(1-p2)*(1-p1*lambda))*(alpha+delta+mu);
    
    r = roots([1 -a1 -a2 -a3]);
    A = [exp(r(1)*t0) exp(r(2)*t0) exp(r(3)*t0); 
        exp(r(1)*t0)*r(1) exp(r(2)*t0)*r(2) exp(r(3)*t0)*r(3); 
        exp(r(1)*t0)*r(1)^2 exp(r(2)*t0)*r(2)^2 exp(r(3)*t0)*r(3)^2];
    B = [x0(3); x0(2); x0(1)];
    
    sol = linsolve(A, B);
    t = t0:tf;
    y = sol(1)*exp(r(1)*t)+sol(2)*exp(r(2)*t)+sol(3)*exp(r(3)*t);
    y(y<0) = 0;
end

function y = solution_sigma(x0, t0, tf)
    global alpha gamma mu p2 sigma pe ps delta lambda p1 k
    
    a1 = -(lambda*(1-p2)+gamma+delta+alpha+sigma+3*mu-k*pe);
    a2 = (-(alpha+delta+mu)*(lambda*(1-p2)+gamma+mu)-(sigma-k*pe+mu)*...
        (lambda*(1-p2)+gamma+alpha+delta+2*mu)+sigma*k*ps*(p2+(1-p2)*(1-p1*lambda)));
    a3 = (k*pe-sigma-mu)*(alpha+delta+mu)*(lambda*(1-p2)+gamma+mu)+...
        sigma*k*ps*(p2+(1-p2)*(1-p1*lambda))*(alpha+delta+mu);
    
    r = roots([1 -a1 -a2 -a3]);
    A = [1 1 1; 
        r(1) r(2) r(3); 
        r(1)^2 r(2)^2 r(3)^2];
    
    d = det(A);
    t = 1:length(t0:tf);
    a = x0(3)^2*((r(2)*r(3)^2-r(3)*r(2)^2).*exp(r(1)*t)+(r(3)*r(1)^2-r(1)*r(3)^2).*exp(r(2)*t)+...
    (r(1)*r(2)^2-r(2)*r(1)^2).*exp(r(3)*t)).^2/d^2;
    b = x0(2)^2*((r(2)^2-r(3)^2).*exp(r(1)*t)+(r(3)^2-r(1)^2).*exp(r(2)*t)+(r(1)^2-r(2)^2).*exp(r(3)*t)).^2./d^2;
    c = x0(1)^2*((r(3)-r(2)).*exp(r(1)*t)+(r(1)-r(3)).*exp(r(2)*t)+(r(2)-r(1)).*exp(r(3)*t)).^2./d^2;
    y = sqrt(a+b+c);
    y(y<0) = 0;
end

function figure_result(i,x_real, y_real, x_sim, y_sim, y_sig1, y_sig2)
    figure(i)
    hold on; grid on;
    scatter(x_real,y_real, '.', 'k')
    plot(x_sim,y_sim, '-', 'color', 'b');
    plot(x_sim,y_sig1, '--','color','r');
    plot(x_sim,y_sig2,'--','color','r');
    xlabel("time (days)")
    if i == 1
        ylabel("# people in quarantine")
        l = legend({'$\bar{Q}(k-\tau_2)$','$\bar{Q}(t-\tau_2)$',...
            '$\bar{Q}(t-\tau_2)\pm\sigma_{\bar{Q}}(t-\tau_2)$'});
        title("People in Quarantine per Day in the Period 17/07/2020-29/01/2021 ")
        %title("People in Quarantine per Day the Period 01/03/2020-16/07/2020 ")
        
    elseif i==2
        ylabel("# death people")
        l = legend({'$\bar{dp}(k-\tau_2)$','$\bar{dp}(t-\tau_2)$',...
            '$\bar{dp}(t-\tau_2)\pm\sigma_{\bar{dp}}(t-\tau_2)$'});
        title("Death People per Day in the Period 17/07/2020-29/01/2021 ")
        title("Death People per Day in the Period 01/03/2020-16/07/2020 ")
    else
        ylabel("# recovered people")
        l = legend({'$\bar{r}(k-\tau_2)$','$\bar{r}(t-\tau_2)$',...
            '$\bar{r}(t-\tau_2)\pm\sigma_{\bar{r}}(t-\tau_2)$'});
        title("Recovered People per Day in the Period 17/07/2020-29/01/2021 ")
        title("Recovered People per Day in the Period 01/03/2020-16/07/2020 ")
    end
    set(l,'interpreter','Latex')
end
