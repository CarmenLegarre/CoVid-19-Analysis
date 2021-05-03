all_data = load("range2.mat");
clear p;
all_data = all_data.range;
collected_data = [];
figure(5)
hold on;
error_m = [];
for i = 1:length(all_data)
    data = all_data{i};
    intervals(i) = {find_interval([tp, test], data)};
    collected_data = [collected_data; data];
    [fitobject,gof1]=fit(data(:, 1),data(:, 2),'poly1');
    coeff = coeffvalues(fitobject);
    error = confint(fitobject);
    error_m = [error_m; diff(error(:, 1))/2 diff(error(:, 2))/2];
    x = min(data(:, 1)):max(data(:, 1));
    hold on;
    plot(x, coeff(1)*x+coeff(2), "b", "LineWidth", 1, 'DisplayName','This one');
    p(i) = plot(x, coeff(1)*x+coeff(2), "b", "LineWidth", 1, 'DisplayName','This one');
    l(i) = "y_"+num2str(i)+"="+num2str(round(coeff(1)*100)/100)+"x_"+num2str(i)+"+"+num2str(round(coeff(2))); 
    p2 = scatter(data(:, 1), data(:, 2), ".", 'MarkerEdgeColor',[1 0.411764705882353 0.16078431372549],"SizeData", 80);
end

%not_collected_data = delete(not_collected_data, [x y]);
not_collected_data = delete([tp, test], collected_data);

p3=scatter(not_collected_data(:, 1), not_collected_data(:, 2),".", 'MarkerEdgeColor',[0 0 0],"SizeData", 80);

legend([p, p2, p3],[l(1:end), "collected data", "not collected data"])


display(error_m)



function new_data = delete(data, data_e)
new_data1 = data(:, 1);
new_data2 = data(:, 2);
new_data3 = data(:, 2)./data(:, 1);
    for i = 1:length(data_e(:, 1))
        [v, index] = min(abs(new_data3-data_e(i, 2)/data_e(i, 1)));
        new_data1(index) = 0;
        new_data2(index) = 0;
    end

new_data = [new_data1(new_data1~=0) new_data2(new_data2~=0)];
length(new_data)
end

function interval = find_interval(data, data_f)

    new_data1 = data(:, 1);
    new_data2 = data(:, 2);
    new_data3 = data(:, 2)./data(:, 1);
    for i = 1:length(data_f(:, 1))
        [v, index] = min(abs(new_data3-data_f(i, 2)/data_f(i, 1)));
        interval(i) = index;
    end

end