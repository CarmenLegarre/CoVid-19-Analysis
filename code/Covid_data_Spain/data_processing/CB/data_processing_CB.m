%clear all
[reduced_dates, test, tp, ph, pi, dp, r] = generate_data();
dim = length(reduced_dates);
createfigure1(reduced_dates, tp);
positivity_ratio = tp./test;
createfigure2(reduced_dates(10:dim), positivity_ratio(10:dim));
createfigure3(reduced_dates, [ph(:), pi(:), dp(:)]);
createfigure4(reduced_dates, [test(:), tp(:)], [ph(:), pi(:), dp(:)],...
    [test(:), tp(:)]);
