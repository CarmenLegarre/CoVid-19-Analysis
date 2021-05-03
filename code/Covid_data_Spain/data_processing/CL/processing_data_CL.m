clear all;
[tests_date, tests, tp_date, tp, ph_date, ph, pi_date, pi, dp_date, dp, r_date, r, n_ph, n_uci] = generate_data_CL();
createfigure1(tp_date, tp);
createfigure2(tp_date, tp./tests);
createfigure3(ph_date, [ph(:), pi(:), dp(:)]);
createfigure4([tests_date, tp_date], [tests, tp], ...
    [ph_date, pi_date, dp_date], [ph, pi, dp]);