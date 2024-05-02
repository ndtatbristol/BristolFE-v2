clear
close all

leg_txt = {};

el_size = 1e-3;
% E1 = 210e9;
% rho1 = 8000;
% k1 = E1 * el_size; 

%Following approx for steel water
k1 = 1e11;
m1 = 1e-3;
k2 = 1e-3;
m2 = 1e-16;

%Following steel both sides
% k1 = 1e11;
% m1 = 1e-3;
% k2 = 1 / m1;
% m2 = 1 / k1;

w1 = sqrt(k1 / m1);
w2 = sqrt(k2 / m2);
crit_dt = min(2 ./ [w1, w2]);
dt = logspace(-3,3,101)' * crit_dt;



figure

%Uncoupled case first to show critical time step
uncoupled = zeros(numel(dt), 2);
for i = 1:numel(dt)
    uncoupled(i, 1) = fn_max_EV_of_simple_interface(k1, m1, 0, k1, m1, 0, dt(i), 0);
    uncoupled(i, 2) = fn_max_EV_of_simple_interface(k2, m2, 0, k2, m2, 0, dt(i), 0);
end
leg_txt(end + 1) = {'Uncoupled1'};
leg_txt(end + 1) = {'Uncoupled2'};
loglog(dt / crit_dt, uncoupled); hold on;

% %With normal interface element
% basic_interface = zeros(size(dt));
% for i = 1:numel(dt)
%     basic_interface(i) = fn_max_EV_of_simple_interface(k1, m1, 0, k2, m2, 0, dt(i), 1);
% end
% loglog(dt / crit_dt, basic_interface); hold on;
% leg_txt(end + 1) = {'Basic interface'};

%With normal interface element and flipped sign on one side
basic_interface = zeros(size(dt));
for i = 1:numel(dt)
    basic_interface(i) = fn_max_EV_of_simple_interface(k1, m1, 0, -k2, -m2, 0, dt(i), 1);
end
loglog(dt / crit_dt, basic_interface, 'r.'); hold on;
leg_txt(end + 1) = {'Basic interface flip'};

% %With normal interface element and flipped sign on other side
% basic_interface = zeros(size(dt));
% for i = 1:numel(dt)
%     basic_interface(i) = fn_max_EV_of_simple_interface(-k1, -m1, 0, k2, m2, 0, dt(i), 1);
% end
% loglog(dt / crit_dt, basic_interface, 'ro'); hold on;
% leg_txt(end + 1) = {'Basic interface flip2'};

% With reduced coupling interface element
% fixed_reduced_coupling = zeros(size(dt));
% for i = 1:numel(dt)
%     fixed_reduced_coupling(i) = fn_max_EV_of_simple_interface(k1, m1, 0, -k2, -m2, 0, dt(i), 0.1);
% end
% loglog(dt / crit_dt, fixed_reduced_coupling); hold on;
% leg_txt(end + 1) = {'Reduced coupling'};

% % Added damping based on material on same side
% c1 = 8 * sqrt(k1 * m1);
% c2 = 8 * sqrt(k2 * m2);
% max_EV_same_damping = zeros(size(dt));
% for i = 1:numel(dt)
%     max_EV_same_damping(i) = fn_max_EV_of_simple_interface(k1, m1, c1, -k2, -m2, -c2, dt(i), 1);
% end
% loglog(dt / crit_dt, max_EV_same_damping); hold on;
% leg_txt(end + 1) = {'Fixed same damping'};

% %Added damping based on material on other side (radiation into other medium) -
% %this seems to be correct!
% c1 = 8 / sqrt(k2 * m2);
% c2 = 8 / sqrt(k1 * m1);
% fixed_opp_damping = zeros(size(dt));
% for i = 1:numel(dt)
%     fixed_opp_damping(i) = fn_max_EV_of_simple_interface(k1, m1, c1, k2, m2, c2, dt(i), 1);
% end
% loglog(dt / crit_dt, fixed_opp_damping); hold on;
% leg_txt(end + 1) = {'Fixed opp. damping'};

legend(leg_txt);

