% % Double tank
 close all;
 clear all;
 [X,U,J_array,U_star,T] = double_tank_problem_proj();
 global r;
 plot(T,U+1,T,U_star+1);
 title('Relaxed Control and U\_star');
 legend('U','U*')
 saveas(gcf,'1_relaxed_control.png');
 plot(T,X(2,:),T,r,T,X(1,:));
 title('State vs Time');
 legend('x2','ref','x1');
 saveas(gcf,'1_state_vs_time.png');
 
 plot(J_array);
 title('Cost vs iteration');
 saveas(gcf,'1_cost.png');


% Robot path planning 
 close all;
 clear all;
 [X,U,J_array,U_star,T] = robot_path_planning();
 plot(T,U(1,:),T,U(2,:),T,1-U(1,:)-U(2,:));
 title('Relaxed mode');
 legend('fg','fcl','fcc')
 saveas(gcf,'2_relaxed_mode.png');
 plot(T,U_star(1,:),T,U_star(2,:),T, 1-U_star(1,:)-U_star(2,:));
 title('Optimal mode');
 legend('fg','fcl','fcc')
 saveas(gcf,'2_optimal_mode.png');
 plot(X(1,:),X(2,:));
 title('State Trajectory');
 
 saveas(gcf,'2_state_traj.png');
 
 plot(J_array);
 title('Cost vs iteration');
 saveas(gcf,'2_cost.png');

close all;
clear all;

[X,U,J_array,U_star,T] = poweraware_pathplanning(1)

plot(U(1,:),U(2,:),U(3,:),U(4,:),U(5,:),U(6,:),U(7,:),U(8,:),U(9,:),U(10,:),U(11,:),U(12,:));
title('Phase Potrait: Velocities');
legend('agent1','agent2','agent3','agent4','agent5','agent6')
saveas(gcf,'3_C1_pase_potrait_vel.png');
plot(X(1,:),X(2,:),X(3,:),X(4,:),X(5,:),X(6,:),X(7,:),X(8,:),X(9,:),X(10,:),X(11,:),X(12,:));
title('Phase Potrait: Postion');
legend('agent1','agent2','agent3','agent4','agent5','agent6')
saveas(gcf,'3_C1_pase_potrait_pos.png');
plot(J_array);
title('Cost vs iter');
saveas(gcf,'3_C1_cost.png');


close all;
clear all;

[X,U,J_array,U_star,T] = poweraware_pathplanning(7)

plot(U(1,:),U(2,:),U(3,:),U(4,:),U(5,:),U(6,:),U(7,:),U(8,:),U(9,:),U(10,:),U(11,:),U(12,:));
title('Phase Potrait: Velocities');
legend('agent1','agent2','agent3','agent4','agent5','agent6')
saveas(gcf,'3_C7_pase_potrait_vel.png');
plot(X(1,:),X(2,:),X(3,:),X(4,:),X(5,:),X(6,:),X(7,:),X(8,:),X(9,:),X(10,:),X(11,:),X(12,:));
title('Phase Potrait: Postion');
legend('agent1','agent2','agent3','agent4','agent5','agent6')
saveas(gcf,'3_C7_pase_potrait_pos.png');
plot(J_array);
title('Cost vs iter');
saveas(gcf,'3_C7_cost.png');
