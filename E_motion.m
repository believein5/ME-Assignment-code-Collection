clear all; clc; close all;

% 初始参数
K = 1.48;
deg2rad = pi/180;
rad2deg = 180/pi;



% 长度 mm
L_AC = 370; % A点在C点正上方
L_AB = 110.63;
L_CD = 534.45; 
L_DE = L_CD * 0.26;
CE_Vertical = 550; % C点到滑块E的垂直距离

% 角度范围
theta_deg = linspace(0, 360, 1000);

% 预分配phi数组
phi = zeros(size(theta_deg));   % CD与竖直夹角
% 初始化v_E数组
v_E = zeros(size(theta_deg));
% 预分配加速度数组
a_E = zeros(size(theta_deg));

% 速度参数
A_rpm = 56;
w_A = A_rpm * 2 * pi / 60;  % 转换为弧度/秒
v_B = w_A * L_AB;

% 加速度参数
a_B_n = w_A^2 * L_AB;
a_B_t = 0;
a_B = sqrt(a_B_n^2 + a_B_t^2);

% 时间参数
T = 60 / A_rpm;  % 曲柄转动周期（秒）
dt = T / (length(theta_deg) - 1);  % 时间步长


% 预分配位置数组
A_pos = zeros(length(theta_deg), 2);  % A点位置（固定）
B_pos = zeros(length(theta_deg), 2);  % B点位置
C_pos = zeros(length(theta_deg), 2);  % C点位置（固定）
D_pos = zeros(length(theta_deg), 2);  % D点位置
E_pos = zeros(length(theta_deg), 2);  % E点位置

% 设置固定点
A_pos = repmat([0, 0], length(theta_deg), 1);
C_pos = repmat([0, -L_AC], length(theta_deg), 1);



% 循环计算每个theta对应的phi和v_E
for i = 1:length(theta_deg)
    theta = theta_deg(i);
    theta_rad = theta * deg2rad;


    % B点坐标（曲柄末端）
    B_pos(i,:) = [L_AB * -sind(theta), L_AB * -cosd(theta)];


    % 求phi
    phi(i) = atand(sind(theta)/(L_AC/L_AB-cosd(theta)));
    phi_rad = phi(i) * deg2rad;


    % D点坐标（摇杆末端）
    % 注意：phi是CD与竖直方向的夹角
    D_pos(i,:) = [C_pos(i,1) - L_CD * sin(phi_rad), ...
                  C_pos(i,2) + L_CD * cos(phi_rad)];
    
    % gamma 连杆DE与水平面的夹角
    gamma = asind((CE_Vertical - L_CD*cosd(phi(i)))/L_DE);
    gamma_rad = gamma * deg2rad;    

    % 计算B点分速度
    v_B_CD = v_B * cos(deg2rad * (phi(i) + theta));
    % B点到C点的距离（r_B）
    r_B = sqrt(L_AC^2 + L_AB^2 - 2 * L_AC * L_AB * cosd(theta));

    % E点坐标（滑块）
    % 假设E点沿水平导轨运动
    E_pos(i,:) = [D_pos(i,1) - L_DE * cos(gamma_rad), ...
                  C_pos(i,2) + CE_Vertical];  % E点在固定高度
    
    % 计算CD杆的角速度
    w_CD = v_B_CD / r_B; 
    
    % D点速度
    v_D = w_CD * L_CD; 
    

    % E点速度
    v_E(i) = v_D / sind(90-gamma) * sind(90 + gamma -phi(i)); % 正弦定理

end

% 通过数值微分计算加速度
% 使用中心差分法
for i = 2:length(v_E)-1
    a_E(i) = (v_E(i+1) - v_E(i-1)) / (2*dt);
end
% 边界处理
a_E(1) = (v_E(2) - v_E(1)) / dt;
a_E(end) = (v_E(end) - v_E(end-1)) / dt;

% 或者使用MATLAB的gradient函数（更简洁）
% a_E = gradient(v_E, dt);


% ============================
% 创建动画
% ============================

figure('Position', [100, 100, 1200, 800]);

% 子图1：机构运动动画
subplot(2, 3, [1, 2, 4, 5]);
hold on;
axis equal;
grid on;
box on;
xlabel('X (mm)', 'FontSize', 12);
ylabel('Y (mm)', 'FontSize', 12);
title('杆组机构运动仿真', 'FontSize', 14);

% 设置坐标轴范围
x_min = min([min(B_pos(:,1)), min(D_pos(:,1)), min(E_pos(:,1))]) - 50;
x_max = max([max(B_pos(:,1)), max(D_pos(:,1)), max(E_pos(:,1))]) + 50;
y_min = min([min(B_pos(:,2)), min(D_pos(:,2)), min(E_pos(:,2))]) - 50;
y_max = max([max(B_pos(:,2)), max(D_pos(:,2)), max(E_pos(:,2))]) + 50;
axis([x_min, x_max, y_min, y_max]);

% 绘制固定点和导轨
plot(A_pos(1,1), A_pos(1,2), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
plot(C_pos(1,1), C_pos(1,2), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
% 绘制E点导轨（水平线）
plot([x_min, x_max], [CE_Vertical, CE_Vertical], 'k--', 'LineWidth', 1);

% 添加文字标签
text(A_pos(1,1), A_pos(1,2)-20, 'A', 'FontSize', 12, 'HorizontalAlignment', 'center');
text(C_pos(1,1), C_pos(1,2)+20, 'C', 'FontSize', 12, 'HorizontalAlignment', 'center');
text(x_max-100, CE_Vertical+20, 'E点导轨', 'FontSize', 10);

% 预分配图形对象句柄
h_AB = line([A_pos(1,1), B_pos(1,1)], [A_pos(1,2), B_pos(1,2)], ...
            'Color', 'b', 'LineWidth', 3);
h_BC = line([B_pos(1,1), C_pos(1,1)], [B_pos(1,2), C_pos(1,2)], ...
            'Color', 'g', 'LineWidth', 2, 'LineStyle', '--');
h_CD = line([C_pos(1,1), D_pos(1,1)], [C_pos(1,2), D_pos(1,2)], ...
            'Color', 'r', 'LineWidth', 3);
h_DE = line([D_pos(1,1), E_pos(1,1)], [D_pos(1,2), E_pos(1,2)], ...
            'Color', 'm', 'LineWidth', 2);

% 绘制各点
h_B = plot(B_pos(1,1), B_pos(1,2), 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b');
h_D = plot(D_pos(1,1), D_pos(1,2), 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'g');
h_E = plot(E_pos(1,1), E_pos(1,2), 'ks', 'MarkerSize', 10, 'MarkerFaceColor', 'y');

% 添加B、D、E点标签
h_B_text = text(B_pos(1,1), B_pos(1,2)+15, 'B', 'FontSize', 10);
h_D_text = text(D_pos(1,1), D_pos(1,2)+15, 'D', 'FontSize', 10);
h_E_text = text(E_pos(1,1), E_pos(1,2)+15, 'E', 'FontSize', 10);

% 轨迹跟踪
h_B_trace = plot(B_pos(1,1), B_pos(1,2), 'b:', 'LineWidth', 0.5);
h_D_trace = plot(D_pos(1,1), D_pos(1,2), 'r:', 'LineWidth', 0.5);
h_E_trace = plot(E_pos(1,1), E_pos(1,2), 'k-', 'LineWidth', 1);

% 信息显示
info_text = text(x_min+50, y_max-50, '', 'FontSize', 11, 'BackgroundColor', 'w');

% 子图2：速度曲线（实时更新）
subplot(2, 3, 3);
h_vel = plot(theta_deg(1), v_E(1), 'b-', 'LineWidth', 2);
xlabel('Theta (degrees)', 'FontSize', 10);
ylabel('Velocity v_E (mm/s)', 'FontSize', 10);
title('速度曲线', 'FontSize', 12);
grid on;
xlim([0, 360]);
ylim([min(v_E)*1.1, max(v_E)*1.1]);
hold on;
h_current_vel = plot(theta_deg(1), v_E(1), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');

% 子图3：加速度曲线（实时更新）
subplot(2, 3, 6);
h_acc = plot(theta_deg(1), a_E(1), 'r-', 'LineWidth', 2);
xlabel('Theta (degrees)', 'FontSize', 10);
ylabel('Acceleration a_E (mm/s²)', 'FontSize', 10);
title('加速度曲线', 'FontSize', 12);
grid on;
xlim([0, 360]);
ylim([min(a_E)*1.1, max(a_E)*1.1]);
hold on;
h_current_acc = plot(theta_deg(1), a_E(1), 'mo', 'MarkerSize', 8, 'MarkerFaceColor', 'm');

% 动画循环
animation_speed = 2;  % 动画速度因子
for i = 1:animation_speed:length(theta_deg)
    % 更新主图
    subplot(2, 3, [1, 2, 4, 5]);
    
    % 更新杆件位置
    set(h_AB, 'XData', [A_pos(i,1), B_pos(i,1)], ...
              'YData', [A_pos(i,2), B_pos(i,2)]);
    set(h_BC, 'XData', [B_pos(i,1), C_pos(i,1)], ...
              'YData', [B_pos(i,2), C_pos(i,2)]);
    set(h_CD, 'XData', [C_pos(i,1), D_pos(i,1)], ...
              'YData', [C_pos(i,2), D_pos(i,2)]);
    set(h_DE, 'XData', [D_pos(i,1), E_pos(i,1)], ...
              'YData', [D_pos(i,2), E_pos(i,2)]);
    
    % 更新各点位置
    set(h_B, 'XData', B_pos(i,1), 'YData', B_pos(i,2));
    set(h_D, 'XData', D_pos(i,1), 'YData', D_pos(i,2));
    set(h_E, 'XData', E_pos(i,1), 'YData', E_pos(i,2));
    
    % 更新标签位置
    set(h_B_text, 'Position', [B_pos(i,1), B_pos(i,2)+15]);
    set(h_D_text, 'Position', [D_pos(i,1), D_pos(i,2)+15]);
    set(h_E_text, 'Position', [E_pos(i,1), E_pos(i,2)+15]);
    
    % 更新轨迹
    set(h_B_trace, 'XData', B_pos(1:i,1), 'YData', B_pos(1:i,2));
    set(h_D_trace, 'XData', D_pos(1:i,1), 'YData', D_pos(1:i,2));
    set(h_E_trace, 'XData', E_pos(1:i,1), 'YData', E_pos(1:i,2));
    
    % 更新信息文本
    info_str = sprintf('θ = %.1f°\nφ = %.1f°\nv_E = %.1f mm/s\na_E = %.1f mm/s²', ...
                      theta_deg(i), phi(i), v_E(i), a_E(i));
    set(info_text, 'String', info_str);
    
    % 更新速度曲线
    subplot(2, 3, 3);
    set(h_vel, 'XData', theta_deg(1:i), 'YData', v_E(1:i));
    set(h_current_vel, 'XData', theta_deg(i), 'YData', v_E(i));
    
    % 更新加速度曲线
    subplot(2, 3, 6);
    set(h_acc, 'XData', theta_deg(1:i), 'YData', a_E(1:i));
    set(h_current_acc, 'XData', theta_deg(i), 'YData', a_E(i));
    
    % 更新标题显示当前角度
    subplot(2, 3, [1, 2, 4, 5]);
    title(sprintf('杆组机构运动仿真 (θ = %.1f°)', theta_deg(i)), 'FontSize', 14);
    
    % 刷新图形
    drawnow;
    
    % 控制动画速度
    pause(0.01);
end

% 显示最终结果
fprintf('=== 运动分析结果 ===\n');
fprintf('曲柄转速: %.1f rpm (ω = %.3f rad/s)\n', A_rpm, w_A);
fprintf('v_E 范围: %.2f 到 %.2f mm/s\n', min(v_E), max(v_E));
fprintf('v_E 平均值: %.2f mm/s\n', mean(v_E));
fprintf('a_E 范围: %.2f 到 %.2f mm/s²\n', min(a_E), max(a_E));
fprintf('a_E 平均值: %.2f mm/s²\n', mean(a_E));
fprintf('Phi 范围: %.2f 到 %.2f 度\n', min(phi), max(phi));

% 添加静态分析图
figure('Position', [100, 100, 1200, 400]);

% 子图1：速度曲线
subplot(1, 3, 1);
plot(theta_deg, v_E, 'b-', 'LineWidth', 2);
xlabel('Theta (degrees)', 'FontSize', 12);
ylabel('Velocity v_E (mm/s)', 'FontSize', 12);
title('速度 v_E vs. Theta', 'FontSize', 14);
grid on;
xlim([0, 360]);
% 标记极值点
[max_vE, max_idx] = max(v_E);
[min_vE, min_idx] = min(v_E);
hold on;
plot(theta_deg(max_idx), max_vE, 'ro', 'MarkerSize', 10, 'LineWidth', 2);
plot(theta_deg(min_idx), min_vE, 'go', 'MarkerSize', 10, 'LineWidth', 2);
legend('v_E', sprintf('Max: %.1f', max_vE), sprintf('Min: %.1f', min_vE), 'Location', 'best');

% 子图2：加速度曲线
subplot(1, 3, 2);
plot(theta_deg, a_E, 'r-', 'LineWidth', 2);
xlabel('Theta (degrees)', 'FontSize', 12);
ylabel('Acceleration a_E (mm/s²)', 'FontSize', 12);
title('加速度 a_E vs. Theta', 'FontSize', 14);
grid on;
xlim([0, 360]);
% 标记极值点
[max_aE, max_a_idx] = max(a_E);
[min_aE, min_a_idx] = min(a_E);
hold on;
plot(theta_deg(max_a_idx), max_aE, 'mo', 'MarkerSize', 10, 'LineWidth', 2);
plot(theta_deg(min_a_idx), min_aE, 'co', 'MarkerSize', 10, 'LineWidth', 2);
legend('a_E', sprintf('Max: %.1f', max_aE), sprintf('Min: %.1f', min_aE), 'Location', 'best');

% 子图3：phi角度曲线
subplot(1, 3, 3);
plot(theta_deg, phi, 'k-', 'LineWidth', 2);
xlabel('Theta (degrees)', 'FontSize', 12);
ylabel('Phi (degrees)', 'FontSize', 12);
title('Phi vs. Theta', 'FontSize', 14);
grid on;
xlim([0, 360]);

disp('动画完成！');