%%% 凸轮编程仿真程序 牛头刨床凸轮
clear all; clc; close all;

% initial
% 滚子摆动盘型凸轮 旋转方向为顺时针
L = 132;                % 从动件杆长(mm)
w_A_rpm = 56;           % 凸轮转速(rpm)
theta_K = 34.84;        % 刨头行程角(度)
phi_max = 15;           % 最大摆角(度)
phi_allowed_go = 40;    % 推程许用压力角(度)
phi_allowed_back = 50;  % 回程许用压力角(度)
r_roller = 15;          % 滚子半径(mm)
r_tool = 8;             % 刀具半径(mm)

% 设计参数
a = 180;                % 中心距(mm)
r0 = 80;                % 基圆半径(mm)

% 凸轮角度分配（根据初始问题描述）
delta_t = 72.58;        % 推程运动角(度)
delta_h = 72.58;        % 回程运动角(度)
delta_s = 360 - delta_t - delta_h;  % 停歇角(度)

fprintf('======================================\n');
fprintf('牛头刨床间歇送料凸轮机构设计\n');
fprintf('凸轮转速: %.0f rpm\n', w_A_rpm);
fprintf('推程角: %.2f°\n', delta_t);
fprintf('回程角: %.2f°\n', delta_h);
fprintf('停歇角: %.2f°\n', delta_s);
fprintf('最大摆幅: %.1f°\n', phi_max);
fprintf('======================================\n');

% 1. 计算初始角ψ₀
psi0 = acos((a^2 + L^2 - r0^2) / (2*a*L));
psi0_deg = rad2deg(psi0);
fprintf('初始角 ψ₀ = %.2f°\n', psi0_deg);

% 2. 定义运动规律（摆线运动）

% 凸轮转角数组（度）
theta_deg = linspace(0, 360, 1001);  % 0到360度
theta_rad = deg2rad(theta_deg);     % 转换为弧度

% 凸轮角度分配（弧度）
theta_rise_start = 0;                    % 推程起始角
theta_rise_end = deg2rad(delta_t);       % 推程结束角
theta_fall_end = deg2rad(delta_t + delta_h);  % 回程结束角
theta_dwell_end = 2*pi;                  % 停歇结束角

fprintf('推程阶段: %.2f° 到 %.2f°\n', rad2deg(theta_rise_start), rad2deg(theta_rise_end));
fprintf('回程阶段: %.2f° 到 %.2f°\n', rad2deg(theta_rise_end), rad2deg(theta_fall_end));
fprintf('停歇阶段: %.2f° 到 %.2f°\n', rad2deg(theta_fall_end), rad2deg(theta_dwell_end));

% 最大摆角（弧度）
Phi_max_rad = deg2rad(phi_max);  

% 推程和回程的角度范围（弧度）
rise_norm = theta_rise_end - theta_rise_start;  % 推程总角
fall_norm = theta_fall_end - theta_rise_end;    % 回程总角

% 初始化运动参数数组（弧度）
psi_rad = zeros(size(theta_rad));      % 摆角（弧度）
dpsi_dtheta = zeros(size(theta_rad));  % 角速度（弧度/弧度）

% 使用循环计算摆角
for i = 1:length(theta_rad)
    current_theta = theta_rad(i);
    
    if current_theta <= theta_rise_start
        % 推程前休止
        psi_rad(i) = 0;
        dpsi_dtheta(i) = 0;
        
    elseif (current_theta > theta_rise_start) && (current_theta <= theta_rise_end)
        % 推程阶段
        theta_rel = current_theta - theta_rise_start;  % 相对转角
        xi = theta_rel / rise_norm;  % 归一化转角
        
        % 摆线运动公式
        psi_rad(i) = Phi_max_rad * (xi - sin(2*pi*xi)/(2*pi));
        dpsi_dtheta(i) = (Phi_max_rad / rise_norm) * (1 - cos(2*pi*xi));
        
    elseif (current_theta > theta_rise_end) && (current_theta <= theta_fall_end)
        % 回程阶段
        theta_rel = current_theta - theta_rise_end;  % 相对转角
        xi = theta_rel / fall_norm;  % 归一化转角
        
        % 摆线运动公式
        psi_rad(i) = Phi_max_rad * (1 - xi + sin(2*pi*xi)/(2*pi));
        dpsi_dtheta(i) = (Phi_max_rad / fall_norm) * (cos(2*pi*xi) - 1);
        
    else
        % 回程后休止
        psi_rad(i) = 0;
        dpsi_dtheta(i) = 0;
    end
end

% 转换为度（用于显示）
psi_deg = rad2deg(psi_rad);

% 3. 压力角计算
% 压力角公式: tan(alpha) = |(L * dpsi/dtheta - a * cos(psi0+psi)) / (a * sin(psi0+psi))|
pressure_angle = zeros(size(theta_rad));

for i = 1:length(theta_rad)
    current_psi = psi_rad(i);
    current_dpsi = dpsi_dtheta(i);
    
    % 计算角度参数
    angle_sum = psi0 + current_psi;
    
    % 分子和分母
    numerator = L * (1 + current_dpsi) - a * cos(angle_sum); % 逆摆式凸轮 +
    denominator = a * sin(angle_sum);
    
    % 避免除零错误
    if abs(denominator) < 1e-6
        if numerator > 0
            pressure_angle(i) = pi/2;  % 正无穷，接近90度
        else
            pressure_angle(i) = -pi/2; % 负无穷，接近-90度
        end
    else
        tan_alpha = numerator / denominator;
        pressure_angle(i) = atan(tan_alpha);
    end
    
end

pressure_angle_deg = rad2deg(abs(pressure_angle));

% 找出推程和回程的最大压力角
rise_indices = (theta_rad >= theta_rise_start) & (theta_rad <= theta_rise_end);
fall_indices = (theta_rad > theta_rise_end) & (theta_rad <= theta_fall_end);

max_pressure_rise = max(pressure_angle_deg(rise_indices));
max_pressure_fall = max(pressure_angle_deg(fall_indices));

fprintf('推程最大压力角: %.2f° (许用: %.2f°)\n', max_pressure_rise, phi_allowed_go);
fprintf('回程最大压力角: %.2f° (许用: %.2f°)\n', max_pressure_fall, phi_allowed_back);

if max_pressure_rise <= phi_allowed_go && max_pressure_fall <= phi_allowed_back
    fprintf('✓ 压力角校核通过！\n');
else
    fprintf('✗ 压力角校核不通过！\n');
end

% 4. 绘制运动规律曲线
figure(1);
clf;
set(gcf, 'Position', [100, 100, 900, 600]);

subplot(3,1,1);
plot(theta_deg, psi_deg, 'b-', 'LineWidth', 1.5);
xlabel('凸轮转角 (度)');
ylabel('摆角 (度)');
title('摆角曲线');
grid on;
xlim([0, 360]);
% 添加区域标记
hold on;
fill([0, delta_t, delta_t, 0], [0, 0, phi_max, phi_max], 'r', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
fill([delta_t, delta_t+delta_h, delta_t+delta_h, delta_t], [phi_max, phi_max, 0, 0], 'g', 'FaceAlpha', 0.1, 'EdgeColor', 'none');

subplot(3,1,2);
plot(theta_deg, dpsi_dtheta, 'r-', 'LineWidth', 1.5);
xlabel('凸轮转角 (度)');
ylabel('dψ/dθ (无量纲)');
title('角速度曲线');
grid on;
xlim([0, 360]);
% 添加区域标记
hold on;
fill([0, delta_t, delta_t, 0], [min(dpsi_dtheta)*1.1, min(dpsi_dtheta)*1.1, max(dpsi_dtheta)*1.1, max(dpsi_dtheta)*1.1], 'r', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
fill([delta_t, delta_t+delta_h, delta_t+delta_h, delta_t], [min(dpsi_dtheta)*1.1, min(dpsi_dtheta)*1.1, max(dpsi_dtheta)*1.1, max(dpsi_dtheta)*1.1], 'g', 'FaceAlpha', 0.1, 'EdgeColor', 'none');

subplot(3,1,3);
plot(theta_deg, pressure_angle_deg, 'g-', 'LineWidth', 1.5);
hold on;

% 绘制许用压力角线
plot([0, 360], [phi_allowed_go, phi_allowed_go], 'r--', 'LineWidth', 1);
plot([0, 360], [phi_allowed_back, phi_allowed_back], 'm--', 'LineWidth', 1);

% 分别找出推程和回程的最大压力角及其位置
% 推程最大压力角
[max_rise_val, max_rise_idx] = max(pressure_angle_deg(rise_indices));
% 注意：max_rise_idx是rise_indices中为true的索引，需要转换为原始theta_deg的索引
rise_true_indices = find(rise_indices);
rise_max_idx = rise_true_indices(max_rise_idx);
rise_max_theta = theta_deg(rise_max_idx);

% 回程最大压力角
[max_fall_val, max_fall_idx] = max(pressure_angle_deg(fall_indices));
fall_true_indices = find(fall_indices);
fall_max_idx = fall_true_indices(max_fall_idx);
fall_max_theta = theta_deg(fall_max_idx);

% 标记推程最大压力角（红色圆）
plot(rise_max_theta, max_rise_val, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
text(rise_max_theta, max_rise_val+2, sprintf('推程最大: %.2f°', max_rise_val), ...
    'HorizontalAlignment', 'center', 'FontSize', 9);

% 标记回程最大压力角（蓝色方框）
plot(fall_max_theta, max_fall_val, 'bs', 'MarkerSize', 8, 'MarkerFaceColor', 'b');
text(fall_max_theta, max_fall_val+2, sprintf('回程最大: %.2f°', max_fall_val), ...
    'HorizontalAlignment', 'center', 'FontSize', 9);

% 也可以标记整体最大压力角
[max_val, max_idx] = max(pressure_angle_deg);
if max_val == max_rise_val
    % 整体最大在推程
    text(theta_deg(max_idx), max_val-4, sprintf('整体最大: %.2f°', max_val), ...
        'HorizontalAlignment', 'center', 'FontSize', 9, 'Color', 'r');
else
    % 整体最大在回程
    text(theta_deg(max_idx), max_val-4, sprintf('整体最大: %.2f°', max_val), ...
        'HorizontalAlignment', 'center', 'FontSize', 9, 'Color', 'b');
end

xlabel('凸轮转角 (度)');
ylabel('压力角 (度)');
title('压力角曲线');
grid on;
xlim([0, 360]);
legend('压力角', '推程许用值', '回程许用值', '推程最大值', '回程最大值', 'Location', 'best');



% 5. 凸轮轮廓绘制
fprintf('\n====== 开始绘制凸轮轮廓 ======\n');

% 反转法绘制凸轮轮廓
% 对于每个凸轮转角θ，将整个机构以-θ角速度绕凸轮轴心旋转（使凸轮静止）

% 初始化理论轮廓坐标数组
x_theory = zeros(size(theta_rad));  % 理论轮廓X坐标
y_theory = zeros(size(theta_rad));  % 理论轮廓Y坐标

% 初始化实际轮廓坐标数组
x_actual = zeros(size(theta_rad));  % 实际轮廓X坐标
y_actual = zeros(size(theta_rad));  % 实际轮廓Y坐标

% 计算理论轮廓线（滚子中心轨迹）
for i = 1:length(theta_rad)
    x_theory(i) = (a - L * cos(psi0 + psi_rad(i))) * cos(theta_rad(i)) + ...
        L * sin(psi0 + psi_rad(i)) * sin(theta_rad(i)) ;
    y_theory(i) = -(a - L * cos(psi0 + psi_rad(i))) * sin(theta_rad(i)) + ...
        L * sin(psi0 + psi_rad(i)) * cos(theta_rad(i)) ;
end
% 计算理论轮廓线的导数（用于求法线方向）
dx_dtheta = zeros(size(theta_rad));
dy_dtheta = zeros(size(theta_rad));

for i = 1:length(theta_rad)
    current_theta = theta_rad(i);
    current_psi = psi_rad(i);
    current_dpsi = dpsi_dtheta(i);
    
    % 导数计算（对理论轮廓公式求导）
    % 令: u = a - L*cos(ψ0+ψ), v = L*sin(ψ0+ψ)
    % 则: x = u*cosθ + v*sinθ, y = -u*sinθ + v*cosθ
    % 对θ求导:
    % dx/dθ = du/dθ*cosθ - u*sinθ + dv/dθ*sinθ + v*cosθ
    % dy/dθ = -du/dθ*sinθ - u*cosθ + dv/dθ*cosθ - v*sinθ
    % 其中:
    % du/dθ = L*sin(ψ0+ψ)*dψ/dθ
    % dv/dθ = L*cos(ψ0+ψ)*dψ/dθ
    
    sin_theta = sin(current_theta);
    cos_theta = cos(current_theta);
    sin_psi_sum = sin(psi0 + current_psi);
    cos_psi_sum = cos(psi0 + current_psi);
    
    u = a - L * cos_psi_sum;
    v = L * sin_psi_sum;
    du_dtheta = L * sin_psi_sum * current_dpsi;
    dv_dtheta = L * cos_psi_sum * current_dpsi;
    
    dx_dtheta(i) = du_dtheta * cos_theta - u * sin_theta + dv_dtheta * sin_theta + v * cos_theta;
    dy_dtheta(i) = -du_dtheta * sin_theta - u * cos_theta + dv_dtheta * cos_theta - v * sin_theta;
end

for i = 1:length(theta_rad)
    % 计算切向量模长
    v_norm = sqrt(dx_dtheta(i)^2 + dy_dtheta(i)^2);
    
    % 避免除零
    if v_norm < 1e-6
        % 特殊点处理
        if i == 1
            % 使用前向差分
            Nx = (y_theory(i+1) - y_theory(i)) / norm([x_theory(i+1)-x_theory(i), y_theory(i+1)-y_theory(i)]);
            Ny = -(x_theory(i+1) - x_theory(i)) / norm([x_theory(i+1)-x_theory(i), y_theory(i+1)-y_theory(i)]);
        elseif i == length(theta_rad)
            % 使用后向差分
            Nx = (y_theory(i) - y_theory(i-1)) / norm([x_theory(i)-x_theory(i-1), y_theory(i)-y_theory(i-1)]);
            Ny = -(x_theory(i) - x_theory(i-1)) / norm([x_theory(i)-x_theory(i-1), y_theory(i)-y_theory(i-1)]);
        else
            % 使用中心差分
            Nx = (y_theory(i+1) - y_theory(i-1)) / norm([x_theory(i+1)-x_theory(i-1), y_theory(i+1)-y_theory(i-1)]);
            Ny = -(x_theory(i+1) - x_theory(i-1)) / norm([x_theory(i+1)-x_theory(i-1), y_theory(i+1)-y_theory(i-1)]);
        end
    else
        % 内法线单位向量
        Nx = dy_dtheta(i) / v_norm;
        Ny = -dx_dtheta(i) / v_norm;
    end
    
    % 确保内法线指向凸轮内部（朝向凸轮轴心）
    % 计算从理论轮廓点到凸轮轴心的向量
    vec_to_center = [-x_theory(i), -y_theory(i)];
    if norm(vec_to_center) > 1e-6
        vec_to_center = vec_to_center / norm(vec_to_center);
        dot_product = Nx * vec_to_center(1) + Ny * vec_to_center(2);
        
        % 如果点积为负，说明法线方向指向外部，需要取反
        if dot_product < 0
            Nx = -Nx;
            Ny = -Ny;
        end
    end
    
    % 计算实际轮廓点
    x_actual(i) = x_theory(i) + r_roller * Nx;
    y_actual(i) = y_theory(i) + r_roller * Ny;
end

% 计算基圆
theta_base = linspace(0, 2*pi, 361);
x_base = r0 * cos(theta_base);
y_base = r0 * sin(theta_base);

% 验证公式正确性（检查初始位置）
fprintf('\n====== 验证公式正确性 ======\n');
fprintf('初始位置 (θ=0°, ψ=0°)：\n');

% 计算初始位置的理论轮廓点
theta_initial = 0;
psi_initial = 0;
cos_theta = cos(theta_initial);
sin_theta = sin(theta_initial);
cos_psi_sum = cos(psi0 + psi_initial);
sin_psi_sum = sin(psi0 + psi_initial);

x_initial = (a - L * cos_psi_sum) * cos_theta + L * sin_psi_sum * sin_theta;
y_initial = -(a - L * cos_psi_sum) * sin_theta + L * sin_psi_sum * cos_theta;

% 计算距离凸轮轴心的距离
r_initial = sqrt(x_initial^2 + y_initial^2);

fprintf('理论计算得到的初始半径: %.4f mm\n', r_initial);
fprintf('基圆半径 r0: %.4f mm\n', r0);
fprintf('误差: %.6f mm\n', abs(r_initial - r0));

if abs(r_initial - r0) < 1e-3
    fprintf('✓ 公式正确：初始半径等于基圆半径\n');
else
    fprintf('⚠ 公式可能有误：初始半径不等于基圆半径\n');
end


% 6. 绘制凸轮轮廓图
figure(2);
clf;
set(gcf, 'Position', [100, 50, 1200, 800]);

% 子图1：理论轮廓和实际轮廓对比
subplot(2, 2, 1);
plot(x_theory, y_theory, 'b-', 'LineWidth', 1.5, 'DisplayName', '理论轮廓');
hold on;
plot(x_actual, y_actual, 'r-', 'LineWidth', 1.5, 'DisplayName', '实际轮廓');
plot(x_base, y_base, 'g--', 'LineWidth', 1, 'DisplayName', '基圆');
plot(0, 0, 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'k', 'DisplayName', '凸轮轴心');
axis equal;
grid on;
xlabel('X (mm)');
ylabel('Y (mm)');
title('凸轮轮廓对比');
legend('Location', 'best');
xlim([-a-50, a+50]);
ylim([-a-50, a+50]);

% 子图4：轮廓半径变化
subplot(2, 2, 2);
% 计算轮廓点到凸轮轴心的距离
radius_theory = sqrt(x_theory.^2 + y_theory.^2);
radius_actual = sqrt(x_actual.^2 + y_actual.^2);

plot(theta_deg, radius_theory, 'b-', 'LineWidth', 1.5, 'DisplayName', '理论轮廓半径');
hold on;
plot(theta_deg, radius_actual, 'r-', 'LineWidth', 1.5, 'DisplayName', '实际轮廓半径');
plot([0, 360], [r0, r0], 'g--', 'LineWidth', 1.5, 'DisplayName', '基圆半径');

% 标记最小值
[~, idx_min_actual] = min(radius_actual);
[~, idx_min_theory] = min(radius_theory);
plot(theta_deg(idx_min_actual), min(radius_actual), 'ro', 'MarkerSize', 8, ...
    'MarkerFaceColor', 'r', 'DisplayName', '实际最小');
plot(theta_deg(idx_min_theory), min(radius_theory), 'bo', 'MarkerSize', 8, ...
    'MarkerFaceColor', 'b', 'DisplayName', '理论最小');

xlabel('凸轮转角 (度)');
ylabel('半径 (mm)');
title('轮廓半径变化');
grid on;
legend('Location', 'best');
xlim([0, 360]);

% 子图5：轮廓闭合检查
subplot(2, 2, 3);
% 检查轮廓是否闭合
gap_start_end = norm([x_actual(1)-x_actual(end), y_actual(1)-y_actual(end)]);
fprintf('\n轮廓闭合检查：\n');
fprintf('起始点和终止点的距离: %.6f mm\n', gap_start_end);

if gap_start_end < 1e-3
    fprintf('✓ 轮廓完全闭合\n');
elseif gap_start_end < 0.1
    fprintf('✓ 轮廓基本闭合（微小间隙）\n');
else
    fprintf('⚠ 轮廓未闭合，有%.3f mm间隙\n', gap_start_end);
end

% 绘制起始点和终止点
plot(x_actual, y_actual, 'b-', 'LineWidth', 1.5);
hold on;
plot(x_actual(1), y_actual(1), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r', ...
    'DisplayName', '起始点');
plot(x_actual(end), y_actual(end), 'gs', 'MarkerSize', 10, 'MarkerFaceColor', 'g', ...
    'DisplayName', '终止点');

if gap_start_end > 1e-3
    % 绘制连接线
    plot([x_actual(1), x_actual(end)], [y_actual(1), y_actual(end)], 'r--', ...
        'LineWidth', 1, 'DisplayName', '间隙');
    text((x_actual(1)+x_actual(end))/2, (y_actual(1)+y_actual(end))/2, ...
        sprintf('%.3f mm', gap_start_end), 'HorizontalAlignment', 'center');
end

axis equal;
grid on;
xlabel('X (mm)');
ylabel('Y (mm)');
title('轮廓闭合检查');
legend('Location', 'best');

% 子图6：3D凸轮视图
subplot(2, 2, 4);
% 创建极坐标网格
[TH, R] = meshgrid(linspace(0, 2*pi, 100), linspace(min(radius_actual), max(radius_actual), 50));
% 将极坐标转换为直角坐标
[X, Y] = pol2cart(TH, R);

% 创建Z坐标（凸轮厚度）
Z = zeros(size(X));
thickness = 20;  % 凸轮厚度(mm)
Z = Z + thickness/2;  % 凸轮上表面

% 绘制3D凸轮
surf(X, Y, Z, 'FaceColor', [0.8, 0.8, 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.8);
hold on;

% 绘制基圆柱
[X_base, Y_base, Z_base] = cylinder(r0, 50);
Z_base = Z_base * thickness;
surf(X_base, Y_base, Z_base, 'FaceColor', [0.6, 0.9, 0.6], 'EdgeColor', 'none', 'FaceAlpha', 0.5);

% 绘制实际轮廓线（作为3D凸轮的边缘）
plot3(x_actual, y_actual, ones(size(x_actual))*thickness, 'r-', 'LineWidth', 2);
plot3(x_actual, y_actual, zeros(size(x_actual)), 'r-', 'LineWidth', 2);

axis equal;
grid on;
xlabel('X (mm)');
ylabel('Y (mm)');
zlabel('Z (mm)');
title('3D凸轮视图');
view(30, 30);
light('Position', [1, 1, 1], 'Style', 'infinite');
lighting gouraud;

%% 7. 输出轮廓数据
fprintf('\n====== 轮廓数据统计 ======\n');
fprintf('理论轮廓最大半径: %.2f mm\n', max(radius_theory));
fprintf('实际轮廓最大半径: %.2f mm\n', max(radius_actual));
fprintf('理论轮廓最小半径: %.2f mm\n', min(radius_theory));
fprintf('实际轮廓最小半径: %.2f mm\n', min(radius_actual));
fprintf('基圆半径: %.2f mm\n', r0);
fprintf('滚子半径: %.2f mm\n', r_roller);

% 检查实际轮廓最小半径是否接近基圆半径
if abs(min(radius_theory) - r0) < 1
    fprintf('✓ 理论轮廓最小半径接近基圆半径，轮廓正确\n');
else
    fprintf('⚠ 注意：理论轮廓最小半径与基圆半径有差异\n');
    fprintf('   差异: %.2f mm\n', min(radius_actual) - r0);
end

% 计算最小曲率半径（检查过切风险）
fprintf('\n====== 过切风险检查 ======\n');
curvature_radii = zeros(length(theta_rad)-2, 1);
for i = 2:length(theta_rad)-1
    x1 = x_actual(i-1); y1 = y_actual(i-1);
    x2 = x_actual(i);   y2 = y_actual(i);
    x3 = x_actual(i+1); y3 = y_actual(i+1);
    
    % 计算三角形面积
    area = abs((x2-x1)*(y3-y1) - (x3-x1)*(y2-y1)) / 2;
    
    % 计算边长
    a_len = sqrt((x2-x1)^2 + (y2-y1)^2);
    b_len = sqrt((x3-x2)^2 + (y3-y2)^2);
    c_len = sqrt((x1-x3)^2 + (y1-y3)^2);
    
    if area < 1e-6
        curvature_radii(i-1) = Inf;
    else
        curvature_radii(i-1) = (a_len * b_len * c_len) / (4 * area);
    end
end

min_curvature = min(curvature_radii);
fprintf('最小曲率半径: %.2f mm\n', min_curvature);

if min_curvature > r_roller
    fprintf('✓ 曲率半径大于滚子半径，无过切风险\n');
else
    fprintf('⚠ 警告：最小曲率半径小于滚子半径，有过切风险！\n');
end

% 保存轮廓数据
fprintf('\n保存轮廓数据...\n');
theta_dense = linspace(0, 2*pi, 3601);  % 每0.1度一个点
x_actual_dense = interp1(theta_rad, x_actual, theta_dense, 'spline');
y_actual_dense = interp1(theta_rad, y_actual, theta_dense, 'spline');

% 保存为CSV格式
T = table(theta_dense', x_actual_dense', y_actual_dense', ...
    'VariableNames', {'Theta_rad', 'X_mm', 'Y_mm'});
writetable(T, 'cam_profile_corrected.csv');
fprintf('轮廓数据已保存到 cam_profile_corrected.csv\n');

fprintf('\n======================================\n');
fprintf('凸轮轮廓绘制完成（使用正确公式）！\n');
fprintf('======================================\n');


