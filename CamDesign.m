% 参数设置
r0 = 30; % 基圆半径 (mm)
h = 15; % 最大升程 (mm)
theta_deg = 0:0.1:360; % 高精度转角采样（0.1°步长）
theta_rad = deg2rad(theta_deg);
Phi = [270, 40, 10, 40]; % [近休止°,推程°,远休止°,回程°]
n = length(theta_deg);
% 初始化位移数组
s = zeros(size(theta_deg));
ds_dtheta = zeros(1, length(theta_deg));
% 计算位移s(θ)（摆线运动规律）
for i = 1:length(theta_deg)
theta = theta_deg(i);
if theta < 270 % 近休止
s(i) = 0;
elseif theta < 310 % 推程阶段
t = (theta - 270) / 40;
s(i) = h * (t - (1/(2*pi)) * sin(2*pi*t));
elseif theta < 320 % 远休止
s(i) = h;
else % 回程阶段
t = (theta - 320) / 40;
s(i) = h * (1 - t + (1/(2*pi)) * sin(2*pi*t));
end
end
%% 绘制位移曲线
figure;
hold on; grid on; box on;
% 分阶段绘制曲线
idx_near = theta_deg < 270;
plot(theta_deg(idx_near), s(idx_near), 'Color', [0.5,0.5,0.5], 'LineWidth', 2);
idx_rise = (theta_deg >= 270) & (theta_deg < 310);
plot(theta_deg(idx_rise), s(idx_rise), 'Color', [0,0.7,0], 'LineWidth', 2);
idx_upper = (theta_deg >= 310) & (theta_deg < 320);
plot(theta_deg(idx_upper), s(idx_upper), 'Color', [0,0,1], 'LineWidth', 2);
idx_return = theta_deg >= 320;
plot(theta_deg(idx_return), s(idx_return), 'Color', [1,0.5,0], 'LineWidth', 2);
% 标注
xline(270, '--', 'Color', [0.3,0.3,0.3], 'LineWidth', 1);
xline(310, '--', 'Color', [0.3,0.3,0.3], 'LineWidth', 1);
xline(320, '--', 'Color', [0.3,0.3,0.3], 'LineWidth', 1);
text(135, 1, '近休止区 (270°)', 'HorizontalAlignment','center', 'Color', [0.5,0.5,0.5])
text(290, 8, '推程阶段', 'Rotation',20, 'Color', [0,0.7,0])
text(315, 15, '远休止', 'HorizontalAlignment','center', 'Color', [0,0,1])
text(340, 8, '回程阶段', 'Rotation',-20, 'Color', [1,0.5,0])
plot(310, 15, 'ro', 'MarkerFaceColor','r') % 最大升程点
% 图形修饰
title('平底推杆位移曲线 s(θ)');
xlabel('凸轮转角θ (°)'); 
ylabel('推杆位移s (mm)');
xlim([0,360]); 
ylim([-1,16]);
set(gca, 'XTick', 0:30:360);
set(gca, 'FontSize', 10);
legend('近休止', '推程', '远休止', '回程', 'Location', 'northwest');
%% 包络线计算
cam_points = zeros(n, 2);
for i = 1:n
phi = theta_rad(i);
R = r0 + s(i);
dR = ds_dtheta(i); % 已转换为角度导数
% 构造方程组系数矩阵
A = [cos(phi), sin(phi);
-sin(phi), cos(phi)];
b = [R; 
dR]; % 注意此处导数已包含dθ转换
% 解线性方程组
point = A\b;
cam_points(i,:) = point';
end
%% 后处理优化
% 去除重复点（容差1e-4）
[~, unique_idx] = uniquetol(cam_points, 1e-4, 'ByRows', true);
cam_points = cam_points(sort(unique_idx), :);

figure;
hold on; axis equal; grid on;
% 1. 绘制基圆
theta_base = linspace(0, 2*pi, 200);
plot(r0*cos(theta_base), r0*sin(theta_base), 'k--', 'LineWidth', 1.2);
% 2. 分阶段绘制凸轮廓线
% 获取处理后的角度数组
theta_processed = theta_deg(sort(unique_idx));
% 创建各阶段逻辑索引
idx_near = theta_processed < 270; % 近休止
idx_rise = (theta_processed >= 270) & (theta_processed < 310); % 推程
idx_upper = (theta_processed >= 310) & (theta_processed < 320); % 远休止
idx_return = theta_processed >= 320; % 回程
plot(cam_points(idx_near,1), cam_points(idx_near,2), 'Color',[0.5 0.5 0.5], 'LineWidth',2); % 灰色-近休止
plot(cam_points(idx_rise,1), cam_points(idx_rise,2), 'Color',[0 0.7 0], 'LineWidth',2); % 绿色-推程
plot(cam_points(idx_upper,1), cam_points(idx_upper,2), 'Color',[0 0 1], 'LineWidth',2); % 蓝色-远休止
plot(cam_points(idx_return,1), cam_points(idx_return,2), 'Color',[1 0.5 0], 'LineWidth',2);% 橙色-回程
% 3. 标注关键角度位置
marker_size = 12;
% 推程起点(270°)
theta_mark = 270;
[~, idx] = min(abs(theta_processed - theta_mark));
plot(cam_points(idx,1), cam_points(idx,2), '^', 'Color','k',...
'MarkerFaceColor','g', 'MarkerSize',marker_size);
text(cam_points(idx,1)+3, cam_points(idx,2), '推程起点',...
'FontSize',9, 'VerticalAlignment','bottom');
% 远休止起点(310°)
theta_mark = 310;
[~, idx] = min(abs(theta_processed - theta_mark));
plot(cam_points(idx,1), cam_points(idx,2), 's', 'Color','k',...
'MarkerFaceColor','b', 'MarkerSize',marker_size);
text(cam_points(idx,1), cam_points(idx,2)+3, '远休起点',...
'FontSize',9, 'HorizontalAlignment','center');
% 回程起点(320°)
theta_mark = 320;
[~, idx] = min(abs(theta_processed - theta_mark));
plot(cam_points(idx,1), cam_points(idx,2), 'v', 'Color','k',...
'MarkerFaceColor',[1 0.5 0], 'MarkerSize',marker_size);
text(cam_points(idx,1)-5, cam_points(idx,2), '回程起点',...
'FontSize',9, 'HorizontalAlignment','right');
% 图形修饰
title('平底推杆凸轮轮廓（分阶段标注）');
xlabel('x (mm)'); ylabel('y (mm)');
legend('基圆','近休止','推程','远休止','回程','Location','best');
xlim([-50,50]); ylim([-50,50]);