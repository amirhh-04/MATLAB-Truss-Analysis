clc; clear; close all;
data = cell(1, 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%   Inputs Section    %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [x y z]
data{1, 1}.nodes = [1.5 0 2
                    4.5 0 2
                    0 0 0
                    3 0 0
                    6 0 0
                    1.5 3 2
                    4.5 3 2
                    0 3 0
                    3 3 0
                    6 3 0];
nodeCount = size(data{1, 1}.nodes, 1);

% [node1 node2 A(m²) E(Pa) yieldStress(Pa)]
data{1, 1}.elements = [1 2 0.0005 200e9 20e7
                       1 3 0.0005 200e9 20e7
                       1 4 0.0005 200e9 20e7
                       2 4 0.0005 200e9 20e7
                       2 5 0.0005 200e9 20e7
                       3 4 0.0005 200e9 20e7
                       4 5 0.0005 200e9 20e7
                       6 7 0.0005 200e9 20e7
                       6 8 0.0005 200e9 20e7
                       6 9 0.0005 200e9 20e7
                       7 9 0.0005 200e9 20e7
                       7 10 0.0005 200e9 20e7
                       8 9 0.0005 200e9 20e7
                       9 10 0.0005 200e9 20e7
                       3 8 0.0005 200e9 20e7
                       1 6 0.0005 200e9 20e7
                       4 9 0.0005 200e9 20e7
                       2 7 0.0005 200e9 20e7
                       5 10 0.0005 200e9 20e7
                       4 8 0.0005 200e9 20e7
                       4 10 0.0005 200e9 20e7
                       1 7 0.0005 200e9 20e7];
elementCount = size(data{1, 1}.elements, 1);

% [nodeNumber F(N)-x F(N)-y F(N)-z]
data{1, 1}.forces = [1 0 0 -25000
                     2 0 0 -25000
                     6 0 0 -25000
                     7 0 0 -25000];
forceNodesCount = size(data{1, 1}.forces, 1);

% [nodeNumber Dx Dy Dz] (D => 0[free],1[close])
data{1, 1}.supports = [3 1 1 1
                       5 1 1 1
                       8 1 1 1
                       10 1 1 1];
supportCount = size(data{1, 1}.supports, 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%  Main Code Section  %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data{1, 1}.element_lengths = zeros(elementCount, 1);
data{1, 1}.element_angles = zeros(elementCount, 3);

k_data = zeros(3 * nodeCount);
k0_data = k_data;
for i = 1:elementCount
    [n_1, n_2] = getNodesByElement(i, data);
    [element_lengths, element_angles_x, element_angles_y, element_angles_z, area, young_modulus, delta_f_mat] = getElementInformation(i, n_1, n_2, data);

    data{1, 1}.element_lengths(i) = element_lengths;

    data{1, 1}.element_angles(i, 1) = element_angles_x;
    data{1, 1}.element_angles(i, 2) = element_angles_y;
    data{1, 1}.element_angles(i, 3) = element_angles_z;

    k{i} = area * young_modulus / element_lengths * [delta_f_mat -delta_f_mat; -delta_f_mat delta_f_mat];

    [n_1, n_2] = getNodesByElement(i, data);

    k0_data(3 * n_1 - 2:3 * n_1, 3 * n_1 - 2:3 * n_1) = k{i}(1:3, 1:3);
    k0_data(3 * n_1 - 2:3 * n_1, 3 * n_2 - 2:3 * n_2) = k{i}(1:3, 4:6);
    k0_data(3 * n_2 - 2:3 * n_2, 3 * n_1 - 2:3 * n_1) = k{i}(4:6, 1:3);
    k0_data(3 * n_2 - 2:3 * n_2, 3 * n_2 - 2:3 * n_2) = k{i}(4:6, 4:6);

    k_data = k_data + k0_data;
    k0_data(:, :) = 0;
end

force_0 = zeros(3 * nodeCount, 1);
for i = 1:forceNodesCount
    [force_node, fx, fy, fz] = getForceInfo(i, data);
    force_0((3 * force_node) - 2:(3 * force_node), 1) = [fx; fy; fz];
end

num_data = [];
for i = 1:supportCount
    [support_node, Dx, Dy, Dz, num] = getSupportInfo(i, data, num_data);
    num_data = num(num ~= 0);
    uu_zero = [num_data(:)];
end

k_small_data = k_data;
forces_data = force_0;
k_small_data(:, uu_zero) = [];
k_small_data(uu_zero, :) = [];
forces_data(uu_zero, :) = [];
u_all_data = (k_small_data^-1) * forces_data;
uu_all = 1:(3 * nodeCount);
uu_nonzero = uu_all;
uu_nonzero(uu_zero) = [];
u_data(uu_all, 1) = 0;
u_data(uu_nonzero, 1) = u_all_data;
force = k_data * u_data;

element_res = cell(elementCount, 1);
data{1, 1}.newNodesLoc = data{1, 1}.nodes;
for i = 1:elementCount
    c_o_s_x = cos(data{1, 1}.element_angles(i, 1)*pi/180);
    c_o_s_y = cos(data{1, 1}.element_angles(i, 2)*pi/180);
    c_o_s_z = cos(data{1, 1}.element_angles(i, 3)*pi/180);

    [n_1, n_2] = getNodesByElement(i, data);

    delta = [-c_o_s_x -c_o_s_y -c_o_s_z c_o_s_x c_o_s_y c_o_s_z] * [u_data((3*n_1) - 2); u_data((3*n_1) - 1); u_data(3*n_1); u_data((3*n_2) - 2); u_data((3*n_2) - 1); u_data(3*n_2)];

    n1_new_loc_x = u_data((3 * n_1) - 2);
    n1_new_loc_y = u_data((3 * n_1) - 1);
    n1_new_loc_z = u_data(3 * n_1);

    n2_new_loc_x = u_data((3 * n_2) - 2);
    n2_new_loc_y = u_data((3 * n_2) - 1);
    n2_new_loc_z = u_data(3 * n_2);

    data{1, 1}.newNodesLoc(n_1, :) = data{1, 1}.newNodesLoc(n_1, :) + [n1_new_loc_x, n1_new_loc_y, n1_new_loc_z];
    data{1, 1}.newNodesLoc(n_2, :) = data{1, 1}.newNodesLoc(n_2, :) + [n2_new_loc_x, n2_new_loc_y, n2_new_loc_z];

    area = data{1, 1}.elements(i, 3);
    young_modulus = data{1, 1}.elements(i, 4);

    p = area * young_modulus / data{1, 1}.element_lengths(i) * delta;

    element_res{i, 1}.delta_len = delta;
    element_res{i, 1}.force = p;
    element_res{i, 1}.stress = p/area;

    [delta_len, force, stress, strain, stress_utilisation] = getResultByElement(i, data, element_res);
    element_res{i, 1}.stress_utilisation = stress_utilisation;
    element_res{i, 1}.strain = strain;

    fprintf(' -------------------- Element (%g) [n: %g,%g] --------------------\n   Delta: %g (um) \n   Force: %g (N) \n   Stress: %g (MPa) \n   Strain: %g \n   Stress Utilisation: %g%%  \n', i, n_1, n_2, delta_len, force, stress, strain, stress_utilisation);
end
fprintf('\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%   Figure Section    %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nodes = data{1, 1}.nodes;
new_nodes = data{1, 1}.newNodesLoc;
elements = data{1, 1}.elements;
supports = data{1, 1}.supports;

figure;
scatter3(nodes(:, 1), nodes(:, 2), nodes(:, 3), 120, [0.8500 0.3250 0.0980], 'filled');
hold on;

for i = 1:elementCount
    [node1, node2] = getNodesByElement(i, data);    
    [x1, y1, z1] = getNodePos(node1, false, data);
    [x2, y2, z2] = getNodePos(node2, false, data);
    [new_x1, new_y1, new_z1] = getNodePos(node1, true, data);
    [new_x2, new_y2, new_z2] = getNodePos(node2, true, data);
     
    broken = false;
    if element_res{i, 1}.stress_utilisation > 100; color = [0.6350 0.0780 0.1840]; broken = true; else; color = [0 0.4470 0.7410]; end;
    
    plot3([new_x1, new_x2], [new_y1, new_y2], [new_z1, new_z2], 'Color', "magenta", 'LineWidth', 1);
    plot3([x1, x2], [y1, y2], [z1, z2], 'Color', color, 'LineWidth', 4);

    text(x2, y2, z2+.23, sprintf('%d', node2), 'FontSize', 15, 'Color', "#279501");
    text(x1, y1, z1+.23, sprintf('%d', node1), 'FontSize', 15, 'Color', "#279501");

    mid_x = (x1 + x2) / 2;
    mid_y = (y1 + y2) / 2;
    mid_z = (z1 + z2) / 2;

    text(mid_x, mid_y, mid_z, sprintf('%d', i), 'FontSize', 9, 'Color', "#FC0B8E");
end

for i = 1:supportCount
    [support_node, Dx, Dy, Dz, num] = getSupportInfo(i, data, false);
    [x, y, z] = getNodePos(support_node, false, data);
    
    if Dx == 1 && Dy == 1 && Dz == 1
        color = 'k';
    elseif Dx == 0 && Dy == 0 && Dz == 0
        color = 'b';
    elseif Dx == 1 && Dy == 0 && Dz == 0
        color = 'g';
    elseif Dx == 0 && Dy == 1 && Dz == 0
        color = 'c';
    elseif Dx == 0 && Dy == 0 && Dz == 1
        color = 'm';
    elseif Dx == 1 && Dy == 1 && Dz == 0
        color = 'w';
    elseif Dx == 0 && Dy == 1 && Dz == 1
        color = [0.4660 0.6740 0.1880];
    elseif Dx == 1 && Dy == 0 && Dz == 1
        color = [0.8500 0.3250 0.0980]	;
    end
    scatter3(x, y, z, 120, color, 'filled');
end

title('STRUCTURE: ');
grid on;
xlabel('X-Axis');
ylabel('Y-Axis');
zlabel('Z-Axis');
axis equal;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%  Functions Section  %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [delta_len, force, stress, strain, stress_utilisation] = getResultByElement(i, data, element_res)
    el_data = element_res{i, 1};
    yield_stress = data{1, 1}.elements(i, 5);
    el_len = data{1, 1}.element_lengths(i);

    delta_len = meterTomicroMetres(el_data.delta_len);
    force = el_data.force;
    stress = pascalsToMegapascals(el_data.stress);
    strain = (delta_len) / el_len;
    stress_utilisation = abs(el_data.stress / yield_stress) * 100;
end

function [support_node, Dx, Dy, Dz, num] = getSupportInfo(i, data, num)
    support_node = data{1, 1}.supports(i, 1);
    Dx = data{1, 1}.supports(i, 2);
    Dy = data{1, 1}.supports(i, 3);
    Dz = data{1, 1}.supports(i, 4);
    
    if ~(num); return; end;

    if Dx == 1; num(1, (3 * support_node) - 2) = (3 * support_node) - 2; end;
    if Dy == 1; num(1, (3 * support_node) - 1) = (3 * support_node) - 1; end;
    if Dz == 1; num(1, (3 * support_node)) = (3 * support_node); end;
end

function [force_node, fx, fy, fz] = getForceInfo(i, data)
    force_node = data{1, 1}.forces(i, 1);
    force_len = data{1, 1}.forces(i, 2);

    %force_angle_x = data{1, 1}.forces(i, 3) / 180 * pi;
    %force_angle_y = data{1, 1}.forces(i, 4) / 180 * pi;
    %force_angle_z = data{1, 1}.forces(i, 5) / 180 * pi;

    %fx = force_len * cos(force_angle_x);
    %fy = force_len * cos(force_angle_y);
    %fz = force_len * cos(force_angle_z);

    fx = data{1, 1}.forces(i, 2);
    fy = data{1, 1}.forces(i, 3);
    fz = data{1, 1}.forces(i, 4);
end

function [element_lengths, element_angles_x, element_angles_y, element_angles_z, area, young_modulus, delta_f_mat] = getElementInformation(i, n_1, n_2, data)
    [x_1, y_1, z_1] = getNodePos(n_1, false, data);
    [x_2, y_2, z_2] = getNodePos(n_2, false, data);
    delta_x = x_2 - x_1;
    delta_y = y_2 - y_1;
    delta_z = z_2 - z_1;
    
    element_lengths = sqrt((delta_x^2) + (delta_y^2) + (delta_z^2));
    element_angles_x = acos(delta_x/element_lengths) * (180/pi);
    element_angles_y = acos(delta_y/element_lengths) * (180/pi);
    element_angles_z = acos(delta_z/element_lengths) * (180/pi);
    
    area = data{1, 1}.elements(i, 3);
    young_modulus = data{1, 1}.elements(i, 4);
    
    x = element_angles_x * (pi/180);
    y = element_angles_y * (pi/180);
    z = element_angles_z * (pi/180);
    c_o_s_x = cos(x);
    c_o_s_y = cos(y);
    c_o_s_z = cos(z);
    delta_f_mat = [c_o_s_x*c_o_s_x c_o_s_x*c_o_s_y c_o_s_x*c_o_s_z ; c_o_s_y*c_o_s_x c_o_s_y*c_o_s_y c_o_s_y*c_o_s_z ; c_o_s_z*c_o_s_x c_o_s_z*c_o_s_y c_o_s_z*c_o_s_z];
end

function [x_node, y_node, z_node] = getNodePos(node, newPos, data)
    if newPos == true
        new_nodes = data{1, 1}.newNodesLoc;
        x_node = new_nodes(node, 1);
        y_node = new_nodes(node, 2);
        z_node = new_nodes(node, 3);
        return
    end
    nodes = data{1, 1}.nodes;
    x_node = nodes(node, 1);
    y_node = nodes(node, 2);
    z_node = nodes(node, 3);
end

function [node_1, node_2] = getNodesByElement(element, data)
    node_1 = data{1, 1}.elements(element, 1);
    node_2 = data{1, 1}.elements(element, 2);
end

function microMetres = meterTomicroMetres(meter)
    microMetres = meter * 1e6;
end

function megapascals = pascalsToMegapascals(pascals)
    megapascals = pascals / 1e6;
end