clc; clear; close all;

data = cell(1, 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%   Inputs Section    %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [x y z]
data{1, 1}.nodes = [0 0
                    3.8 9.12/7
                    3.8 0
                    7 2.4
                    7 0
                    10.2 9.12/7
                    10.2 0
                    14 0];
nodeCount = size(data{1, 1}.nodes, 1);

% [node1 node2 A(m²) E(Pa) yieldStress(Pa)]
data{1, 1}.elements = [1 2 0.0005 200e9 20e7
                       1 3 0.0005 200e9 20e7
                       2 3 0.0005 200e9 20e7
                       2 4 0.0005 200e9 20e7
                       3 4 0.0005 200e9 20e7
                       3 5 0.0005 200e9 20e7
                       4 5 0.0005 200e9 20e7
                       4 7 0.0005 200e9 20e7
                       4 6 0.0005 200e9 20e7
                       6 7 0.0005 200e9 20e7
                       6 8 0.0005 200e9 20e7
                       5 7 0.0005 200e9 20e7
                       7 8 0.0005 200e9 20e7];
elementCount = size(data{1, 1}.elements, 1);

% [nodeNumber F(N) nodeAngle]
data{1, 1}.forces = [1 5700 -90
                     2 10500 -90
                     4 9600 -90
                     6 10500 -90
                     8 5700 -90];
forceNodesCount = size(data{1, 1}.forces, 1);

% [nodeNumber type(1,2) orien(0,1,2)]
data{1, 1}.supports = [1 2 0
                       8 1 2];
supportCount = size(data{1, 1}.supports, 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%  Main Code Section  %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data{1, 1}.element_lengths = zeros(elementCount, 1);
data{1, 1}.element_angles = zeros(elementCount, 1);

k_data = zeros(2 * nodeCount);
k0_data = k_data;
for i = 1:elementCount
    [n_1, n_2] = getNodesByElement(i, data);
    [element_lengths, element_angles, area, young_modulus, delta_f_mat] = getElementInformation(i, n_1, n_2, data);

    data{1, 1}.element_lengths(i) = element_lengths;
    data{1, 1}.element_angles(i) = element_angles;
    k{i} = area * young_modulus / element_lengths * delta_f_mat;
    [n_1, n_2] = getNodesByElement(i, data);

    k0_data(2*n_1-1:2*n_1, 2*n_1-1:2*n_1) = k{i} (1:2, 1:2);
    k0_data(2*n_1-1:2*n_1, 2*n_2-1:2*n_2) = k{i} (1:2, 3:4);
    k0_data(2*n_2-1:2*n_2, 2*n_1-1:2*n_1) = k{i} (3:4, 1:2);
    k0_data(2*n_2-1:2*n_2, 2*n_2-1:2*n_2) = k{i} (3:4, 3:4);

    k_data = k_data + k0_data;
    k0_data(:, :) = 0;
end

force_0 = zeros(2*nodeCount, 1);
for i = 1:forceNodesCount
    [force_node, fx, fy] = getForceInfo(i, data);
    force_0((2 * force_node) - 1:(2 * force_node), 1) = [fx; fy];
end

counter = 0;
for i = 1:supportCount
    [support_node, support_type, support_orien] = getSupportInfo(i, data);

    if support_type == 1
        if support_orien == 1
            counter = counter + 1;
            uu_zero(counter) = (2 * support_node) - 1;
        elseif support_orien == 2
            counter = counter + 1;
            uu_zero(counter) = 2 * support_node;
        end
    elseif support_type == 2
        counter = counter + 2;
        uu_zero(counter-1:counter) = (2 * support_node) - 1:(2 * support_node);
    end
end

k_small_data = k_data;
forces_data = force_0;
k_small_data(:, uu_zero) = [];
k_small_data(uu_zero, :) = [];
forces_data(uu_zero, :) = [];
u_all_data = (k_small_data^-1) * forces_data;
uu_all = 1:(2 * nodeCount);
uu_nonzero = uu_all;
uu_nonzero(uu_zero) = [];
u_data(uu_all, 1) = 0;
u_data(uu_nonzero, 1) = u_all_data;
force = k_data * u_data;

element_res = cell(elementCount, 1);
data{1, 1}.newNodesLoc = data{1, 1}.nodes;
for i = 1:elementCount
    s_i_n = sin(data{1, 1}.element_angles(i));
    c_o_s = cos(data{1, 1}.element_angles(i));

    [n_1, n_2] = getNodesByElement(i, data);

    delta = [-c_o_s -s_i_n c_o_s s_i_n] * [u_data(2*n_1-1); u_data(2*n_1); u_data(2*n_2-1); u_data(2*n_2)];

    n1_new_loc_x = u_data((2 * n_1) - 1);
    n1_new_loc_y = u_data(2 * n_1);

    n2_new_loc_x = u_data((2 * n_2) - 1);
    n2_new_loc_y = u_data(2 * n_2);

    data{1, 1}.newNodesLoc(n_1, :) = data{1, 1}.newNodesLoc(n_1, :) + [n1_new_loc_x, n1_new_loc_y];
    data{1, 1}.newNodesLoc(n_2, :) = data{1, 1}.newNodesLoc(n_2, :) + [n2_new_loc_x, n2_new_loc_y];

    area = data{1, 1}.elements(i, 3);
    young_modulus = data{1, 1}.elements(i, 4);

    p = area * young_modulus / data{1, 1}.element_lengths(i) * delta;

    element_res{i, 1}.delta_len = delta;
    element_res{i, 1}.force = p;
    element_res{i, 1}.stress = p/area;

    [delta_len, force, stress, strain, stress_utilisation] = getResultByElement(i, data, element_res);
    element_res{i, 1}.stress_utilisation = stress_utilisation;
    element_res{i, 1}.strain = strain;

    % Static Analys %
    %fprintf('\n -------------------- Element (%g) [n: %g,%g] --------------------\n  Force: %g (N) \n', i, n_1, n_2, force);

    fprintf('\n -------------------- Element (%g) [n: %g,%g] --------------------\n   Delta: %g (um) \n   Force: %g (N) \n   Stress: %g (MPa) \n   Strain: %g \n   Stress Utilisation: %g%%', i, n_1, n_2, delta_len, force, stress, strain, stress_utilisation);
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
scatter(nodes(:, 1), nodes(:, 2), 120, [0.8500 0.3250 0.0980], 'filled');
hold on;
for i = 1:elementCount
    [node1, node2] = getNodesByElement(i, data);    
    [x1, y1] = getNodePos(node1, false, data);
    [x2, y2] = getNodePos(node2, false, data);
    [new_x1, new_y1] = getNodePos(node1, true, data);
    [new_x2, new_y2] = getNodePos(node2, true, data);
     
    broken = false;
    if element_res{i, 1}.stress_utilisation > 100; color = [0.6350 0.0780 0.1840]; broken = true; else; color = [0 0.4470 0.7410]; end;
    
    line([new_x1, new_x2], [new_y1, new_y2], 'Color', "magenta", 'LineWidth', 2);
    line([x1, x2], [y1, y2], 'Color', color, 'LineWidth', 3);

    
    text(x2, y2, sprintf('%d', node2), 'FontSize', 15, 'Color', "#0BEED9");
    text(x1, y1, sprintf('%d', node1), 'FontSize', 15, 'Color', "#0BEED9");

    mid_x = (x1 + x2) / 2;
    mid_y = (y1 + y2) / 2;

    text(mid_x, mid_y, sprintf('%d', i), 'FontSize', 11, 'Color', "#000000");
end

for i = 1:supportCount
    [node, support_type, support_orien] = getSupportInfo(i, data);
    [x, y] = getNodePos(node, false, data);
    
    if support_type == 1
        color = 'b';
        h = 0.2;
        w = h * sqrt(2);
        x_points = [x-w/2, x, x+w/2];
        y_points = [y-h/2, y+h/2, y-h/2];
    elseif support_type == 2
        if support_orien == 1
            color = 'y';
        else 
            color = 'k';
        end
        h = 0.2;
        w = h * sqrt(2);
        x_points = [x-w/2, x, x+w/2];
        y_points = [y-h/2, y+h/2, y-h/2];
    end

    patch(x_points, y_points, color, 'EdgeColor', 'none');
end

title('STRUCTURE: ');
grid on;
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

function [support_node, support_type, support_orien] = getSupportInfo(i, data)
    support_node = data{1, 1}.supports(i, 1);
    support_type = data{1, 1}.supports(i, 2);
    support_orien = data{1, 1}.supports(i, 3);
end

function [force_node, fx, fy] = getForceInfo(i, data)
    force_node = data{1, 1}.forces(i, 1);
    force_len = data{1, 1}.forces(i, 2);
    force_angle = data{1, 1}.forces(i, 3) / 180 * pi;

    fx = force_len * cos(force_angle);
    fy = force_len * sin(force_angle);
end

function [element_lengths, element_angles, area, young_modulus, delta_f_mat] = getElementInformation(i, n_1, n_2, data)
    [x_1, y_1] = getNodePos(n_1, false, data);
    [x_2, y_2] = getNodePos(n_2, false, data);
    delta_x = x_2 - x_1;
    delta_y = y_2 - y_1;
    
    element_lengths = sqrt((delta_x^2) +(delta_y^2));
    element_angles = atan2(delta_y, delta_x);
    
    area = data{1, 1}.elements(i, 3);
    young_modulus = data{1, 1}.elements(i, 4);
    
    s_i_n = sin(element_angles);
    c_o_s = cos(element_angles);

    delta_f_mat = [c_o_s^2 s_i_n*c_o_s -c_o_s^2 -s_i_n*c_o_s;s_i_n*c_o_s s_i_n^2 -s_i_n*c_o_s -s_i_n^2;-c_o_s^2 -s_i_n*c_o_s c_o_s^2 s_i_n*c_o_s;-s_i_n*c_o_s -s_i_n^2 s_i_n*c_o_s s_i_n^2];
end

function [x_node, y_node] = getNodePos(node, newPos, data)
    if newPos == true
        new_nodes = data{1, 1}.newNodesLoc;
        x_node = new_nodes(node, 1);
        y_node = new_nodes(node, 2);
        return
    end
    nodes = data{1, 1}.nodes;
    x_node = nodes(node, 1);
    y_node = nodes(node, 2);
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