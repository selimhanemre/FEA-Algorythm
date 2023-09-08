% AE 464 Project
% Selimhan Dagtas 2248052
% Mehmet Can Yalcin 2265080

clear;clc
coord = [0 0; 1 0; 2 0; 3 0;
         0 1; 1 1; 2 1; 3 1;
         0 2; 1 2; 2 2; 3 2];

x_coord = coord(:,1);
y_coord = coord(:,2);

input = menu("Choose mesh type","Triangular mesh","Quadrilateral mesh","Analytic Solution","Combined Solution");

if input == 1 || input == 4 % Triangular
    connectivity_triangular = [1 2 6;  6 5 1;
        2 3 7;  7 6 2;
        3 4 8;  8 7 3;
        5 6 10; 10 9 5;
        6 7 11; 11 10 6;
        7 8 12; 12 11 7;];

    Size_triangular = size(connectivity_triangular);
    numberOfElements_triangular = Size_triangular(1);
    xe = zeros(numberOfElements_triangular,3); ye = zeros(numberOfElements_triangular,3);
    part_triangular = cell(1,numberOfElements_triangular);

    for e = 1:numberOfElements_triangular


        xe(e,:) = [x_coord(connectivity_triangular(e,1)),x_coord(connectivity_triangular(e,2)),x_coord(connectivity_triangular(e,3))];
        ye(e,:) = [y_coord(connectivity_triangular(e,1)),y_coord(connectivity_triangular(e,2)),y_coord(connectivity_triangular(e,3))];

        part_triangular{e} = elem_stiffness_triangular(xe(e,:),ye(e,:));

    end

    K_global_triangular=zeros(12,12);
    for e = 1 : numberOfElements_triangular
        for i=1:length(connectivity_triangular(1,:))
            for j=1:length(connectivity_triangular(1,:))
                I_triangular = connectivity_triangular(e,i);J_triangular = connectivity_triangular(e,j);
                K_global_triangular(I_triangular,J_triangular) = K_global_triangular(I_triangular,J_triangular)+part_triangular{e}(i,j);
            end
        end

    end


    syms T1_triangular T2_triangular T3_triangular T5_triangular T6_triangular T7_triangular Q4_triangular ...
        Q8_triangular Q9_triangular Q10_triangular Q11_triangular Q12_triangular

    T_triangular = [T1_triangular T2_triangular T3_triangular 0 T5_triangular ...
        T6_triangular T7_triangular 0 100 100*cos(pi/6) 100*cos(2*pi/6) 0];

    Q_triangular = [0 0 0 Q4_triangular 0 0 0 Q8_triangular ...
        Q9_triangular Q10_triangular Q11_triangular Q12_triangular];

    BCS_triangular = [T_triangular;Q_triangular];

    F_init_triangular=cell(1,12);
    for i=1:12

        F_init_triangular{i}=sym(K_global_triangular(i,:).*T_triangular);

    end

    for i=1:12
        F_second_triangular(i,:)=[F_init_triangular{i}];
    end
    for i=1:12
        F_triangular(i)=sum(F_second_triangular(i,:));
    end

    F_triangular=F_triangular-Q_triangular;

    Answer_triangular = (solve(F_triangular(1,1),F_triangular(1,2),F_triangular(1,3),F_triangular(1,4),...
        F_triangular(1,5),F_triangular(1,6),F_triangular(1,7),F_triangular(1,8),F_triangular(1,9),...
        F_triangular(1,10),F_triangular(1,11),F_triangular(1,12)));

    % Edits to make it look good

    TabularAnswer_triangular=struct2table(Answer_triangular);

    Names_triangular = TabularAnswer_triangular.Properties.VariableNames;
    Answers_triangular=double(struct2array(Answer_triangular));


end

if input == 2 || input == 4 % Quadrilateral
    connectivity_quad = [1 2 6 5;
        2 3 7 6;
        3 4 8 7;
        5 6 10 9;
        6 7 11 10;
        7 8 12 11;];

    Size_quad = size(connectivity_quad);
    numberOfElements_quad = Size_quad(1);
    part_quad = cell(1,numberOfElements_quad);

    for i = 1:numberOfElements_quad
        part_quad{i} = round(double(subs((elem_stiffness_quad(x_coord(connectivity_quad(i,1)),x_coord(connectivity_quad(i,3)), ...
            y_coord(connectivity_quad(i,1)), y_coord(connectivity_quad(i,3)))))),3);
    end

    K_global_quad=zeros(12,12);
    for e = 1 : numberOfElements_quad
        for i=1:length(connectivity_quad(1,:))
            for j=1:length(connectivity_quad(1,:))
                I = connectivity_quad(e,i);J = connectivity_quad(e,j);
                K_global_quad(I,J) = K_global_quad(I,J)+part_quad{e}(i,j);
            end
        end

    end


    syms T1_quad T2_quad T3_quad T5_quad T6_quad T7_quad Q4_quad Q8_quad Q9_quad Q10_quad Q11_quad Q12_quad

    T_quad = [T1_quad T2_quad T3_quad 0 T5_quad T6_quad T7_quad 0 100 100*cos(pi/6) 100*cos(2*pi/6) 0];

    Q_quad = [0 0 0 Q4_quad 0 0 0 Q8_quad Q9_quad Q10_quad Q11_quad Q12_quad];

    BCS_quad = [T_quad;Q_quad];

    F_init_quad=cell(1,12);
    for i=1:12

        F_init_quad{i}=sym(K_global_quad(i,:).*T_quad);

    end

    for i=1:12
        F_second_quad(i,:)=[F_init_quad{i}];
    end
    for i=1:12
        F_quad(i)=sum(F_second_quad(i,:));
    end

    F_quad=F_quad-Q_quad;

    Answer_quad = (solve(F_quad(1,1),F_quad(1,2),F_quad(1,3),F_quad(1,4),F_quad(1,5),...
        F_quad(1,6),F_quad(1,7),F_quad(1,8),F_quad(1,9),F_quad(1,10),F_quad(1,11),F_quad(1,12)));

    % Edits to make it look good

    TabularAnswer_quad=struct2table(Answer_quad);

    Names_quad = TabularAnswer_quad.Properties.VariableNames;
    Answers_quad=double(struct2array(Answer_quad));

end

if input == 3 || input == 4
    Temp_analytic=zeros;
    for i=1:12
        Temp_analytic(i) = round(100*(cosh(pi/6*y_coord(i)) *cos(pi/6*x_coord(i))) / (cosh(pi/3)) ,3);
    end

end



% Graphs & Tables
NodeId=1:12;

if input == 1 || input == 4

    Temp_triangular = [Answers_triangular(7) Answers_triangular(8) Answers_triangular(9) 0;
        Answers_triangular(10) Answers_triangular(11) Answers_triangular(12) 0;
        100 100*cos(pi/6) 100*cos(2*pi/6) 0];

    figure
    contourf(Temp_triangular)
    colorbar;
    grid on;
    xlabel('x-axis');
    ylabel('y-axis');
    title('Triangular Mesh Temperature Distribution')
    Temp_triangular2 = [Answers_triangular(7) Answers_triangular(8) Answers_triangular(9) 0 ...
        Answers_triangular(10) Answers_triangular(11) Answers_triangular(12) 0 ...
        100 100*cos(pi/6) 100*cos(2*pi/6) 0];
    Q_triangular2 = [0 0 0 Answers_triangular(4) 0 0 0 Answers_triangular(5) ...
        Answers_triangular(6) Answers_triangular(1) Answers_triangular(2) Answers_triangular(3)];
    figure
    plot(y_coord,Q_triangular2)
    grid on
    xlabel('y coordinates')
    ylabel('Heat flux')
    title('Triangular Mesh Heat Flux Variation')
end

if input == 2 || input == 4

    Temp_quad= [Answers_quad(7) Answers_quad(8) Answers_quad(9) 0;
        Answers_quad(10) Answers_quad(11) Answers_quad(12) 0;
        100 100*cos(pi/6) 100*cos(2*pi/6) 0];

    figure
    contourf(Temp_quad)
    colorbar;
    grid on;
    xlabel('x-axis');
    ylabel('y-axis');
    title('Quadrilateral Mesh Temperature Distribution')

    Temp_quad2= [Answers_quad(7) Answers_quad(8) Answers_quad(9) 0 ...
        Answers_quad(10) Answers_quad(11) Answers_quad(12) 0 ...
        100 100*cos(pi/6) 100*cos(2*pi/6) 0];
    Q_quad2 = [0 0 0 Answers_quad(4) 0 0 0 Answers_quad(5) Answers_quad(6)...
        Answers_quad(1) Answers_quad(2) Answers_quad(3)];

    figure
    plot(y_coord,Q_quad2)
    grid on
    xlabel('y coordinates')
    ylabel('Heat flux')
    title('Quadrilateral Mesh Heat Flux Variation')

end


if input == 3 || input == 4

    Temp_analytic1= [Temp_analytic(1) Temp_analytic(2) Temp_analytic(3) Temp_analytic(4);
        Temp_analytic(5) Temp_analytic(6) Temp_analytic(7) Temp_analytic(8);
        Temp_analytic(9) Temp_analytic(10) Temp_analytic(11) Temp_analytic(12)];

    figure
    contourf(Temp_analytic1)
    colorbar;
    grid on;
    xlabel('x-axis');
    ylabel('y-axis');
    title('Analytic solution Temperature Distribution')

    Temp_analytic2= [Temp_analytic(1) Temp_analytic(2) Temp_analytic(3) Temp_analytic(4)...
        Temp_analytic(5) Temp_analytic(6) Temp_analytic(7) Temp_analytic(8)...
        Temp_analytic(9) Temp_analytic(10) Temp_analytic(11) Temp_analytic(12)];
end
if input == 4
    varname = {'Node_ID';'Triangular_Mesh';'Quadrilateral_Mesh';'Analytic_Solution'};
    TabularForm= table(NodeId',Temp_triangular2',Temp_quad2',Temp_analytic2','VariableNames',varname);

    disp(TabularForm)
end



function [K] = elem_stiffness_triangular(x,y)
A_e = abs((x(1)*(y(2)-y(3)))+(x(2)*(y(3)-y(1)))+(x(3)*(y(1)-y(2))))/2;
B_e = (1/(A_e*2))*[(y(2)-y(3)) (y(3)-y(1)) (y(1)-y(2)); (x(3)-x(2)) (x(1)-x(3)) (x(2)-x(1))];
K = B_e.'*B_e*A_e;
end

function [K] = elem_stiffness_quad(x1,x2,y1,y2)
A_e = (x2-x1)*(y2-y1);
syms x y
N_e =  (1/A_e)*[(x-x2)*(y-y2) -(x-x1)*(y-y2) (x-x1)*(y-y1) -(x-x2)*(y-y1)];
B_e = [diff(N_e,x) ; diff(N_e,y)];
I = B_e.'*B_e;
K = I;

for i=1:length(K(1,:))
    for j=1:length(K(:,1))
        K(i,j) = integral2(matlabFunction(K(i,j)),x1,x2,y1,y2);
    end
end


end

