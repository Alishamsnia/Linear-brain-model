%% Input

clc;
close all;
clear all;

% first we consider a plane for inputs ax+by+dz+e=0
a = 2;
b = -4;
d = 8;
e = -0;

% we define 100 random points on this plane
rand_points = [100,1];
x = 10 * rand(rand_points);
y = 10 * rand(rand_points);
% we define z according to x & y
z = (1/d) * (-a*x - b*y - e);

% we show scatter of this points 
figure("Name",'input')
input = [x y z];
scatter3(x,y,z,"MarkerFaceColor",'b');
title('input for 100 points on 2x+3y-4z-1=0 plane')
xlabel('x');
ylabel('y');
zlabel('z');
%% 
%% Trainig data with LMS algorythm
w_1 = rand(3,3);
iteration = 20:20:100; %matrix of number of iteration
step = 0.0001; % learning step
for i = 1:length(iteration)
    %weight matrix
    w = zeros(3,3,iteration(i));
    w( : , : , 1 ) = step * w_1; % initial weight

    %preallocating output and error
    output = zeros(100 , 3 , iteration(i));
    error = zeros(100 , 3 , iteration(i));
    for j = 1 : iteration(i)
        output( : , : , j ) = input( : , :) * w ( : , : , j );
        error( : , : , j ) = output( : , : , j ) - input( : , : );
        w( : , : , j+1 ) =  w( : , : , j ) + step *  input.' * error( : , : , j ); % Updating weight matrix
    end
    number_iteration = sprintf('%d',iteration(i));
    weight_matrix = w( : , : , iteration(i));
    MAError(i) = mae(error(:,:,iteration(i)));
    figure("Name",sprintf('iteration:%d',iteration(i)));
    scatter3(x,y,z,"MarkerFaceColor",'b');
    hold on
    scatter3(output( : , 1 , iteration(i)) , output( : , 2 , iteration(i)) , output( : , 3 , iteration(i)),'Marker','o',"MarkerEdgeColor",'r' )
    title(sprintf('number of Iterations = %d',iteration(i)))
    xlabel('x')
    ylabel('y')
    zlabel('z')
    legend('training inputs', 'training outputs')
    hold off
end


%%
% *Testing data*

%  Testing input
rand_points = [100,1];
x_test = 20 * rand(rand_points);
y_test = 20 * rand(rand_points);
% we define z according to x & y
z_test = (1/d) * (-a*x_test - b*y_test - e);
input_test = [x_test y_test z_test]

iteration = 20:20:100 % Matrix of numbers of iterations
MAError_t = zeros(1,length(iteration));

for i= 1:length(iteration)
    
    % Preallocating output and error
    output_test = zeros(100,3,iteration(i));
    error_test = zeros(100,3,iteration(i));
    
    output_test(:,:,iteration(i)) = input_test(:,:) * weight_matrix(:,:,iteration(i))  ;
    error_test(:,:,iteration(i)) = input_test(:,:) - output_test(:,:,iteration(i));
    
    figure(i+6)
    scatter3(x_test , y_test , z_test , 'marker', 'o' , "MarkerFaceColor", 'g', 'MarkerEdgeColor','g')
    hold on
    scatter3(output_test(:,1,iteration(i)),output_test(:,2,iteration(i)), output_test(:,3,iteration(i)), 'marker', 'o' ,'MarkerEdgeColor','c')
    title(sprintf('number of Iterations = %d',iteration(i)))
    xlabel('x')
    ylabel('y')
    zlabel('z')
    legend('test inputs', 'test outputs')
    hold off
    MAError_t(i) = mae(error_test(:,:,iteration(i)));
    
end

figure('name','validation')
plot(iteration,MAError)
hold on
plot(iteration,MAError_t)
grid on
title('Validation and Training Diagram')
xlabel('Number of Iterations')
ylabel('Mean Absolute Error')
legend('Training', 'Validation')







