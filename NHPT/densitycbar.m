clear my_cbar

% my_cbar = ones(255, 3);
% 
% my_cbar(1:85, 3) = (84:-1:0)/84;
% my_cbar(86:255, 3) = 0;
% my_cbar(86:170, 2) = (84:-1:0)/84;
% my_cbar(171:255, 2) = 0;
% my_cbar(171:255, 1) = (84:-1:0)/84;
% colormap(my_cbar);
% or just invert cubehelix...

% or just invert cubehelix...

my_cbar = cubehelix([], [2.75,-0.1,2.5,1]);
my_cbar = my_cbar(end:-1:1,:);


colormap(my_cbar);