% This script implements the CFS algorithm to solve a problem specified in
% problem.m
addpath('lib')

%% Initialize the problem
Multi_car_problem

%% The Iteration
MAX_ITER = 100;
%figure; %hold on
tic
mini_distance = 5;
A = zeros(nstep + 4*(nstep-1), dim * nstep);
b = zeros(nstep + 4*(nstep-1), 1);
fig_1 = figure()
for k = 1:MAX_ITER
    % The constraints
    for i = 1 : nstep
        x_ref_1 = refpath(dim * (i - 1) + 1 : dim * i);
        if i < nstep
            x_ref_2 = refpath(dim * (i) + 1 : dim * (i + 1));
        end
        x_ref = get_ref(x_ref_1(1:2), x_ref_1(3:4), mini_distance);
        
        A_step = 2 * x_ref' * I_2;
        A(-4 + 5*i, dim * (i - 1) + 1 : dim * i) = -A_step;
        b(-4 + 5*i) = -mini_distance^2 - 1/2 * A_step * x_ref;
        if i < nstep
            [k1,b1] = get_line(x_ref_1(1), x_ref_1(2), x_ref_2(1), x_ref_2(2));
            [k2,b2] = get_line(x_ref_1(3), x_ref_1(4), x_ref_2(3), x_ref_2(4));
            
            flag = 0; %If flag = 1, then the reference is in a cross.
            
            %Set constraint for car 2 based on reference path of car 1;
            %If two points of car 2 are all in one side, do nothing
            if k1 * x_ref_1(3:4) + b1 > 0 && k1 * x_ref_2(3:4) + b1 > 0
                ;
            elseif k1 * x_ref_1(3:4) + b1 < 0 && k1 * x_ref_2(3:4) + b1 < 0
                ;

            %If the two cars' trajectory cross, let car 2 passes while car
            %1 remain in the side where the first point of car 1 is in. 
            elseif k1 * x_ref_1(3:4) + b1 > 0 && k1 * x_ref_2(3:4) + b1 < 0 && (k2 * x_ref_1(1:2) + b2) * (k2 * x_ref_2(1:2) + b2) < 0
                flag = 1;
                A(-2 + 5*i, dim * (i - 1) + dim/2 + 1 : dim * i) = -k1;  %<0
                A(0 + 5*i, dim * i + dim/2 + 1 : dim * (i + 1)) = k1;  %>0
                b(-2 + 5*i) = b1;
                b(0 + 5*i) = -b1;
                
                if k2 * x_ref_1(1:2) + b2 > 0
                    A(-3 + 5*i, dim * (i - 1) + 1 : dim * (i - 1) + dim/2) = -k2;
                    A(-1 + 5*i, dim * i + 1 : dim * i + dim/2) = -k2;  %>0 
                    b(-3 + 5*i) = b2;
                    b(-1 + 5*i) = b2;
                elseif k2 * x_ref_1(1:2) + b2 < 0
                    A(-3 + 5*i, dim * (i - 1) + 1 : dim * (i - 1) + dim/2) = k2;
                    A(-1 + 5*i, dim * i + 1 : dim * i + dim/2) = k2;  %>0 
                    b(-3 + 5*i) = -b2;
                    b(-1 + 5*i) = -b2;
                end
                 
            elseif k1 * x_ref_1(3:4) + b1 < 0 && k1 * x_ref_2(3:4) + b1 > 0 && (k2 * x_ref_1(1:2) + b2) * (k2 * x_ref_2(1:2) + b2) < 0
                flag = 1;
                A(-2 + 5*i, dim * (i - 1) + dim/2 + 1 : dim * i) = k1;  %<0
                A(0 + 5*i, dim * i + dim/2 + 1 : dim * (i + 1)) = -k1;  %>0
                b(-2 + 5*i) = -b1;
                b(0 + 5*i) = b1;
                
                if k2 * x_ref_1(1:2) + b2 > 0
                    A(-3 + 5*i, dim * (i - 1) + 1 : dim * (i - 1) + dim/2) = -k2;
                    A(-1 + 5*i, dim * i + 1 : dim * i + dim/2) = -k2;  %>0 
                    b(-3 + 5*i) = b2;
                    b(-1 + 5*i) = b2;
                elseif k2 * x_ref_1(1:2) + b2 < 0
                    A(-3 + 5*i, dim * (i - 1) + 1 : dim * (i - 1) + dim/2) = k2;
                    A(-1 + 5*i, dim * i + 1 : dim * i + dim/2) = k2;  %>0 
                    b(-3 + 5*i) = -b2;
                    b(-1 + 5*i) = -b2;
                end
            %If car 2 cross car 1's two points line, but car 1 doesn't cross car 2's
            %two points line. No cross constraints are need. Just set b as
            %1 to satisfy 0 < 1 as the row of A is zero.
            else
                ;
            end
        
        end
    end
    
    
    
%     Lstack = []; Sstack = [];
%     for i=1:nstep
%         for j=1:nobj
%             poly = obs{j}.poly+obs{j}.v*ones(1,4)*dt*i;
%             [L,S,d] = d2poly(refpath((i-1)*dim+1:i*dim)',poly');
%             Lstack = [Lstack;zeros(1,(i-1)*dim) L zeros(1,(nstep-i)*dim) zeros(1,nstep)];
%             Sstack = [Sstack;S-margin];
%         end
%     end
    

    
    %A = refpath * I_2
    % QP
    Qe = Qref+Qabs; %Shouldn't the x for Qref and Qabs be different?
    %soln = quadprog(Qe,-Qref*oripath,Lstack(:,1:size(Qe,1)),Sstack,Aeq(:,1:size(Qe,1)),beq);
    
    %tic
    soln = quadprog(Qe,-Qref*oripath,A,b,Aeq(:,1:size(Qe,1)),beq); %zeros(length(refpath), 1)
    %time_ite = toc
    
    pathnew = soln(1:dim*nstep);
    diff = norm(refpath-pathnew);
    
    %This block checks the distance constraint:
    Ax = 1/2 * A * (pathnew) + A * (pathnew - refpath)
    Ax_1 = 1/2 * A * (pathnew)
    DD = 1/2 * A * (pathnew)
    distance = zeros(nstep, 1);
    Dia = []
    for i = 1:nstep
    %     pathnew(4 * (i - 1) + 1 : 2 * i)
    %     pathnew(4 * i + 1 : 2 * (i + 1))
        distance(i) = norm(pathnew(4 * (i - 1) + 1 : 4 * (i - 1) + 2) - pathnew(4 * (i - 1) + 3 : 4 * i));
        Dia = [Dia; pathnew(4 * (i - 1) + 1 : 4 * i)' * I_2 * pathnew(4 * (i - 1) + 1 : 4 * i)];
    end
    
    
    if diff < 0.001*nstep*(dim)
        disp(['converged at step ',num2str(k)]);
        break
    end
    refpath = pathnew;
    %Draw the figures of results after 1 iteration.
    for i = 1:dim/2
        plot(pathnew(2 * (i - 1) + 1 : dim : end),pathnew(2 * i : dim : end),'LineWidth',2)
        hold on
    end
    legend('Car 1', 'Car 2')
    clf(fig_1)
    
end
%% Statistics
time = toc
disp('final cost: ');
%cost_quadratic(pathnew,oripath,Qref,Qabs)

%% Plot final solution
figure(1)
for i = 1:dim/2
    plot(pathnew(2 * (i - 1) + 1 : dim : end),pathnew(2 * i : dim : end),'LineWidth',2)
    hold on
end
title('Trajectory of Two Cars')


distance = zeros(nstep, 1);
for i = 1:nstep
%     pathnew(4 * (i - 1) + 1 : 2 * i)
%     pathnew(4 * i + 1 : 2 * (i + 1))
    distance(i) = norm(pathnew(4 * (i - 1) + 1 : 4 * (i - 1) + 2) - pathnew(4 * (i - 1) + 3 : 4 * i));
end
D_display = sprintf('The minimal distance between two cars is: %.2f.',min(distance));
disp(D_display)

figure(2)
plot(distance)
title('Distance Between Two Cars')




%% Plot obstacles
% for i=1:nobj
%     ob = Polyhedron('V',obs{i}.poly');
%     h = ob.plot('color', [0.8 0.8 0.8], 'EdgeColor','k');
%     alpha(h,0.5);
% end
% axis equal
% grid off
% box on
% legend('Iter1','Iter2','Iter3','Iter4','Iter5')
% axis([-0.5 9.5 -3 3])
%% Draw the gif
gif_flag = draw_gif(pathnew, dim, nstep)
function gif_out = draw_gif(pathnew, dim, nstep)
h = figure;
axis tight manual % this ensures that getframe() returns a consistent size
filename = 'testAnimated.gif';
for n = 1:1:nstep
    % Draw plot for two cars.
    x1 = pathnew(2 * (1 - 1) + 1 + (n - 1) * dim);
    y1 = pathnew(2 * 1 + (n - 1) * dim);
    x2 = pathnew(2 * (2 - 1) + 1 + (n - 1) * dim);
    y2 = pathnew(2 * 2 + (n - 1) * dim);
    plot(x1,y1,'bo-','MarkerSize',5) 
    hold on
    plot(x2,y2,'ro-','MarkerSize',5)
    axis([0 22 0 22])
    legend('Car 1', 'Car 2')
    drawnow 
      % Capture the plot as an image 
      frame = getframe(h); 
      im = frame2im(frame); 
      [imind,cm] = rgb2ind(im,256); 
      % Write to the GIF File 
      if n == 1 
          imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
      else 
          imwrite(imind,cm,filename,'gif','WriteMode','append'); 
      end 
end
gif_out = 1;
end
