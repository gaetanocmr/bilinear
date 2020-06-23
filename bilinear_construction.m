function [A,B,C,K_star] = bilinear_construction(displacement, force, tolerance, animation)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% This function will construct a bilinear curve from a pushover capacity  %
% curve. This uses the "area under the curve method" proposed by NTC2018  %
% Input:                                                                  %
% force = vector containing incremental force steps                       %
% displacement = vector containing incremental displacement steps         %
% tolerance = analysis tolerance expressed in percentage (for example 1)  %
% animation = 1 (print analysis iteration) / 2 (show the final results)   %
% NOTE:                                                                   %
% animation = 1 will save an iteration video named "bilinear_iter.mp4"    %
% animation = 3 NO PLOT													  %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% OUTPUT: A,B,C,K_star three points that defines the bilinear curve,      %
%         and bilinear elastic stiffness                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                                               Gaetano Camarda           %
%                                               V_1.5_BETA                %
%       gaetano.camarda@outlook.com             23/06/2020                %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%$
% Copyright 2020 Gaetano Camarda 
% For licence refer to: https://github.com/gaetanocmr/bilinear/blob/master/LICENSE

FS = [abs(force) abs(displacement)];
tol = tolerance;    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Si ANI
if animation == 1
    if size(FS,1) < 100
        for iter = 1 : 6
            t1 = iter;
            for i = 1 : size(FS,1) - 1
                t = i;
                x = (FS(i+1,1) + FS(i,1))/2;
                y = (FS(i+1,2) + FS(i,2))/2;
                FS_temp(1+(i-1)*2,1) = FS(i,1);
                FS_temp(1+(i-1)*2,2) = FS(i,2);
                FS_temp(2+(i-1)*2,1) = x;
                FS_temp(2+(i-1)*2,2) = y;
            end
        FS = FS_temp;
        clear FS_temp
        end
    end
    Fbu_star =  0.6 * max(FS(:,2));
    [minDistance, d_star_pos] = min(abs(FS(:,2)-ones(size(FS,1),1)*Fbu_star));
    d_star = FS(d_star_pos,1);
    A = [d_star Fbu_star];
    K_star = A(1,2) / A(1,1);
    figure
    plot(FS(:,1),FS(:,2),'b');
    hold on;
    grid on;
    plot(d_star,Fbu_star,'-o','MarkerEdgeColor','r');
    txt = ['   A = (',num2str(d_star),' ; ',num2str(Fbu_star),')'];
    text(d_star,Fbu_star,txt,'FontSize',8);
    curve = animatedline('Color','r','LineWidth',2);
    curve1 = animatedline('Color','r','LineWidth',2);
    title('BiLinear construction')
    xlabel('[m]');
    ylabel('[kN]');
    set(gca,'XLim',[0 max(FS(:,1))*1.2],'YLim',[0 max(FS(:,2))*1.1]);
    frame(1) = getframe(gcf);
    for d = 2 : size(FS,1);

        t = d;
        r_1 = @(x) K_star * x;
        A_r_1 = trapz(FS(1:d,1),K_star * (FS(1:d,1)));
        A_C1 = trapz(FS(1:d,1),FS(1:d,2));
        A1 = A_C1 - A_r_1;
        
        A_C2 = trapz(FS(d-1:end,1),FS(d-1:end,2));
        [minDistance, a_pos] = min(abs(FS(:,2)-ones(size(FS,1),1)*r_1(FS(d,1))));
        A_r_2 = trapz(FS(d-1:end,1),ones(size(FS(d-1:end,1),1),1)*FS(a_pos,2));
          
        A2 = A_C2 - A_r_2;
        save_A(d,:) = [A1 A2];

         perc_A1 = (abs(A1) / (abs(A1)+abs(A2))) * 100;
         perc_A2 = (abs(A2) / (abs(A1)+abs(A2))) * 100;
        if A2 < A1
            else if perc_A2 - perc_A1 < tol
                break
            end
        end
        if A1 < A2
            else if perc_A1 - perc_A2 < tol
                break
            end
        end

    addpoints(curve,FS(d,1),r_1(FS(d,1)));
    d_y = FS(d,1);
    pb = plot(d_y,FS(a_pos,2),'-o','MarkerEdgeColor','r');
    txt = ['   B = (',num2str(d_y),' ; ',num2str(FS(a_pos,2)),')'];
    pb_txt = text(d_y,FS(a_pos,2)+10,txt,'FontSize',8);
    pc = plot(FS(end,1),FS(a_pos,2),'-o','MarkerEdgeColor','r');
    txt = ['   C = (',num2str(FS(end,1)),' ; ',num2str(FS(a_pos,2)),')'];
    pc_txt = text(FS(end,1),FS(a_pos,2),txt,'FontSize',8);
    p = plot([FS(d,1),FS(end,1)],[FS(a_pos,2),FS(a_pos,2)],'r','LineWidth',2);
    drawnow;
    pause(0.02);
    frame(d) = getframe(gcf);
    delete(p);
    delete(pb);
    delete(pb_txt);
    delete(pc);
    delete(pc_txt);
    end
    d_y = FS(d,1);
    for i = 1 : 6
        plot(d_y,FS(a_pos,2),'-o','MarkerEdgeColor','r');
        txt = ['   B = (',num2str(d_y),' ; ',num2str(FS(a_pos,2)),')'];
        text(d_y,FS(a_pos,2)+6,txt,'FontSize',8);
        plot(FS(end,1),FS(a_pos,2),'-o','MarkerEdgeColor','r');
        txt = ['   C = (',num2str(FS(end,1)),' ; ',num2str(FS(a_pos,2)),')'];
        text(FS(end,1),FS(a_pos,2),txt,'FontSize',8);
        plot([FS(d,1),FS(end,1)],[FS(a_pos,2),FS(a_pos,2)],'r','LineWidth',2);
        frame(t) = getframe(gcf);
    end
    hold off
    videorinf = VideoWriter('bilinear_iter.mp4','MPEG-4');
    open(videorinf);
    writeVideo(videorinf,frame);
    close(videorinf);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NO ANI
if animation == 2
    if size(FS,1) < 100
        for iter = 1 : 6
        t1 = iter;
            for i = 1 : size(FS,1) - 1
                t = i;
                x = (FS(i+1,1) + FS(i,1))/2;
                y = (FS(i+1,2) + FS(i,2))/2;
                FS_temp(1+(i-1)*2,1) = FS(i,1);
                FS_temp(1+(i-1)*2,2) = FS(i,2);
                FS_temp(2+(i-1)*2,1) = x;
                FS_temp(2+(i-1)*2,2) = y;
            end
            FS = FS_temp;
            clear FS_temp
        end
    end
Fbu_star =  0.6 * max(FS(:,2));
[minDistance, d_star_pos] = min(abs(FS(:,2)-ones(size(FS,1),1)*Fbu_star));
d_star = FS(d_star_pos,1);
A = [d_star Fbu_star];
K_star = A(1,2) / A(1,1);
        
     for d = 2 : size(FS,1);

        t = d;
        r_1 = @(x) K_star * x;
        A_r_1 = trapz(FS(1:d,1),K_star * (FS(1:d,1)));
        A_C1 = trapz(FS(1:d,1),FS(1:d,2));
        A1 = A_C1 - A_r_1;
        
        A_C2 = trapz(FS(d-1:end,1),FS(d-1:end,2));
        [minDistance, a_pos] = min(abs(FS(:,2)-ones(size(FS,1),1)*r_1(FS(d,1))));
        A_r_2 = trapz(FS(d-1:end,1),ones(size(FS(d-1:end,1),1),1)*FS(a_pos,2));
          
        A2 = A_C2 - A_r_2;
        save_A(d,:) = [A1 A2];

         perc_A1 = (abs(A1) / (abs(A1)+abs(A2))) * 100;
         perc_A2 = (abs(A2) / (abs(A1)+abs(A2))) * 100;
        if A2 < A1
            else if perc_A2 - perc_A1 < tol
                break
            end
        end
        if A1 < A2
            else if perc_A1 - perc_A2 < tol
                break
            end
        end
end
    figure
    plot([0,FS(d,1)],[0,r_1(FS(d,1))],'r','LineWidth',2);
    hold on
    plot([FS(d,1),FS(end,1)],[FS(a_pos,2),FS(a_pos,2)],'r','LineWidth',2)
    hold on
    plot(FS(:,1),FS(:,2),'b')
    hold on
    plot(d_star,Fbu_star,'-o','MarkerEdgeColor','r')
    txt = ['   A = (',num2str(d_star),' ; ',num2str(Fbu_star),')'];
    text(d_star,Fbu_star,txt,'FontSize',8)
    hold on
    d_y = FS(d,1);
    plot(d_y,FS(a_pos,2),'-o','MarkerEdgeColor','r')
    txt = ['   B = (',num2str(d_y),' ; ',num2str(FS(a_pos,2)),')'];
    text(d_y,FS(a_pos,2)+6,txt,'FontSize',8)
    hold on
    plot(FS(end,1),FS(a_pos,2),'-o','MarkerEdgeColor','r')  
    txt = ['   C = (',num2str(FS(end,1)),' ; ',num2str(FS(a_pos,2)),')'];
    text(FS(end,1),FS(a_pos,2),txt,'FontSize',8)
    title('BiLinear construction')
    xlabel('[m]');
    ylabel('[kN]');
    set(gca,'XLim',[0 max(FS(:,1))*1.2],'YLim',[0 max(FS(:,2))*1.1]);
    hold off
    grid on
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NO PLOT
if animation == 3
    if size(FS,1) < 100
        for iter = 1 : 6
        t1 = iter;
            for i = 1 : size(FS,1) - 1
                t = i;
                x = (FS(i+1,1) + FS(i,1))/2;
                y = (FS(i+1,2) + FS(i,2))/2;
                FS_temp(1+(i-1)*2,1) = FS(i,1);
                FS_temp(1+(i-1)*2,2) = FS(i,2);
                FS_temp(2+(i-1)*2,1) = x;
                FS_temp(2+(i-1)*2,2) = y;
            end
            FS = FS_temp;
            clear FS_temp
        end
    end
Fbu_star =  0.6 * max(FS(:,2));
[minDistance, d_star_pos] = min(abs(FS(:,2)-ones(size(FS,1),1)*Fbu_star));
d_star = FS(d_star_pos,1);
A = [d_star Fbu_star];
K_star = A(1,2) / A(1,1);
        
     for d = 2 : size(FS,1);

        t = d;
        r_1 = @(x) K_star * x;
        A_r_1 = trapz(FS(1:d,1),K_star * (FS(1:d,1)));
        A_C1 = trapz(FS(1:d,1),FS(1:d,2));
        A1 = A_C1 - A_r_1;
        
        A_C2 = trapz(FS(d-1:end,1),FS(d-1:end,2));
        [minDistance, a_pos] = min(abs(FS(:,2)-ones(size(FS,1),1)*r_1(FS(d,1))));
        A_r_2 = trapz(FS(d-1:end,1),ones(size(FS(d-1:end,1),1),1)*FS(a_pos,2));
          
        A2 = A_C2 - A_r_2;
        save_A(d,:) = [A1 A2];

         perc_A1 = (abs(A1) / (abs(A1)+abs(A2))) * 100;
         perc_A2 = (abs(A2) / (abs(A1)+abs(A2))) * 100;
        if A2 < A1
            else if perc_A2 - perc_A1 < tol
                break
            end
        end
        if A1 < A2
            else if perc_A1 - perc_A2 < tol
                break
            end
        end
     end
    d_y = FS(d,1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END

A = [d_star Fbu_star];
B = [d_y FS(a_pos,2)];
C = [FS(end,1) FS(a_pos,2)];
end

