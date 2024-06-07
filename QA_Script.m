% Created by Justin Cooke
% Purpose of this script is to conduct quadrant analysis on the data from
% the dune field calculations

% Can Reproduce Figs. 7, 8, 10, 11, 12, 14, and 15 with this script.

%% Start

clc;
clear;
close all;

%% Load Data

myDir = dir('./Data/Sept13/x*');

Ndir = length(myDir);

ixyz = 0;
iu = 0;
iw = 0;
iuw = 0;


for i = 1:Ndir
    

    myName = myDir(i).name;
    myFolder = myDir(i).folder;
    
    myPath = strcat(myFolder,'/',myName);
    
    data = load(myPath);
    
    if contains(myPath,'README')
        ixyz = ixyz + 1;
        x(ixyz) = round(data(1,2));
        z_local{ixyz} = data(:,end);
    elseif contains(myPath,'comp(u,0)')
        iu = iu + 1;
        u{iu} = data(:,4:end);
    elseif contains(myPath,'comp(u,2)')
        iw = iw + 1;
        w{iw} = data(:,4:end);
    elseif contains(myPath,'comp(u_rey,1)')
        iuw = iuw + 1;
        uw{iuw} = data(:,4:end);
    end
        
end

clear my* Ndir i* data

%% Create a standard global z

nu = 0.32e-4;
delta_abl = 300;
u_tau_avg = 0.14;
delta_nu = nu/u_tau_avg;

z_global = linspace(1.5,400,100);
zplus_global = z_global./delta_nu;
zdelta_global = z_global./delta_abl;


%% Create Mean and Fluctuating Quantities 

N_u = length(u);
N_w = length(w);
N_uw = length(uw);

for i = 1:N_u

    tempu = u{i};
    U{i} = mean(tempu);
    up{i} = tempu - U{i};
    urms{i} = rms(up{i});
    
    tempw = w{i};
    W{i} = mean(tempw);
    wp{i} = tempw - W{i};
    wrms{i} = rms(up{i});
    
    tempuw = uw{i};
    UW{i} = mean(tempuw);
    
end
    
clear temp*

%% Quadrants

% Q1 == u'>0 & w'>0 (outward motions)
% Q2 == u'<0 & w'>0 (ejections)
% Q3 == u'<0 & w'<0 (inward motions)
% Q4 == u'>0 & w'<0 (sweeps)

% For each x-location, will want to create a product of u' and w' in time
% and at each z

[Nt,Nz] = size(up{1});

for i = 1:N_u

    tempup = up{i};
    tempwp = wp{i};
    
    for j = 1:Nz
        
        tempupz = tempup(:,j);
        tempwpz = tempwp(:,j);
        tempprod{j} = tempupz.*tempwpz;
        
    end
        
    temp = cell2mat(tempprod);
    uw_prod{i} = temp;
    
end

clear temp*

%% Time Stuff

dt = 0.0065*50; %0.0065*100;
tend = dt*Nt;
t = linspace(0,tend,Nt);
T = t*6/delta_abl;


%% Divide the multiplied uw matrix into quadrants

for i = 1:N_uw
   
    tempu = up{i};
    tempw = wp{i};
    
    for j = 1:Nz
        
        tempuz = tempu(:,j);
        tempwz = tempw(:,j);
        
        for k = 1:Nt
            
            if tempuz(k) > 0 && tempwz(k) > 0
                q_store(k) = 1;
            elseif tempuz(k) < 0 && tempwz(k) > 0
                q_store(k) = 2;
            elseif tempuz(k) < 0 && tempwz(k) < 0
                q_store(k) = 3;
            elseif tempuz(k) > 0 && tempwz(k) < 0
                q_store(k) = 4;
            end
            
        end
        
        temp_quad{j} = q_store';
        
    end
    store_quad = cell2mat(temp_quad);
    quadrants{i} = store_quad;
end

clear store_quad temp* q_store 

%% Partition events and times 


for i = 1:N_uw
   
    x_quad = quadrants{i};
    uw_quad = uw_prod{i};
    
    for j = 1:Nz
        
        sum_eject = 0;
        sum_sweep = 0;
        sum_q1    = 0;
        sum_q3    = 0;
        count_eject = 0;
        count_sweep = 0;
        count_q1    = 0;
        count_q3    = 0;
       
        
        for k = 1:Nt
            
            if x_quad(k,j) == 2
                add_me = uw_quad(k,j);
                sum_eject = sum_eject + add_me;
                count_eject = count_eject + 1;
                vec_eject(count_eject) = add_me;
            elseif x_quad(k,j) == 4
                add_me = uw_quad(k,j);
                sum_sweep = sum_sweep + add_me;
                count_sweep = count_sweep + 1;
                vec_sweep(count_sweep) = add_me;
            elseif x_quad(k,j) == 1
                add_me = uw_quad(k,j);
                sum_q1 = sum_q1 + add_me;
                count_q1 = count_q1 + 1;
                vec_q1(count_q1) = add_me;
            elseif x_quad(k,j) == 3
                add_me = uw_quad(k,j);
                sum_q3 = sum_q3 + add_me;
                count_q3 = count_q3 + 1;
                vec_q3(count_q3) = add_me; 
            end
            
            
        end
        
        q1{j} = vec_q1;
        q2{j} = vec_eject;
        q3{j} = vec_q3;
        q4{j} = vec_sweep;
        
        eject(j) = (sum_eject/count_eject);
        sweep(j) = (sum_sweep/count_sweep);
        outward(j) = (sum_q1/count_q1);
        inward(j) = (sum_q3/count_q3);

        clear vec_eject vec_sweep
        
    end
    
    Q1{i} = outward;
    Q2{i} = eject;
    Q3{i} = inward;
    Q4{i} = sweep;
    
    q1_x{i} = q1;
    q2_x{i} = q2;
    q3_x{i} = q3;
    q4_x{i} = q4;
       
end

total_events = count_eject + count_sweep + count_q1 + count_q3;

clear sum_* count_* add_me uw_quad x_quad

%% Summing contributions to <u'w'>

for i = 1:N_u
   
    for j = 1:Nz
       
        temp_q1 = abs(mean(q1_x{i}{j}));
        temp_q2 = abs(mean(q2_x{i}{j}));
        temp_q3 = abs(mean(q3_x{i}{j}));
        temp_q4 = abs(mean(q4_x{i}{j}));
        
        sum_qs = temp_q1 + temp_q2 + temp_q3 + temp_q4;
        contr_q1{i}(j) = temp_q1 / sum_qs;
        contr_q2{i}(j) = temp_q2 / sum_qs;
        contr_q3{i}(j) = temp_q3 / sum_qs;
        contr_q4{i}(j) = temp_q4 / sum_qs;
    end
    
end

for i = 1:N_u
   
    contr_q1{i} = contr_q1{i}';
    contr_q2{i} = contr_q2{i}';
    contr_q3{i} = contr_q3{i}';
    contr_q4{i} = contr_q4{i}';
    
end

q1_c = cell2mat(contr_q1);
q2_c = cell2mat(contr_q2);
q3_c = cell2mat(contr_q3);
q4_c = cell2mat(contr_q4);


clear temp_q* sum_qs clear contr_q*

%% Plot contributions

qColors = ["#daf8e3" "#ffc100" "#00c2c7" "#bf0000"];
zdelta_colors = ['#0c9bba','#44d6a9','#fafa70'];



close all;

figure()
plot(zdelta_global,q1_c(:,1),'ksquare-','MarkerFaceColor',qColors(1),...
        'MarkerSize',8,'LineWidth',2); hold on
    plot(zdelta_global,q2_c(:,1),'k^--','MarkerFaceColor',qColors(2),...
        'MarkerSize',8,'LineWidth',2); 
    plot(zdelta_global,q3_c(:,1),'ko:','MarkerFaceColor',qColors(3),...
        'MarkerSize',8,'LineWidth',2);
    plot(zdelta_global,q4_c(:,1),'kv-.','MarkerFaceColor',qColors(4),...
        'MarkerSize',8,'LineWidth',2);
set(gca,'FontName','SansSerif','FontSize',18);
    title('Alkali Flat','Interpreter','Latex','FontSize',24);
    xlabel('$z/\delta$',...
        'Interpreter','Latex','FontSize',24,'FontName','SansSerif');
    ylabel('$S_i$',...
        'Interpreter','Latex','FontSize',24,'FontName','SansSerif');
    xlim([0 1]);
    ylim([0 0.75]);
    
    
    
figure()
for i = 2:N_u-1
   
    nexttile;
    plot(zdelta_global,q1_c(:,i),'ksquare-','MarkerFaceColor',qColors(1),...
        'MarkerSize',8,'LineWidth',2); hold on
    plot(zdelta_global,q2_c(:,i),'k^--','MarkerFaceColor',qColors(2),...
        'MarkerSize',8,'LineWidth',2); 
    plot(zdelta_global,q3_c(:,i),'ko:','MarkerFaceColor',qColors(3),...
        'MarkerSize',8,'LineWidth',2);
    plot(zdelta_global,q4_c(:,i),'kv-.','MarkerFaceColor',qColors(4),...
        'MarkerSize',8,'LineWidth',2);
    makeATitle = strcat('$\hat{x} = $',num2str(round(x(i)) - 1850),' m'); 
    set(gca,'FontName','SansSerif','FontSize',18);
    title(makeATitle,'Interpreter','Latex','FontSize',24);
    xlabel('$z/\delta$',...
        'Interpreter','Latex','FontSize',24,'FontName','SansSerif');
    ylabel('$S_i$',...
        'Interpreter','Latex','FontSize',24,'FontName','SansSerif');
%     if i == N_u
%         legend('Outward Motions','Ejections','Inward Motions','Sweeps',...
%             'Location','EastOutside');
%     end
    xlim([0 1]);
    ylim([0 0.75]);
end


%% 

for i = 1:N_uw
   
    temp_q2 = q2_x{i};
    temp_q4 = q4_x{i};
    
    for j = 1:Nz

        q2_mean(j) = mean(temp_q2{j});
        q4_mean(j) = mean(temp_q4{j});
    
               
    end
    
    Q2_mean{i} = q2_mean;
    Q4_mean{i} = q4_mean;
    
end



%% Scatter Plots


close all;

for i = 2:1:N_u-1
   figure(i)
   p1 = scatter(up{i}(:,5),wp{i}(:,5),'^','MarkerEdgeColor','#2436db'); hold on
   p2 = scatter(up{i}(:,8),wp{i}(:,8),'square','MarkerEdgeColor','#fa007f');
   p3 = scatter(up{i}(:,23),wp{i}(:,23),'o','MarkerEdgeColor','#ff9a00');
   l1 = xline(0);
   l2 = yline(0);
   set(gca,'FontName','SansSerif','FontSize',18);
   xlabel('$u^\prime [m/s]$','Interpreter','Latex','FontName','SansSerif',...
       'FontSize',24);
   ylabel('$w^\prime [m/s]$','Interpreter','Latex','FontName','SansSerif',...
       'FontSize',24);
   grid on;
%    legend([p1,p2,p3],'z/\delta = 0.06','z/\delta = 0.1','z/\delta = 0.3');
   myTitle = strcat('$\hat{x} = $',num2str(round(x(i))-1850),' m');
   title(myTitle,'Interpreter','Latex','FontName','SansSerif',...
       'FontSize',24);
   ylim([-1 1]);
   xlim([-1 1]);
end



%% Partitioning quadrant events over the time 

% 1 = zdelta 0.01
% 5 = zdelta 0.05
% 8 = zdelta 0.1
% 23 = zdelta 0.3
% 38 = zdelta 0.5

%%%------------------------PICK Z Location HERE----------------------%%%
zselection = 8;

clear temp;

for i = 1:N_u

    i1 = find(quadrants{i}(:,zselection) == 1);
    frac_q1 = length(i1)/Nt;

    i2 = find(quadrants{i}(:,zselection) == 2);
    frac_q2 = length(i2)/Nt;

    i3 = find(quadrants{i}(:,zselection) == 3);
    frac_q3 = length(i3)/Nt;

    i4 = find(quadrants{i}(:,zselection) == 4);
    frac_q4 = length(i4)/Nt;

    temp{i} = [frac_q1 frac_q2 frac_q3 frac_q4];

end

temp = temp';
frac = cell2mat(temp);

%% Doing the above for all Z locations now

clear temp;

for j = 1:Nz

    for i = 1:N_u

        i1 = find(quadrants{i}(:,j) == 1);
        frac_q1 = length(i1)/Nt;

        i2 = find(quadrants{i}(:,j) == 2);
        frac_q2 = length(i2)/Nt;

        i3 = find(quadrants{i}(:,j) == 3);
        frac_q3 = length(i3)/Nt;

        i4 = find(quadrants{i}(:,j) == 4);
        frac_q4 = length(i4)/Nt;

        temp{i} = [frac_q1 frac_q2 frac_q3 frac_q4];

    end


z_frac{j} = (temp);

end

%% Creating matrices within the cell zFrac

for i = 1:Nz
    zFrac{i} = cell2mat(z_frac{i}');
end


%% Bar Chart for z = 0.01, 0.06, 0.1, and 0.3

close all;

theseZs = [1,5,8,23];
XDunes = [1,2,3,4,5,6,7,8,9,10];
SDunes = {'Alkali Flat','$\hat{x}_1$','$\hat{x}_2$','$\hat{x}_3$','$\hat{x}_4$',...
    '$\hat{x}_5$','$\hat{x}_6$','$\hat{x}_7$','$\hat{x}_8$','$\hat{x}_9$'};

figure();
for i = 1:4
    subplot(4,1,i);
    b = bar(XDunes,zFrac{theseZs(i)}(1:end-1,:));
    set(gca,'xtick',XDunes,'XTickLabel',SDunes,'TickLabelInterpreter','latex');
    b(1).FaceColor = "#daf8e3";
    b(2).FaceColor = "#ffc100";
    b(3).FaceColor = "#00c2c7";
    b(4).FaceColor = "#bf0000";
%     legend('Q1','Q2','Q3','Q4');
    createATitle = strcat('$z/\delta$ = ',num2str(round(zdelta_global(theseZs(i)),2)));
    title(createATitle,'Interpreter','Latex');
    ylabel('$N_Q$/$N_tot$','Interpreter','Latex');
    ylim([0 1]);
    set(gca,'FontSize',16);
end






%% Appending zeros to the end of the quadrants for integration bounds purposes

pad = zeros(1,Nz);

for i = 1:N_u
   
    tmp = quadrants{i};
    tmp0 = [tmp; pad];
    bnd_qds{i} = tmp0;
    
end

clear tmp*

%% Time duration bounds for single Z-Location

% select z/delta location for analysis
%zselection = 23;

for i = 1:N_u

    data = bnd_qds{i}(:,zselection); 
    
    mystep = 0;
    q1i = 0;
    q2i = 0;
    q3i = 0;
    q4i = 0;

    clear q*Bounds

    % while loop will loop over the whole data set and find the start and 
    % end points of each event, then bin by type
    while mystep <= Nt + 1 

        mystep = mystep + 1;
        comp = data(mystep);

        if mystep == 1
            event = data(mystep);
            A = mystep;
        end

        if comp == event

        else
            B = mystep - 1;
            if event == 1
                q1i = q1i + 1;
                q1Bounds{q1i} = [A B];
            elseif event == 2
                q2i = q2i + 1;
                q2Bounds{q2i} = [A B];
            elseif event == 3
                q3i = q3i + 1;
                q3Bounds{q3i} = [A B];
            elseif event == 4
                q4i = q4i + 1;
                q4Bounds{q4i} = [A B];
            end
            A = mystep;
            event = data(A);
        end

        if event == 0
            break;
        end
    end

    
    tmpq1 = cell2mat(q1Bounds');
    tmpq2 = cell2mat(q2Bounds');
    tmpq3 = cell2mat(q3Bounds');
    tmpq4 = cell2mat(q4Bounds');
    xbounds{i} = {tmpq1 tmpq2 tmpq3 tmpq4};
    
end

clear tmpq*

%% Now I will find the avg time-duration of each quadrant

% xbounds is set up like so: there are 10 cells, each representative of the
% 10 x-locations at which we probe. Then within each of these cells are an
% additional 4 cells which contain N_AB x 2 matrices with the A and B
% values for each event. Cell 1 corresponds to quad 1, and so on and so 
% forth here

% First we will add a value to each B end point

for i = 1:N_u    
    for j = 1:4
        XB{1,i}{1,j} = [xbounds{1,i}{1,j}(:,1) xbounds{1,i}{1,j}(:,end)...
            + 1];
    end
end

% Now we will go through and take the difference between A and B and
% multiply by dt 

k_a = 3; % avg dune height

for i = 1:N_u
   
    for j = 1:4
       
        thisA = XB{1,i}{1,j}(:,1);
        thisB = XB{1,i}{1,j}(:,end);
        
        BmAdt = (thisB - thisA).*dt;
        BmAdtm = mean(BmAdt);
        
        if j == 1
            T_q1 = BmAdtm;
        elseif j == 2
            T_q2 = BmAdtm;
        elseif j == 3
            T_q3 = BmAdtm;
        else
            T_q4 = BmAdtm;
        end
        
    end
    
    % Non-dim avg time duration of each quadrant event
    T_Q{i} = [T_q1 T_q2 T_q3 T_q4]*(U{i}(zselection)/k_a);
    
    
end

T_Qmat = cell2mat(T_Q');

%% Now need some figures to share the results

close all;

X = [1,2,3,4,5,6,7,8,9,10,11];
S = {'$\hat{x}_0$','$\hat{x}_1$','$\hat{x}_2$','$\hat{x}_3$','$\hat{x}_4$',...
    '$\hat{x}_5$','$\hat{x}_6$','$\hat{x}_7$','$\hat{x}_8$','$\hat{x}_9$',...
    '$\hat{x}_{10}$'};

figure();
subplot(2,1,1);
b = bar(X,T_Qmat);
set(gca,'xtick',X,'XTickLabel',S,'TickLabelInterpreter','latex');
b(1).FaceColor = "#daf8e3";
b(2).FaceColor = "#ffc100";
b(3).FaceColor = "#00c2c7";
b(4).FaceColor = "#bf0000";
myTitle = strcat('$z/\delta$ = ',num2str(round(zdelta_global(zselection),2)));
title(myTitle,'Interpreter','latex');
legend('Q1','Q2','Q3','Q4','Location','NorthWest');
ylabel('Normalized Average Event Duration $T^*_Q$','Interpreter','Latex');
set(gca,'FontSize',16);

subplot(2,1,2);
b_ = bar(X,frac);
set(gca,'xtick',X,'XTickLabel',S,'TickLabelInterpreter','latex');
b_(1).FaceColor = "#daf8e3";
b_(2).FaceColor = "#ffc100";
b_(3).FaceColor = "#00c2c7";
b_(4).FaceColor = "#bf0000";
legend('Q1','Q2','Q3','Q4');
ylabel('$\%$ Time in Each Quadrant','Interpreter','Latex');
ylim([0 1]);
set(gca,'FontSize',16);

%% Want to determine 'Impulse' as in Bristow 2021

% Will need to do some type of integration between each event 


for i = 1:N_u 
   
    thisUWProd = abs(uw_prod{i}(:,zselection));
        
    for j = 1:4
       thisBounds = XB{i}{j};
       boundLength = length(thisBounds); %XB{i}{j});
       
       for k = 1:boundLength
            Aqk = XB{i}{j}(k,1);
            Bqk = XB{i}{j}(k,2);
            if Bqk == Nt+1
                Bqk = Nt;
                Aqk = Nt-1;
            end
            temp = thisUWProd(Aqk:Bqk);
%             if Aqk == 1
%                 Aqk = 0;
%             end
            intOver = linspace(dt*Aqk,dt*Bqk,length(temp));
            myInt(k) = trapz(intOver,temp);
       end
       AvgInt(j) = mean(myInt);
       clear myInt
    end
    Jq{i} = (U{i}(zselection)/(k_a*u_tau_avg*u_tau_avg))*AvgInt;
    clear AvgInt
end

clear thisUWProd thisBounds boundLength Aqk Bqk temp AvgInt myInt intOver

JQ = cell2mat(Jq');

%% Now Let's plot everything

close all;

figure();
subplot(3,1,1);
b_ = bar(X,frac);
set(gca,'xtick',X,'XTickLabel',S,'TickLabelInterpreter','latex');
b_(1).FaceColor = "#daf8e3";
b_(2).FaceColor = "#ffc100";
b_(3).FaceColor = "#00c2c7";
b_(4).FaceColor = "#bf0000";
legend('Q1','Q2','Q3','Q4','Location','EastOutside');
myTitle = strcat('$z/\delta$ = ',num2str(round(zdelta_global(zselection),2)));
title(myTitle,'Interpreter','latex');
ylabel('$N_Q$/$N_tot$','Interpreter','Latex');
ylim([0 0.5]);
set(gca,'FontSize',16);

subplot(3,1,2);
b = bar(X,JQ);
set(gca,'xtick',X,'XTickLabel',S,'TickLabelInterpreter','latex');
b(1).FaceColor = "#daf8e3";
b(2).FaceColor = "#ffc100";
b(3).FaceColor = "#00c2c7";
b(4).FaceColor = "#bf0000";
legend('Q1','Q2','Q3','Q4','Location','EastOutside');
ylim([0 15]);
ylabel('$J^*_Q$','Interpreter','Latex');
set(gca,'FontSize',16);



subplot(3,1,3);
b = bar(X,T_Qmat);
set(gca,'xtick',X,'XTickLabel',S,'TickLabelInterpreter','latex');
b(1).FaceColor = "#daf8e3";
b(2).FaceColor = "#ffc100";
b(3).FaceColor = "#00c2c7";
b(4).FaceColor = "#bf0000";
legend('Q1','Q2','Q3','Q4','Location','EastOutside');
ylim([0 12]);
ylabel('$T^*_Q$','Interpreter','Latex');
set(gca,'FontSize',16);


%% Finding T_Q and J_Q for all z at each of the 11 streamwise locations



for i = 1:N_u

    for j = 1:Nz
        data = bnd_qds{i}(:,j); 

        mystep = 0;
        q1i = 0;
        q2i = 0;
        q3i = 0;
        q4i = 0;

        clear q*Bounds

        % while loop will loop over the whole data set and find the start and 
        % end points of each event, then bin by type
        while mystep <= Nt + 1 

            mystep = mystep + 1;
            comp = data(mystep);

            if mystep == 1
                event = data(mystep);
                A = mystep;
            end

            if comp == event

            else
                B = mystep - 1;
                if event == 1
                    q1i = q1i + 1;
                    q1Bounds{q1i} = [A B];
                elseif event == 2
                    q2i = q2i + 1;
                    q2Bounds{q2i} = [A B];
                elseif event == 3
                    q3i = q3i + 1;
                    q3Bounds{q3i} = [A B];
                elseif event == 4
                    q4i = q4i + 1;
                    q4Bounds{q4i} = [A B];
                end
                A = mystep;
                event = data(A);
            end

            if event == 0
                break;
            end
        end


        tmpq1 = cell2mat(q1Bounds');
        tmpq2 = cell2mat(q2Bounds');
        tmpq3 = cell2mat(q3Bounds');
        tmpq4 = cell2mat(q4Bounds');
        Z_xbounds{j} = {tmpq1 tmpq2 tmpq3 tmpq4};
    end
    
    ZX_bounds{i} = Z_xbounds;
    clear Z_xbounds
end

clear tmpq*

%% Now I will find the avg time-duration of each quadrant

% ZXbounds is set up like so: there are 11 cells, each representative of the
% 11 x-locations at which we probe. Then within each of these cells are an
% additional 100 cells with 4 cells that contain N_AB x 2 matrices with the 
% A and B values for each event. Cell 1 corresponds to quad 1, and so on 
% and so forth here

% First we will add a value to each B end point

for i = 1:N_u
    for k = 1:Nz
        for j = 1:4
            ZXB{1,i}{1,k}{1,j} = ...
                [ZX_bounds{1,i}{1,k}{1,j}(:,1) ZX_bounds{1,i}{1,k}{1,j}(:,end) + 1];
        end
    end
end

% Now we will go through and take the difference between A and B and
% multiply by dt 

k_a = 3; % avg dune height

for i = 1:N_u
   
    for k = 1:Nz
        for j = 1:4

            thisA = ZXB{1,i}{1,k}{1,j}(:,1);
            thisB = ZXB{1,i}{1,k}{1,j}(:,end);

            BmAdt = (thisB - thisA).*dt;
            BmAdtm = mean(BmAdt);

            if j == 1
                T_q1 = BmAdtm;
            elseif j == 2
                T_q2 = BmAdtm;
            elseif j == 3
                T_q3 = BmAdtm;
            else
                T_q4 = BmAdtm;
            end

        end
    
    % Non-dim avg time duration of each quadrant event
    ZX_T_Q{i}{k} = [T_q1 T_q2 T_q3 T_q4]*(U{i}(k)/k_a);
    
    end
    
end

%% Want to determine 'Impulse' for all Z now

% Will need to do some type of integration between each event 

clear thisUWProd thisBounds boundLength Aqk Bqk temp intOver myInt AvgInt ...
    ZJq ZXJq

for i = 1:N_u 
   
    for kk = 1:Nz
        thisUWProd = abs(uw_prod{i}(:,kk));

        for j = 1:4
           thisBounds = ZXB{i}{kk}{j};
           boundLength = length(thisBounds);

           for k = 1:boundLength
                Aqk = ZXB{i}{kk}{j}(k,1);
                Bqk = ZXB{i}{kk}{j}(k,2);
                if Bqk == Nt+1
                    Bqk = Nt;
                    Aqk = Nt-1;
                end
                temp = thisUWProd(Aqk:Bqk);
    %             if Aqk == 1
    %                 Aqk = 0;
    %             end
                intOver = linspace(dt*Aqk,dt*Bqk,length(temp));
                myInt(k) = trapz(intOver,temp);
           end
           AvgInt(j) = mean(myInt);
           clear myInt
        end
        
        ZJq{kk} = (U{i}(kk)/(k_a*u_tau_avg*u_tau_avg))*AvgInt;
        clear AvgInt
    end
    
    ZXJq{i} = ZJq;
    
end

%% Clean up data



for i = 1:N_u
    testM = ZXJq{i}';
    testc2m = cell2mat(testM);
    
    ZXJq_mat{i} = testc2m;
    
    testN = ZX_T_Q{i}';
    testNc2m = cell2mat(testN);
    
    ZXTq_mat{i} = testNc2m;
end

% Now things are set as 11 cells of streamwise locations with 100x4
% matrices in each cell, giving the value of T and J for each quadrant at
% each wall-normal location

% Now I am only going to have Q2 and Q4

for i = 1:N_u
   
    ZXJq2{i} = ZXJq_mat{i}(:,2);
    ZXTq2{i} = ZXTq_mat{i}(:,2);
    
    ZXJq4{i} = ZXJq_mat{i}(:,4);
    ZXTq4{i} = ZXTq_mat{i}(:,4);
end

clear thisUWProd thisBounds boundLength Aqk Bqk temp AvgInt myInt intOver
clear myLegend* mystep myZ myZ2 i1 i2 i3 i4 frac_q* A temp* test* this*

%% Creating Dual Bar Charts of Tq and Jq at z = 0.01, 0.06, 0.1, and 0.3

close all;

count = 0;
secCount = 0;

for j = 1:4
   for i = 1:N_u
    TqTemp{i} = ZXTq_mat{i}(theseZs(j),:);
    JqTemp{i} = ZXJq_mat{i}(theseZs(j),:);
   end
   
   TqPlot{j} = cell2mat(TqTemp');
   JqPlot{j} = cell2mat(JqTemp');
   
end

figure();
for i = 1:8
    
    subplot(4,2,i);
    if rem(i,2) ~= 0
        count = count + 1;
        b = bar(XDunes,TqPlot{count}(1:end-1,:));
        set(gca,'xtick',XDunes,'XTickLabel',SDunes,'TickLabelInterpreter','latex');
        b(1).FaceColor = "#daf8e3";
        b(2).FaceColor = "#ffc100";
        b(3).FaceColor = "#00c2c7";
        b(4).FaceColor = "#bf0000";
%         legend('Q1','Q2','Q3','Q4','Location','EastOutside');
        createATitle = strcat('$z/\delta$ = ',num2str(round(zdelta_global(theseZs(count)),2)));
        title(createATitle,'Interpreter','Latex');
        ylim([0 12]);
        ylabel('$T^*_Q$','Interpreter','Latex');
        set(gca,'FontSize',16);
    else
        secCount = secCount + 1;
        b = bar(XDunes,JqPlot{secCount}(1:end-1,:));
        set(gca,'xtick',XDunes,'XTickLabel',SDunes,'TickLabelInterpreter','latex');
        b(1).FaceColor = "#daf8e3";
        b(2).FaceColor = "#ffc100";
        b(3).FaceColor = "#00c2c7";
        b(4).FaceColor = "#bf0000";
%         legend('Q1','Q2','Q3','Q4','Location','EastOutside');
        createATitle = strcat('$z/\delta$ = ',num2str(round(zdelta_global(theseZs(secCount)),2)));
        title(createATitle,'Interpreter','Latex');
        ylim([0 20]);
        ylabel('$J^*_Q$','Interpreter','Latex');
        set(gca,'FontSize',16);
    end
end


%% Just Impulse on a Single Plot 

close all;

delta_ibl = [0.1,0.06356,0.1665,0.1665,0.1373,0.1372,0.2968,0.4362,...
    0.4362,0.5289].*delta_abl;

iblColors = ["#000000", "#daf8e3", "#97ebdb", "#00c2c7", "#0086ad", "#005582", ...
    "#ffc100", "#ff9a00", "#ff7400", "#bf0000" ];

figure();
for i = 1:length(delta_ibl)
   
    if i == 1 || i == 2
        plot(z_global./30,ZXJq2{i},'-','Color',...
            iblColors(i),'LineWidth',2); hold on
        plot(z_global./30,ZXJq4{i},'-','Color',...
            iblColors(i),'LineWidth',2); 
    else
        plot(z_global./(delta_ibl(i)),ZXJq2{i},'-','Color',...
                iblColors(i),'LineWidth',2); hold on
        plot(z_global./(delta_ibl(i)),ZXJq4{i},'--','Color',...
                iblColors(i),'LineWidth',2); hold on
    end
    grid on;
    ylabel('$J^*_Q$','Interpreter','Latex','FontName','SansSerif');
    ax1 = gca;
    ax1.YColor = 'black';
    ax1.YLim = [0,20];
    ax1.FontSize = 16;
    ax1.FontName = 'SansSerif';
    xlabel('$z/\hat{\delta}$','Interpreter','Latex','FontName','SansSerif');
    xlim([0,5]);
    
end

figure();
for i = 1:length(delta_ibl)
   
     plot(z_global./delta_abl,ZXJq2{i},'-','Color',...
            iblColors(i),'LineWidth',2); hold on
    plot(z_global./delta_abl,ZXJq4{i},'--','Color',...
            iblColors(i),'LineWidth',2); hold on
    grid on;
    ylabel('$J^*_Q$','Interpreter','Latex','FontName','SansSerif');
    ax1 = gca;
    ax1.YColor = 'black';
    ax1.YLim = [0,20];
    ax1.FontSize = 16;
    ax1.FontName = 'SansSerif';
    xlabel('$z/\delta$','Interpreter','Latex','FontName','SansSerif');
    xlim([0,400/300]);
    
end

%% Plot Q2 and Q4 Separately

close all;

figure();
subplot(1,2,1);
for i = 1:length(delta_ibl)
   
    if i == 1 || i == 2
        plot(z_global./30,ZXJq2{i},'-','Color',...
            iblColors(i),'LineWidth',2); hold on
    else
        plot(z_global./(delta_ibl(i)),ZXJq2{i},'-','Color',...
                iblColors(i),'LineWidth',2); hold on
    end
end
grid on;
ylabel('$J^*_{II}$','Interpreter','Latex','FontName','SansSerif');
ax1 = gca;
ax1.YColor = 'black';
ax1.YLim = [0,20];
ax1.FontSize = 16;
ax1.FontName = 'SansSerif';
xlabel('$z/\hat{\delta}$','Interpreter','Latex','FontName','SansSerif');
xlim([0,3]);
subplot(1,2,2);
for i = 1:length(delta_ibl)
   
    if i == 1 || i == 2
        plot(z_global./30,ZXJq4{i},'--','Color',...
            iblColors(i),'LineWidth',2); hold on;
    else
        plot(z_global./(delta_ibl(i)),ZXJq4{i},'--','Color',...
                iblColors(i),'LineWidth',2); hold on
    end
end
grid on;
ylabel('$J^*_{IV}$','Interpreter','Latex','FontName','SansSerif');
ax1 = gca;
ax1.YColor = 'black';
ax1.YLim = [0,20];
ax1.FontSize = 16;
ax1.FontName = 'SansSerif';
xlabel('$z/\hat{\delta}$','Interpreter','Latex','FontName','SansSerif');
xlim([0,3]);
 

figure();
subplot(1,2,1);
for i = 1:length(delta_ibl)
   
    if i == 1 || i == 2
        plot(z_global./delta_abl,ZXJq2{i},'-','Color',...
            iblColors(i),'LineWidth',2); hold on
    else
        plot(z_global./delta_abl,ZXJq2{i},'-','Color',...
                iblColors(i),'LineWidth',2); hold on
    end
end
grid on;
ylabel('$J^*_{II}$','Interpreter','Latex','FontName','SansSerif');
ax1 = gca;
ax1.YColor = 'black';
ax1.YLim = [0,20];
ax1.FontSize = 16;
ax1.FontName = 'SansSerif';
xlabel('$z/\delta$','Interpreter','Latex','FontName','SansSerif');
xlim([0,400/300]);
subplot(1,2,2);
for i = 1:length(delta_ibl)
   
    if i == 1 || i == 2
        plot(z_global./delta_abl,ZXJq4{i},'--','Color',...
            iblColors(i),'LineWidth',2); hold on;
    else
        plot(z_global./delta_abl,ZXJq4{i},'--','Color',...
                iblColors(i),'LineWidth',2); hold on
    end
end
grid on;
ylabel('$J^*_{IV}$','Interpreter','Latex','FontName','SansSerif');
ax1 = gca;
ax1.YColor = 'black';
ax1.YLim = [0,20];
ax1.FontSize = 16;
ax1.FontName = 'SansSerif';
xlabel('$z/\delta$','Interpreter','Latex','FontName','SansSerif');
xlim([0,400/300]);

%% TQ

figure();
subplot(1,2,1);
for i = 1:length(delta_ibl)
   
    if i == 1 || i == 2
        plot(z_global./30,ZXTq2{i},'-','Color',...
            iblColors(i),'LineWidth',2); hold on
    else
        plot(z_global./(delta_ibl(i)),ZXTq2{i},'-','Color',...
                iblColors(i),'LineWidth',2); hold on
    end
end
grid on;
ax1 = gca;
ax1.YColor = 'black';
ax1.YLim = [0,12];
ax1.FontSize = 16;
ax1.FontName = 'SansSerif';
xlabel('$z/\hat{\delta}$','Interpreter','Latex',...
    'FontName','SansSerif','FontSize',36);
ylabel('$T^*_{II}$','Interpreter','Latex',...
    'FontName','SansSerif','FontSize',36);
xlim([0,3]);
subplot(1,2,2);
for i = 1:length(delta_ibl)
   
    if i == 1 || i == 2
        plot(z_global./30,ZXTq4{i},'--','Color',...
            iblColors(i),'LineWidth',2); hold on;
    else
        plot(z_global./(delta_ibl(i)),ZXTq4{i},'--','Color',...
                iblColors(i),'LineWidth',2); hold on
    end
end
grid on;
ax1 = gca;
ax1.YColor = 'black';
ax1.YLim = [0,12];
ax1.FontSize = 24;
ax1.FontName = 'SansSerif';
xlabel('$z/\hat{\delta}$','Interpreter','Latex',...
    'FontName','SansSerif','FontSize',36);
ylabel('$T^*_{IV}$','Interpreter','Latex',...
    'FontName','SansSerif','FontSize',36);
xlim([0,3]);
 

figure();
subplot(1,2,1);
for i = 1:length(delta_ibl)
   
    if i == 1 || i == 2
        plot(z_global./delta_abl,ZXTq2{i},'-','Color',...
            iblColors(i),'LineWidth',2); hold on
    else
        plot(z_global./delta_abl,ZXTq2{i},'-','Color',...
                iblColors(i),'LineWidth',2); hold on
    end
end
grid on;
ax1 = gca;
ax1.YColor = 'black';
ax1.YLim = [0,12];
ax1.FontSize = 24;
ax1.FontName = 'SansSerif';
xlabel('$z/\delta$','Interpreter','Latex','FontName','SansSerif',...
    'FontSize',36);
ylabel('$T^*_{II}$','Interpreter','Latex','FontName','SansSerif',...
    'FontSize',36);
xlim([0,400/300]);
subplot(1,2,2);
for i = 1:length(delta_ibl)
   
    if i == 1 || i == 2
        plot(z_global./delta_abl,ZXTq4{i},'--','Color',...
            iblColors(i),'LineWidth',2); hold on;
    else
        plot(z_global./delta_abl,ZXTq4{i},'--','Color',...
                iblColors(i),'LineWidth',2); hold on
    end
end
grid on;
ax1 = gca;
ax1.YColor = 'black';
ax1.YLim = [0,12];
ax1.FontSize = 24;
ax1.FontName = 'SansSerif';
xlabel('$z/\delta$','Interpreter','Latex','FontName','SansSerif',...
    'FontSize',36);
ylabel('$T^*_{IV}$','Interpreter','Latex','FontName','SansSerif',...
    'FontSize',36);
xlim([0,400/300]);

%%

close all;

iblColors = ["#000000", "#daf8e3", "#97ebdb", "#00c2c7", "#0086ad", "#005582", ...
    "#ffc100", "#ff9a00", "#ff7400", "#bf0000" ];

figure();
for i = 1:length(delta_ibl)
   
    plot(z_global./(delta_ibl(i)),ZXTq2{i},'-','Color',...
            iblColors(i),'LineWidth',2); hold on
    plot(z_global./(delta_ibl(i)),ZXTq4{i},'--','Color',...
            iblColors(i),'LineWidth',2); hold on
    grid on;
    ylabel('$T^*_q(z)$','Interpreter','Latex','FontName','SansSerif');
    ax1 = gca;
    ax1.YColor = 'black';
    ax1.YLim = [0,12];
    ax1.FontSize = 16;
    ax1.FontName = 'SansSerif';
    xlabel('$z/\delta_{IBL}$','Interpreter','Latex','FontName','SansSerif');
    xlim([0,5]);
    
end

figure();
for i = 1:length(delta_ibl)
   
     plot(z_global,ZXTq2{i},'-','Color',...
            myColors(i),'LineWidth',2); hold on
    plot(z_global,ZXTq4{i},'--','Color',...
            myColors(i),'LineWidth',2); hold on
    grid on;
    ylabel('$T^*_q(z)$','Interpreter','Latex','FontName','SansSerif');
    ax1 = gca;
    ax1.YColor = 'black';
    ax1.YLim = [0,12];
    ax1.FontSize = 16;
    ax1.FontName = 'SansSerif';
    xlabel('$z/\delta_{ABL}$','Interpreter','Latex','FontName','SansSerif');
    xlim([0,400]);
    
end


%% Data about the dunes at the center plane

% binned data by km 

xPeak_bin1 = [1950,2084,2194,2297,2374,2403,2464,2557,2721,2813,2935];
zPeak_bin1 = [3.939,6.717,7.428,9.926,5.916,5.698,7.428,6.875,8.823,17.3,12.48];
zPeak_bin1avg = mean(zPeak_bin1);
xPeak_bin1Spac = mean(diff(xPeak_bin1));

xPeak_bin2 = [3013,3086,3208,3360,3482,3574,3663,3892];
zPeak_bin2 = [9.813,12.08,11.69,15.2,12.42,12.64,14.33,12.51];
zPeak_bin2avg = mean(zPeak_bin2);
xPeak_bin2Spac = mean(diff(xPeak_bin2));

xPeak_bin3 = [4011,4102,4180,4331,4534,4771,4940];
zPeak_bin3 = [9.159,12.07,11.06,15.91,13.55,13.62,9.662];
zPeak_bin3avg = mean(zPeak_bin3);
xPeak_bin3Spac = mean(diff(xPeak_bin3));

xPeak_bin4 = [5076,5243,5335,5455,5545,5653,5743,5865];
zPeak_bin4 = [12.19,13.67,6.713,11.11,10.54,10.59,12.04,11.33];
zPeak_bin4avg = mean(zPeak_bin4);
xPeak_bin4Spac = mean(diff(xPeak_bin4));

xPeak_bin5 = [6108,6235,6352,6457,6518,6768];
zPeak_bin5 = [11.09,8.231,8.421,9.664,9.581,11.01];
zPeak_bin5avg = mean(zPeak_bin5);
xPeak_bin5Spac = mean(diff(xPeak_bin5));

xPeak_bin6 = [7018,7127,7200,7347,7431,7550];
zPeak_bin6 = [6.913,9.19,9.617,10.99,8.97,8.824];
zPeak_bin6avg = mean(zPeak_bin6);
xPeak_bin6Spac = mean(diff(xPeak_bin6));


%% END
