function draw_structure(obj, Ls, Rs, thetas)

    display('    preparing to draw structure ...');
    N_per_split = obj.N_per_split;
    
    theta_calc = @(h,R) (h < 2*R).*acos(1-h./2./R) + (h >= 2*R).*pi/2;
    d_calc     = @(h,R) 2.*R.*sin(theta_calc(h,R));
    L_calc     = @(h,R) (h < 2*R).*2.*R.*theta_calc(h,R) + (h >= 2*R).*(pi.*R+h-2.*R);
    diff_calc  = @(h,R)  L_calc(h,R)-d_calc(h,R);
    
    Nsplits = length(Ls);
    Nout = N_per_split^Nsplits;
    
    wg_labels = {};
    for i = (0:Nout-1)
        wg_labels{i+1} = dec2bin(i,Nsplits);
    end

    xpos_list = [0];
    ypos_list = [0];
    count = 1;
    total_lengths = zeros(Nout,1);
    for i = (0:Nout-1)
        label_str = wg_labels{i+1};
        L_count = 0.0;
        pos = [0 0];
        for j = (1:Nsplits)
            h = obj.L/N_per_split^j;  
            R = Rs(j);
            if strcmp(label_str(j), '1')
                L_count = L_count + L_calc(h,R);
                pos(2) = pos(2) + h;
            else
                L_count = L_count + d_calc(obj.L/N_per_split^j,Rs(j));
            end       
            pos(1) = pos(1) + d_calc(h,R);
            xpos_list = [xpos_list, pos(1)];
            ypos_list = [ypos_list, pos(2)];        
            count = count + 1;
        end
        total_lengths(i+1) = L_count;

    end

    length_total = max(xpos_list);
    height_total = max(ypos_list);
    dl = 5e-8;
    thick = 30*obj.beta;
    wg_arm_length = 50e-6;
    half_thick = round(thick/2);
    spc = 400;
    Nx = ceil(length_total/dl) + 4*spc;
    Ny = ceil(height_total/dl) + 2*spc;

    E = zeros(Nx,Ny);

    upd = textprogressbar(Nout, 'barlength', 20, ...
                         'updatestep', 5, ...
                         'startmsg', '        ',...
                         'endmsg', ' Done! ', ...
                         'showbar', true, ...
                         'showremtime', true, ...
                         'showactualnum', true, ...
                         'barsymbol', '=', ...
                         'emptybarsymbol', ' ');

    for i = (0:Nout-1)
        upd(i+1)

        label_str = wg_labels{i+1};

        posx = spc;
        posy = spc;

        for j = (1:Nsplits)
            h = (obj.L/N_per_split^j);
            r = Rs(j);
            dist = d_calc(h,r);

            H = ceil(h/dl);  
            R = ceil(r/dl);
            R_out = R + half_thick;
            R_in = max(R - half_thick,1);

            D = ceil(dist/dl);

            angle = thetas(j);
            Nt = 1000;
            theta_range = linspace(0,angle,Nt);

            if strcmp(label_str(j), '1')
                circle_center1 = [posx, posy + R];
                circle_center2 = [posx + D, posy + H - R];
                for ti = (1:Nt)
                    t = theta_range(ti);
                    pos_circle1 = circle_center1 + ceil([R*sin(t), - R*cos(t)]);
                    pos_circle2 = circle_center2 + ceil([-R*sin(t),  R*cos(t)]);


                    %pos_circle1_out = circle_center1 + ceil([R_out*sin(t), - R_out*cos(t)]);
                    %pos_circle1_in =  circle_center1 + ceil([R_in*sin(t),  - R_in*cos(t)]);
                    %pos_circle2_out = circle_center2 + ceil([-R_out*sin(t), R_out*cos(t)]);
                    %pos_circle2_in =  circle_center2 + ceil([-R_in*sin(t),  R_in*cos(t)]);

                    %E(pos_circle1_in(1):pos_circle1_out(1),pos_circle1_out(2):pos_circle1_in(2)) = 1;
                    %E(pos_circle2_in(1):pos_circle2_out(1),pos_circle2_out(2):pos_circle2_in(2)) = 1;

                    E(pos_circle1(1):pos_circle1(1)+thick,pos_circle1(2):pos_circle1(2)+thick) = 1;
                    E(pos_circle2(1):pos_circle2(1)+thick,pos_circle2(2):pos_circle2(2)+thick) = 1;                
                end

                E(posx+ceil(dist/2/dl):posx+ceil(dist/2/dl)+thick, posy+R:posy+H-R) = 1;            
                posy = posy + H;
            else
                E(posx:posx+D+thick, posy:posy+thick) = 1;
            end 

            posx = posx + D;  
            if (j == Nsplits)
                E(posx:posx+round(wg_arm_length/dl), posy:posy+thick) = 1;
            end
        end
    end

    %%
    figure(); clf;
    imagesc((1:Nx)*dl*1e6, (1:Ny)*dl*1e6, E')
    %axis equal tight
    set(gca,'YDir','normal')
    colormap(flipud(gray))
    xlabel('x position (um)')
    ylabel('y position (um)')
    set(findall(gcf,'type','text'),'FontSize',22,'fontWeight','bold')
    set(gca,'FontSize',22,'fontWeight','bold')

end