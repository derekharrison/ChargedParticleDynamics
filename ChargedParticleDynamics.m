%Charged particle dynamics

m_e = 1;           %mass of negatively charged particle in kg
max_v_e = 0.0;     %maximum component velocity for negatively charged part.
m_p = 1;            %mass of positively charged particle in kg
Np = 0;
Ne = 100;
k = 1000;
Qe = -1;
Qp = 1;
L = 50;
writevideo = true;
max_t = 7.5;
dt = 0.001;

%Initialization
err = 1e-8;
v_x_p = zeros(Np,1);
v_y_p = zeros(Np,1);

min_y = -L;
max_y = L;
min_x = -L;
max_x = L;

East_wall = 1;
West_wall = 2;
North_wall = 3;
South_wall = 4;

%Generate initial velocities for negatively charged particles
v_x_e = 2*max_v_e*rand(Ne,1) - max_v_e;
v_y_e = 2*max_v_e*rand(Ne,1) - max_v_e;

X_p = 2*L*rand(Np,1) - L;
Y_p = 2*L*rand(Np,1) - L;

X_e = 2*L*rand(Ne,1) - L;
Y_e = 2*L*rand(Ne,1) - L;

frameps = 96;
if writevideo==true
    writerObj = VideoWriter('C:\Users\d-w-h\Desktop\Home\Particle_dynamics.avi','Motion JPEG AVI');
    writerObj.FrameRate = frameps;
    open(writerObj);
end


t = 0;
frame_counter = 0;
while t < max_t
 
    %Check collision times with boundary
    coll_time = 1e+30;
    coll_partner_1 = 0;
    coll_partner_2 = 0;
    collision_with_p = false;
    
    %Checking collision time between wall and positively charged particle
    for n = 1:Np   
        coll_time_east = (L - X_p(n)) / v_x_p(n);
        if coll_time >= coll_time_east && coll_time_east > 0
            coll_time = coll_time_east;
            coll_partner_1 = n;
            coll_partner_2 = East_wall;
            collision_with_p = true;
        end
        
        coll_time_west = (-L - X_p(n)) / v_x_p(n);
        if coll_time >= coll_time_west && coll_time_west > 0
            coll_time = coll_time_west;
            coll_partner_1 = n;
            coll_partner_2 = West_wall;
            collision_with_p = true;
        end 
        
        coll_time_north = (L - Y_p(n)) / v_y_p(n);
        if coll_time >= coll_time_north && coll_time_north > 0
            coll_time = coll_time_north;
            coll_partner_1 = n;
            coll_partner_2 = North_wall;
            collision_with_p = true;            
        end  

        coll_time_south = (-L - Y_p(n)) / v_y_p(n);
        if coll_time >= coll_time_south && coll_time_south > 0
            coll_time = coll_time_south;
            coll_partner_1 = n;
            coll_partner_2 = South_wall;
            collision_with_p = true;            
        end            
    end
    
    %Checking collision time between wall and negatively charged particle
    for n = 1:Ne   
        coll_time_east = (L - X_e(n)) / v_x_e(n);
        if coll_time >= coll_time_east && coll_time_east > 0
            coll_time = coll_time_east;
            coll_partner_1 = n;
            coll_partner_2 = East_wall;
            collision_with_p = false;            
        end
        
        coll_time_west = (-L - X_e(n)) / v_x_e(n);
        if coll_time >= coll_time_west && coll_time_west > 0
            coll_time = coll_time_west;
            coll_partner_1 = n;
            coll_partner_2 = West_wall;
            collision_with_p = false;             
        end 
        
        coll_time_north = (L - Y_e(n)) / v_y_e(n);
        if coll_time >= coll_time_north && coll_time_north > 0
            coll_time = coll_time_north;
            coll_partner_1 = n;
            coll_partner_2 = North_wall;
            collision_with_p = false;             
        end  

        coll_time_south = (-L - Y_e(n)) / v_y_e(n);
        if coll_time >= coll_time_south && coll_time_south > 0
            coll_time = coll_time_south;
            coll_partner_1 = n;
            coll_partner_2 = South_wall;
            collision_with_p = false;             
        end            
    end    
    
    if coll_time < dt
        X_p = v_x_p * coll_time * (1 - err) + X_p;
        Y_p = v_y_p * coll_time * (1 - err) + Y_p;
        
        X_e = v_x_e * coll_time * (1 - err) + X_e;
        Y_e = v_y_e * coll_time * (1 - err) + Y_e;
        
        t = t + coll_time
                
        %Update velocities after collision with wall
        if collision_with_p == true
            if coll_partner_2 == East_wall || coll_partner_2 == West_wall
                v_x_p(coll_partner_1) = -v_x_p(coll_partner_1);
            elseif coll_partner_2 == North_wall || coll_partner_2 == South_wall
                v_y_p(coll_partner_1) = -v_y_p(coll_partner_1);
            end
        elseif collision_with_p == false
            if coll_partner_2 == East_wall || coll_partner_2 == West_wall
                v_x_e(coll_partner_1) = -v_x_e(coll_partner_1);
            elseif coll_partner_2 == North_wall || coll_partner_2 == South_wall
                v_y_e(coll_partner_1) = -v_y_e(coll_partner_1);
            end            
        end
        
        
        Fe = zeros(2, Ne);
        %Calculate force acting on negatively charged particles
        for n_e = 1:Ne
            %Calculating force on e particle due to p particles
            for n_p = 1:Np
                re = [X_e(n_e); Y_e(n_e)];
                rp = [X_p(n_p); Y_p(n_p)];
                r = re - rp;
                Fe(:,n_e) = Fe(:,n_e) + k * Qe * Qp * r / (r'*r * sqrt(r'*r) + err);
            end
            %Calculating force on e particles due to other e particles
            for n_e_n = 1:Ne
                if n_e_n ~= n_e
                    re = [X_e(n_e); Y_e(n_e)];
                    re_n = [X_e(n_e_n); Y_e(n_e_n)];
                    r = re - re_n;
                    Fe(:,n_e) = Fe(:,n_e) + k * Qe * Qe * r / (r'*r * sqrt(r'*r) + err);
                end
            end            
        end
        
        Fp = zeros(2, Np);
        %Calculate force acting on positively charged particles
        for n_p = 1:Np
            %Calculating force on e particle due to p particles
            for n_e = 1:Ne
                rp = [X_p(n_p); Y_p(n_p)];                
                re = [X_e(n_e); Y_e(n_e)];
                r = rp - re;
                Fp(:,n_p) = Fp(:,n_p) + k * Qe * Qp * r / (r'*r * sqrt(r'*r) + err);
            end
            %Calculating force on e particles due to other e particles
            for n_p_n = 1:Np
                if n_p_n ~= n_p
                    rp = [X_p(n_p); Y_p(n_p)];
                    rp_n = [X_p(n_p_n); Y_p(n_p_n)];
                    r = rp - rp_n;
                    Fp(:,n_p) = Fp(:,n_p) + k * Qp * Qp * r / (r'*r * sqrt(r'*r) + err);
                end
            end            
        end
        
        %Updating velocities of e particles
        for n_e = 1:Ne
            v_x_e(n_e) = Fe(1,n_e) / m_e * coll_time + v_x_e(n_e);
            v_y_e(n_e) = Fe(2,n_e) / m_e * coll_time + v_y_e(n_e);
        end
        
        %Updating velocities of p particles
        for n_p = 1:Np
            v_x_p(n_p) = Fp(1,n_p) / m_p * coll_time + v_x_p(n_p);
            v_y_p(n_p) = Fp(2,n_p) / m_p * coll_time + v_y_p(n_p);
        end        
        
    elseif coll_time >= dt
        X_p = v_x_p * dt + X_p;
        Y_p = v_y_p * dt + Y_p;
        
        X_e = v_x_e * dt + X_e;
        Y_e = v_y_e * dt + Y_e;
        
        t = t + dt
        
        Fe = zeros(2, Ne);
        %Calculate force acting on negatively charged particles
        for n_e = 1:Ne
            %Calculating force on e particle due to p particles
            for n_p = 1:Np
                re = [X_e(n_e); Y_e(n_e)];
                rp = [X_p(n_p); Y_p(n_p)];
                r = re - rp;
                Fe(:,n_e) = Fe(:,n_e) + k * Qe * Qp * r / (r'*r * sqrt(r'*r) + err);
            end
            %Calculating force on e particles due to other e particles
            for n_e_n = 1:Ne
                if n_e_n ~= n_e
                    re = [X_e(n_e); Y_e(n_e)];
                    re_n = [X_e(n_e_n); Y_e(n_e_n)];
                    r = re - re_n;
                    Fe(:,n_e) = Fe(:,n_e) + k * Qe * Qe * r / (r'*r * sqrt(r'*r) + err);
                end
            end            
        end
        
        Fp = zeros(2, Np);
        %Calculate force acting on positively charged particles
        for n_p = 1:Np
            %Calculating force on e particle due to p particles
            for n_e = 1:Ne
                rp = [X_p(n_p); Y_p(n_p)];                
                re = [X_e(n_e); Y_e(n_e)];
                r = rp - re;
                Fp(:,n_p) = Fp(:,n_p) + k * Qe * Qp * r / (r'*r * sqrt(r'*r) + err);
            end
            %Calculating force on e particles due to other e particles
            for n_p_n = 1:Np
                if n_p_n ~= n_p
                    rp = [X_p(n_p); Y_p(n_p)];
                    rp_n = [X_p(n_p_n); Y_p(n_p_n)];
                    r = rp - rp_n;
                    Fp(:,n_p) = Fp(:,n_p) + k * Qp * Qp * r / (r'*r * sqrt(r'*r) + err);
                end
            end            
        end

        %Updating velocities of e particles
        for n_e = 1:Ne
            v_x_e(n_e) = Fe(1,n_e) / m_e * dt + v_x_e(n_e);
            v_y_e(n_e) = Fe(2,n_e) / m_e * dt + v_y_e(n_e);
        end
        
        %Updating velocities of p particles
        for n_p = 1:Np
            v_x_p(n_p) = Fp(1,n_p) / m_p * dt + v_x_p(n_p);
            v_y_p(n_p) = Fp(2,n_p) / m_p * dt + v_y_p(n_p);
        end         
    end
    
    dummy = floor(t/(2*dt));
    if frame_counter == dummy
        frame_counter = frame_counter + 1;
        plot(X_e, Y_e, 'r.', 'MarkerSize', 10)
        hold on
        plot(X_p, Y_p, 'b.', 'MarkerSize', 10)
        hold off

        set(gca,'Ylim',[min_y max_y])
        set(gca,'Xlim',[min_x max_x])
        set(gca,'nextplot','replacechildren');
        set(gcf,'Renderer','zbuffer')
        frame = getframe(gcf); 
        if writevideo == true      
            writeVideo(writerObj,frame);
        end
    end
end

if writevideo == true
    close(writerObj);
end