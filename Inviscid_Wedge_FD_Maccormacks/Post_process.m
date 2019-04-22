%Script to post-process results from wedge finite-difference code

num_y_dir_nodes = 101;      %Number of nodes chosen to use in y-dir

%Import data [iteration|x|y|u|v|M|p|T|rho]
flow_field = importdata('flow_field.txt');

%Extract x and y-coordinates
x_coord = flow_field.data(:,2);
y_coord = flow_field.data(:,3);

%Determine domain
x_max = max(x_coord);
y_max = max(y_coord);

%Extract u,v,M,p,T,rho
u = flow_field.data(:,4);
v = flow_field.data(:,5);
M = flow_field.data(:,6);
p = flow_field.data(:,7);
T = flow_field.data(:,8);
rho = flow_field.data(:,9);

%Plot Mach, pressure, Temp, and density contour 
PlotContour(x_coord,y_coord,M,x_max,y_max,num_y_dir_nodes,'Mach Contour, 5 Degree Wedge at M_{\infty} = 2')
PlotContour(x_coord,y_coord,p,x_max,y_max,num_y_dir_nodes,'Pressure Contour, 5 Degree Wedge at M_{\infty} = 2')
PlotContour(x_coord,y_coord,T,x_max,y_max,num_y_dir_nodes,'Temperature Contour, 5 Degree Wedge at M_{\infty} = 2')
PlotContour(x_coord,y_coord,rho,x_max,y_max,num_y_dir_nodes,'Density Contour, 5 Degree Wedge at M_{\infty} = 2')

%Plot velocity vector field
PlotVelocityField(x_coord,y_coord,u,v,x_max,y_max,num_y_dir_nodes,'Velocity Vector Field')

%Function to plot velocity vector field
function PlotVelocityField(x_coord,y_coord,u,v,x_max,y_max,num_y_dir_nodes,plot_title)
    plot_x_loc = 150;
    plot_y_loc = 150;
    width = 800;
    height = 500;
    n = 20;
    num_x_samp = floor((length(x_coord) / num_y_dir_nodes)/n);
    num_y_samp = floor(num_y_dir_nodes/n);

    %Sample every 'n' point for plotting to reduce cluttering in plot
    x_samp = x_coord(1:n:end);
    y_samp = y_coord(1:n:end);
    u_samp = u(1:n:end);
    v_samp = v(1:n:end);

    figure
    quiver(x_samp,y_samp,u_samp,v_samp,2.5)
    title(plot_title)
    xlim([0,x_max])
    ylim([0,y_max])
    xlabel('x-location (m)')
    ylabel('y-location (m)')
    set(gcf,'position',[plot_x_loc,plot_y_loc,width,height])    
    
end

%Function to plot contour of selected component
function PlotContour(x_coord, y_coord, profile, x_max, y_max,num_y_dir_nodes,plot_title);
    X = reshape(x_coord,[num_y_dir_nodes,(length(x_coord)/num_y_dir_nodes)]);
    Y = reshape(y_coord,[num_y_dir_nodes,(length(x_coord)/num_y_dir_nodes)]);
    Variable = reshape(profile,[num_y_dir_nodes,(length(x_coord)/num_y_dir_nodes)]);
    
    plot_x_loc = 150;
    plot_y_loc = 150;
    width = 800;
    height = 500;
    
    figure
    contourf(X,Y,Variable)
    title(plot_title)
    xlim([0,x_max])
    ylim([0,y_max])
    colorbar
    xlabel('x-location (m)')
    ylabel('y-location (m)')
    set(gcf,'position',[plot_x_loc,plot_y_loc,width,height])
end