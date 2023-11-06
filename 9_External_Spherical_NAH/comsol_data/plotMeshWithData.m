function plotMeshWithData(ax, nodes, triangles, data, colorLimits, plotcolorbar, viewAngles, axesLimits)
    % Activate the specified axes handle
    axes(ax);
    
    % Enable anti-aliasing
    % ax.GraphicsSmoothing = 'on';
    
    % Plot the mesh without edge lines
    trisurf(triangles, nodes(:,1), nodes(:,2), nodes(:,3), data,'EdgeColor', 'none', 'FaceColor','interp');
    
    % Set colormap
    colormap jet;
    
    % Set color limits
    caxis(colorLimits);
    
    % Add colorbar
    if plotcolorbar
        colorbar;
    end
    
    % Set view angles
    view(viewAngles);
    
    % Set axes limits
    xlim(axesLimits(1,:));
    ylim(axesLimits(2,:));
    zlim(axesLimits(3,:));
    
    % Set axis labels
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    
    % Set title
    title('Mesh with Colored Data');
end




% function ax = plotMeshWithData(nodes, triangles, data, colorLimits, viewAngles, axesLimits)
%     % Create figure and axes
%     fig = figure;
%     ax = axes(fig);
%     
%     % Enable anti-aliasing
%     fig.GraphicsSmoothing = 'on';
% 
%     % Plot the mesh
%     trisurf(triangles, nodes(:, 1), nodes(:, 2), nodes(:, 3), data,'EdgeColor', 'none');
%     
%     % Set colormap
%     colormap jet;
%     % Set color limits
%     caxis(colorLimits);
%     % Add colorbar
%     colorbar;
%     % Set view angles
%     view(viewAngles);
%     
%     % Set axes limits
%     xlim(axesLimits(1,:));
%     ylim(axesLimits(2,:));
%     zlim(axesLimits(3,:));
%     
%     % Set axis labels
%     xlabel('X');
%     ylabel('Y');
%     zlabel('Z');
%     
%     % Set title
% %     title('Mesh with Colored Data');
% end
