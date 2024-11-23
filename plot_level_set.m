function plot_level_set(t, y, plotLevelSet, beforeCorrect, afterCorrect)

    global raw_model

    if plotLevelSet && ~raw_model
        if beforeCorrect
            day = (t(1)/24)+1;
            minx_idx = find(y(:,1) == min(y(:,1)));
            figure();
            hold on;
            xlabel('x', 'FontSize', 13);
            ylabel('xc', 'FontSize', 13);
            title('Day ' + string(day) + ' Level Set Ellipsoid on x,xc Plane',...
                    'FontSize', 13);
            axis([-1.5 1.5, -1.5 1.5]);
            
            for i = 1:30:minx_idx
                xM = reshape(y(i,:),[3,4]);
                x = xM(:,1);
                M = xM(:,2:4);
                opposite_vecs = -1.*M(1:2,1:2);
                span_vectors = [M(1:2,1:2)+x(1:2), opposite_vecs+x(1:2)];
                [A,c] = MinVolEllipse(span_vectors, 0.0001);
                Ellipse_plot(A,c,beforeCorrect,afterCorrect);
                hold on;
            end
%             hold on;
            
        elseif afterCorrect
            for i = 1:30:numel(t)
                xM = reshape(y(i,:),[3,4]);
                x = xM(:,1);
                M = xM(:,2:4);
                opposite_vecs = -1.*M(1:2,1:2);
                span_vectors = [M(1:2,1:2)+x(1:2), opposite_vecs+x(1:2)];
                [A,c] = MinVolEllipse(span_vectors, 0.0001);
                Ellipse_plot(A,c,beforeCorrect,afterCorrect);
%                 hold on;
            end
            hold off;
        end
    end
    
    raw_model = 0;

%     x_result = xM_result(:,1);
%     M_result = xM_result(:,2:4);
% 
%     h = figure;
%     axis tight equal
%     filename = 'LSKF_Ellipsoids.gif';
%     xlabel('x', 'FontSize', 14);
%     ylabel('xc', 'FontSize', 14);
%     title('Level Set Ellipsoid on x,xc Plane', 'FontSize', 14);
%     xlim([-1.25,1.25]);
%     ylim([-1.25,1.25]);
%     hold on;
%     
%     for ij = 1:10:numel(time_points)-1
%         opposite_vecs = -1.*M(1:2,1:2,ij);
%         span_vectors = [M(1:2,1:2,ij)+x(1:2,ij), opposite_vecs+x(1:2,ij)];
%         [A,c] = MinVolEllipse(span_vectors, 0.0001);
%         Ellipse_plot(A,c,isPhase);
%         drawnow;
%     
%         frame = getframe(h);
%         im = frame2im(frame);
%         [imind,cm] = rgb2ind(im,256);
%         
%         if ij == 1 
%           imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
%         else 
%           imwrite(imind,cm,filename,'gif','WriteMode','append'); 
%         end 
%         clf;
%     end
%     hold off;

end