function plot_tracks_scaled(coordinates)
    fig = figure('Position', [130 10 800 800]);
    f_n = fieldnames(coordinates);
    button = '';
    Previous = uicontrol('String','Previous','Position',[30 10 90 30]);
    Next = uicontrol('String','Next','Position',[130 10 90 30]);
    f = 1;
    i = 1;
    while f <= length(f_n)
        while i <= size(coordinates.(f_n{f}).original_x,2)
            Previous.Callback = {@ButtonMode,fig};
            Next.Callback = {@ButtonMode,fig};
            plot(coordinates.(f_n{f}).original_x(:,i), ...
                coordinates.(f_n{f}).original_y(:,i),'Color','b')
            hold on
            plot(coordinates.(f_n{f}).original_x(end,i), ...
                coordinates.(f_n{f}).original_y(end,i),...
                'ko','MarkerFaceColor','k','MarkerSize', 1.75)
            text(.05,.95,[f_n{f},'\_nº',num2str(i)],'Units','normalized')
            axis equal
%             if MaxX > MaxY % Uncomment this to make square plots
%                 axis([-MaxX MaxX -MaxY*(MaxX/MaxY) MaxY*(MaxX/MaxY)]);
%             elseif MaxY > MaxX
%                 axis([-MaxX*(MaxY/MaxX) MaxX*(MaxY/MaxX) -MaxY MaxY]);
%             elseif MaxY == MaxX
%                 axis([-MaxX MaxX -MaxY MaxY]);
%             end
            hold off
            uiwait(fig)
            if contains(button,'Next')
                i = i + 1;
            elseif contains(button,'Previous')
                if (i ~= 1)
                    i = i - 1;
                elseif (i == 1) && (f == 1)
                    msgbox('First entry of the dataset has been reached')
                elseif (i == 1) && (f ~= 1)
                    f = f - 1;
                    i = size(coordinates.(f_n{f}).original_x,2) ;
                end
            end
        end
        f = f + 1;
        i = 1;
    end
    function  ButtonMode(object,~,fig)
        if strcmpi(object.String, 'Previous')
            button = 'Previous';
            uiresume(fig)
        elseif strcmpi(object.String, 'Next')
            button = 'Next';
            uiresume(fig)
        end
    end
end