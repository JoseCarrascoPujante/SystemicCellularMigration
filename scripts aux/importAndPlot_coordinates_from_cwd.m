files = dir('*.xlsx') ;
files = natsortfiles(files) ;
N = length(files) ;
figure
hTracks = gca;

% Custom pixel/mm ratio list for Chemotaxis of Metamoeba leningradensis
ratio_list = [23.44,23.44,23.44,23.44,23.44,23.44,23.44,11.63,11.63,11.63,...
    11.63,11.63,11.63,11.63,11.63,23.44,23.44,23.44,23.44,23.44,23.44,...
    23.44,23.44,23.44,23.44,23.44,23.44,23.44,23.44,23.44,23.44,23.44,...
    23.44,23.44,23.44,23.44,11.63,11.63,11.63,11.63,23.44,23.44,23.44...
    ,11.63,11.63,11.63,11.63,11.63,11.63,11.63,23.44,23.44,38.7651,38.7651,...
    38.7651,38.7651,38.7651,38.7651,38.7651,38.7651,38.7651,38.7651,38.7651,...
    38.7651,38.7651,38.7651,38.7651,38.7651,38.7651,38.7651,38.7651] ;

for i = 30:36
    thisfile = files(i).name ;
    temp_x = readmatrix(thisfile,'Range','E1:E3600') ;
    temp_y = readmatrix(thisfile,'Range','F1:F3600') ;
    % testx = [zeros(length(temp_x), 1)] ;
    % testy = [zeros(length(temp_y), 1)] ;
    % Two separate loops in case x and y have different lengths
%     if mod(temp_x,1) == 0
%         for k1 = 1:length(temp_x)
%             if temp_x(k1) >= 1000000
%                 msgbox('Este valor de X es mayor que 1000, problem')
%             end
%             temp_x(k1) = str2double(separatethousands(temp_x(k1),'.',0)) ;
%         end
%     end
%     if mod(temp_y,1) == 0
%         for k1 = 1:length(temp_y)
%             if temp_y(k1) >= 1000000
%                 msgbox('Este valor de Y es mayor que 1000, problem')
%             end
%             temp_y(k1) = str2double(separatethousands(temp_y(k1),'.',0)) ;
%         end
%     end
%     eval(['x' num2str(i) '=temp_x']);
%     eval(['y' num2str(i) '=temp_y']);
    
    % center
    centered_x = temp_x-temp_x(1) ;
    centered_y = temp_y-temp_y(1) ;
    
    % convert to polar coordinates
%     [teta,rho] = cart2pol(centered_x,centered_y);
%     rho = rho/11.63;
    
    % scale
    scaled_x = centered_x/ratio_list(i);
    scaled_y = centered_y/ratio_list(i);
    
    % Plot trajectory and 'ko' marker at its tip
    disp(i)
    plot(hTracks, scaled_x, scaled_y, 'Color', [0, 0, 0]) ;
    hold on;
    plot(hTracks, scaled_x(end), scaled_y(end), 'ko',  'MarkerFaceColor',  [0, 0, 0], 'MarkerSize', 2) ;
end

%Adjust proportions of 'tracks' figure after each condition   
divx=[-18 18];
divy=[0 0];
plot(hTracks,divx,divy,'k');
hold on;
plot(hTracks,divy,divx,'k');
hold on;
axis([-18 18 -18 18]);
daspect([1 1 1]);
box on;
hold off
% [name0,name2,name3]=fileparts(pwd) ;
% [~,name1] = fileparts(name0) ;
% filename = strcat(name1, 32, name2, name3, '.mat') ;
% clear name0 name1 name2 name3 temp_y temp_x files thisfile N i
% save(filename) ;