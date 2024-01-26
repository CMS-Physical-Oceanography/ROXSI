function  bgpatch2(x,ym,yst,pc,types)
% bgpatch2(x,ym,yst,pc,types)
%       Function to plot envelope of STD around a series of means.
%
% x = x-axis (y-axis for types =2)
% ym = mean  
% yst = standard dev
% pc = 'r' red, 'b' blue, 'k1' black1, 'k2' black2, 'g' green, 'c' cyan
%
% types = 1 or 2 (1 is original (default), 2 is envelope in y) 3/22/11
%
% Difference with bgpatch is that we don't call 'FaceAlpha' to give
% transparency. Weird, but whatever, now can save as .eps.

if pc == 'r' %red
    C = [1 .8 .8];
elseif pc == 'b' %blue
    C = [.8 .8 1]; 
elseif pc == 'k' %black 1
    C = [.8 .8 .8];
elseif pc == 'k1' %black 1
    C = [.8 .8 .8]; 
elseif pc == 'k2' %black 2
    C = [.6 .6 .6]; 
elseif pc == 'g' %green
    C = [.8 1 .8]; 
elseif pc == 'c' %cyan
    C = [.8 1 1]; 
end

% For whether the envelope is in x or in y.
if exist('types','var') == 1
    if types == 1
    hold on
    for k = 1:length(x)-1
    patch([x(k) x(k) x(k+1) x(k+1)],[ym(k)+yst(k) ym(k)-yst(k) ym(k+1)-yst(k+1) ym(k+1)+yst(k+1)],C,'linestyle','none','HandleVisibility','off','facealpha',0.5)
    %[topleft bottomleft bottomright topright] 
    end
    hold off
    else
        xm = ym;
        y = x;
        xst = yst;
        for k = 1:length(xm)-1
        patch([xm(k+1)-xst(k+1) xm(k)-xst(k) xm(k)+xst(k) xm(k+1)+xst(k+1)],[y(k+1) y(k) y(k) y(k+1)],C,'linestyle','none','HandleVisibility','off','facealpha',0.5)
        end
        hold off
    end
else % back to original
hold on
for k = 1:length(x)-1
patch([x(k) x(k) x(k+1) x(k+1)],[ym(k)+yst(k) ym(k)-yst(k) ym(k+1)-yst(k+1) ym(k+1)+yst(k+1)],C,'linestyle','none','HandleVisibility','off','facealpha',0.5) 
end
hold off 
end
set(gca,'Layer','top');
        