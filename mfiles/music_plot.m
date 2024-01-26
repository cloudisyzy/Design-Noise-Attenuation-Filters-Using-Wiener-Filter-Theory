function   music_plot(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10)

%
%   music_plot(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10)
%
%  music_plot:
%     Make a plot with note names on the Y axis. 
%     Same arguments as ordinary plot
%     
%     Mats Bengtsson,  9/2 1998
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This can be simplified in Matlab 5, but some students may use
% Matlab 4.
switch nargin
  case 10,
    plot(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10);
  case 9, 
    plot(p1, p2, p3, p4, p5, p6, p7, p8, p9);
  case 8,
    plot(p1, p2, p3, p4, p5, p6, p7, p8);
  case 7,
    plot(p1, p2, p3, p4, p5, p6, p7);
  case 6,
    plot(p1, p2, p3, p4, p5, p6);
  case 5,
    plot(p1, p2, p3, p4, p5);
  case 4,
    plot(p1, p2, p3, p4);
  case 3,
    plot(p1, p2, p3);
  case 2,
    plot(p1, p2);
  case 1,
    plot(p1);
end
names=['a '; 'a#';'h ';'c ';'c#';'d ';'d#';'e ';'f ';'f#';'g ';'g#'];
ylim=get(gca,'ylim');
ymin=max(-10,ylim(1));
ymax=min(10, ylim(2));
set(gca, 'ylim', [ymin,ymax]);
set(gca, 'ytick', [ymin:ymax]);
set(gca, 'yticklabel', names(1+mod([ymin:ymax],12), :));
