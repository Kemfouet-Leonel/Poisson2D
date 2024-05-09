%vortex_rolls.m Free convection problem between to flat plates
%
% Author: Chloros2 <chloros2@gmx.de>
% Created: 2018-05-13
%
% Copyright (C) 2018 Chloros2 <chloros2@gmx.de>
%
% This program is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see
% <https://www.gnu.org/licenses/>.
%

clear all;
clc

addpath('ffmatlib');

%Reads a FreeFem++ mesh created with the savemesh(Th,"mesh.msh"); command
[p,e,t]=ffreadmesh('poisson2D.msh');
[Vh]=ffreaddata('poisson2D_Vh.txt');
[u]=ffreaddata('poisson2D_u.txt');
[w]=ffreaddata('poisson2D_w.txt');
[h]=ffreaddata('poisson2D_h.txt');
[err]=ffreaddata('poisson2D_err.txt');

% ffpdeplot(p,e,t, 'Mesh','on','XLim',[0 1],'YLim',[0 1]);
% figure;  

 
ffpdeplot(p,e,t, ...
         'VhSeq',Vh, 'Mesh','off', ...
         'Boundary','off', ...
         'XYData',u,'Zstyle','continuous', ...
         'ColorMap',jet(200),'ColorBar','on', ...
         'XLim',[0 1],'YLim',[0 1],'ZLim',[0 0.07] );         
 figure;
 
ffpdeplot(p,e,t, ...
         'VhSeq',Vh, 'Mesh','off', ...
         'Boundary','off', ...
         'XYData',w,'Zstyle','continuous',...
         'ColorMap',jet(200),...
         'ColorBar','on','XLim',[0 1],'YLim',[0 1],'ZLim',[0 0.07]);  
% lighting gouraud;
% view([-47,24]);
% camlight('headlight');
  %   figure; 
   
%ffpdeplot(p,e,t,...%     a = length(err)
      np=size(p,2);
      nt=size(t,2);
  %  h=0:1/(np-1):1;
    length(h)
    length(err)
    plf=polyfit(log(h),log(err),2)
figure;
    plot(log(h), log(err),'r', log(h), plf(1)*log(h),'k-')
%         'VhSeq',Vh, 'Mesh','off',...
%         'Boundary','off',...
%         'XYData',err,'Zstyle','continuous',...
%          'ColorMap','parula',...
%         'ColorBar','on','XLim',[0 1],'YLim',[0 1],'ZLim',[0 1] );   
    
