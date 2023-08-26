function [aveAGI,dev] = secantmethode(structure, v)
%% Metoda siecznych 3D final version dir X
% Author: Piotr ¯urek
% Date of cerate: 01.10.2018
% Up-date: 06.12.209
% Scrpt calculate the anisotropic avarage grain size in X direction
% Based on AGI (Avarage Grain Intersept) methode.
% The method uses a statistical approach.
% The algorithm determines a large number of secants
% Secants intersect the three-dimensional structure
% In place where voxel change value from 0 on 255 and from 255 on 0,
% the algorithm mark a point. Next add secants
% length and then divides this value by the number of all points.
% In results gets AGI coefficient.
% Function downloads folder path and voxel size
% Folder have to include set of images save in .bmp
% Voxel must have the same size in each direction
% Function return grain size and standard deviation

%% Import image

% Function fileparts cuts folder path to 3 part:
% main path, folder name and file extension.
% The path must be entered in this form, e.g.:
% C:\Program Files\Another folders\Folder with images.bmp
%[filepath,name,ext] = fileparts(foldername);
% Count the number of items and saved item names
%path_image = dir(num2str([filepath,'\',name,'\*',ext]));
% Redefining an empty cell matrix equal to the number of images
%Image = cell(1,length(path_image));
% Loop saves bitmap in cells Image array. One image in one cell
%for i = 1:length(path_image)    
 %   filename = strcat(num2str([filepath,'\',name,'\']),path_image(i).name);
  %  Image{i} = imread(filename);  
%end

    % Makes 3D structure from 2D images
    stack = structure;
    % Calculate demonsion 3D structure
    size_stack = size(stack);
              
%% Inputs

% Number of secants
% You can change it as you like (Don't go below 2000)
s = 2000; 
% Randomization of points within the 3D structure.
% Points are the basis for sketching secants

    %if direction == 'x'
        A = [ones(s,1),...
            randi([5,size_stack(1)-5],s,1),...
            randi([5,size_stack(3)-5],s,1)];    
        B = [size_stack(1)*ones(s,1),...
            A(:,2),...
            A(:,3)];
%{        
    elseif  direction == 'y'
        A = [randi([5,size_stack(1)-5],s,1),...
            ones(s,1),...
            randi([5,size_stack(3)-5],s,1)];    
        B = [A(:,1),...
            size_stack(2)*ones(s,1),...
            A(:,3)]; 
    elseif  direction == 'z'
        A = [ones(s,1),...
            randi([5,size_stack(1)-5],s,1),...
            randi([5,size_stack(3)-5],s,1)];    
        B = [size_stack(1)*ones(s,1),...
            A(:,2),...
            A(:,3)];
    elseif  direction == 'none'
       A = [randi([1,size_stack(1)],s,1),...
            randi([1,size_stack(2)],s,1),...
            randi([1,size_stack(3)],s,1)];    
        B = [randi([1,size_stack(1)],s,1),...
            randi([1,size_stack(2)],s,1),...
            randi([1,size_stack(3)],s,1)];
end;
%}
%% Designation of secants

% Redefining an empty cell matrix: X,Y,Z,R.
 X(1:s) = {0};
 Y(1:s) = {0};
 Z(1:s) = {0};
 R(1:s) = {0};
% Bresenham function make a digital line between two points
for i=1:s   
    [X{i}, Y{i}, Z{i}] = bresenham_line3d(A(i,:,:), B(i,:,:));
end

 for i=1:s
    for j=1:length(X{1,i})
        % Cell array R contains information about which voxels 
        % go through each line 
        R{i}(j)=stack(X{1,i}(j), Y{1,i}(j), Z{1,i}(j));
    end 
 end
% Usuwa linie zawieraj¹ce tylko woksele czarne lub bia³e
for i=1:s
    if (R{i}(1,1)==255||R{i}(1,end)==255||...
            sum(R{i})==0||sum(R{i})==length(R{i})*255)
        R{i} = [];
    end
end
empties = cellfun(@isempty,R);
R(empties) = [];
% the loop counts the intersection points on each line
% Redefining an empty cell matrix to save the number of intersection
zm = zeros(length(R),1);
% Redefining an empty cell matrix to save the number of white voxels on the line
count = zeros(length(R),1);
  for i=1:length(R)
    for j=1:length(R{1,i})-1
        if R{i}(1,j)~=R{i}(1,j+1)
            zm(i) = zm(i)+1;
        end
            if R{i}(1,j)==255
                count(i)=count(i)+1;
            end
    end
  end
% Delete lines where algorithm counts only one intercept
for i=1:length(R)
    if zm(i)<2
         R{i} = [];
         zm(i) = inf;
         count(i) = inf;
    end
end
empties = cellfun(@isempty,R);
R(empties) = [];
zm = zm(zm~=inf);
count = count(count~=inf);

%% Calculation AGI and standard deviation
% Constant 1.56 gets from:
% Mendelson M. I. Average Grain Size in Polycrystalline Ceramics. 
% Journal of The American Ceramic Society, 1969.
AGI = 4.5*((v*count)./(0.5*zm));
aveAGI = mean(AGI);
dev = std(AGI);

%% End of script
end