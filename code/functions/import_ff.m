function [theta, phi, Et, Ep, ff_info] = import_ff(source, type)
%
% Read coordinates and E-fields from CST or Wavestudio far-field files. 
% Only point frequencies in CST!
% ------------------------------------------------------------------------
% INPUT  source : string, path to directory OR full filename
%        type   : string, "cst" (default) or "wavestudio" 
%
% OUTPUT theta  : (:,1), theta values in radians
%        phi    : (:,1), phi values in radians
%        Et     : (:,N), E-field theta component, all elements
%        Ep     : (:,N), E-field phi component, all elements
%        ff_info: struct, frequencies and other information
% ------------------------------------------------------------------------
% 06.06.2024 Albert Salmi, Department of Electronics and Nanoengineering,
%                          Aalto University School of Electrical
%                          Engineering
% ------------------------------------------------------------------------
%

arguments
    source (:,:)
    type (:,:)      = "cst"
end

% Check source type and define suitable function handle
if type == "cst"
    importfun = @(fname) import_CST_ff(fname);
    ext = "*.ffs";
elseif type == "wavestudio"
    importfun = @(fname) import_wavestudio_ff(fname);
    ext = "*.txt";
else
    error('Type must be either wavestudio or cst');
end

if ~isfolder(source) % read single file

    [theta, phi, Et, Ep, ff_info] = importfun(source);

else % read all files from directory

    % Read directory
    fol = dir(fullfile(source, ext));
    N = length(fol);

    % Read first file and initialize E-matrices with correct size
    [theta, phi, Et1, Ep1, ff_info1] = importfun(fullfile(fol(1).folder, fol(1).name));
    Et = zeros(length(Et1), N);
    Ep = zeros(length(Ep1), N);
    ff_info_us = cell(1,N); % us = unsorted
    Et(:,1) = Et1;
    Ep(:,1) = Ep1;
    ff_info_us{1} = ff_info1;
    
    % Read other files, assume that theta and phi are same
    for it = 2:N
        [~, ~, Et(:,it), Ep(:,it), ff_info_us{it}] = importfun(fullfile(fol(it).folder, fol(it).name));
    end
    
    % Sort data so that element 1 is at first column
    elements = zeros(N,1);
    for it = 1:N
        elements(it) = ff_info_us{it}.element;
    end
    [~, element_sort_idx] = sort(elements);
    Et = Et(:,element_sort_idx);
    Ep = Ep(:,element_sort_idx);
    ff_info = cell(N,1);
    for it = 1:N
        ff_info{it} = ff_info_us{element_sort_idx(it)};
    end

end

end


function [theta, phi, Et, Ep, ff_info] = import_CST_ff(fname)
%
% Read coordinates and E-fields from CST far-field files. Only point
% frequencies!
% ------------------------------------------------------------------------
% INPUT  fname  : string, filename and path in CST default format
%
% OUTPUT theta  : (:,1), theta values in radians
%        phi    : (:,1), phi values in radians
%        Et     : (:,1), E-field theta component
%        Ep     : (:,1), E-field phi component
%        ff_info: struct, frequencies and other information
% ------------------------------------------------------------------------
% 04.01.2022 Albert Salmi, Department of Electronics and Nanoengineering,
%                          Aalto University School of Electrical
%                          Engineering
% ------------------------------------------------------------------------
%

% Check file type and init ff_info
[path,name,ext] = fileparts(fname);
ff_info.filename = name;
ff_info.filepath = path;
ff_info.frequency = 0;
ff_info.element = 0;  

% CST farfield source file:
if ext == ".ffs"
    
    % Try to get element number from filename
    fname_parts = split(name, '_');
    fname_end_part = fname_parts{end};
    ff_info.element = sscanf(fname_end_part, '[%d]');
    
    % Read headerlines to ff_info struct
    fid = fopen(fname);
    line = '';
    hnum = 0;   % Calculate number of header lines in while-loop
    while strcmp(line, '// >> Phi, Theta, Re(E_Theta), Im(E_Theta), Re(E_Phi), Im(E_Phi): ') == 0
        line = fgetl(fid);
        switch line
            case '// Position'
                line = fgetl(fid);
                ff_info.position = sscanf(line, '%f');
            case '// Radiated/Accepted/Stimulated Power , Frequency '
                line = fgetl(fid);
                ff_info.Prad = sscanf(line, '%f');
                line = fgetl(fid);
                ff_info.Pacc = sscanf(line, '%f');
                line = fgetl(fid);
                ff_info.Pin = sscanf(line, '%f');
                line = fgetl(fid);
                ff_info.frequency = sscanf(line, '%f');
        end
        hnum = hnum+1;
    end
    hnum = hnum+3;
    fclose(fid);
    
    % Read the actual data
    [p,t,reEt,imEt,reEp,imEp] = readvars(fname, 'NumHeaderLines', hnum, 'FileType', 'text');

    theta = deg2rad(t);
    phi = deg2rad(p);

    Et = reEt + 1j*imEt;
    Ep = reEp + 1j*imEp;
    
    
    
% ASCII text file:
elseif ext == ".txt"
    
    % Get frequency and element number
    [fnamedata,n] = sscanf(name, 'farfield (f=%f) [%d]');
    if n==2
        ff_info.frequency = fnamedata(1);
        ff_info.element = fnamedata(2);     
    end

    % Read data
    [t,p,~,mag_Et,phase_Et,mag_Ep,phase_Ep,~] = readvars(fname, 'NumHeaderLines', 2);

    theta = deg2rad(t);
    phi = deg2rad(p);

    phase_Ep = deg2rad(phase_Ep);
    phase_Et = deg2rad(phase_Et);

    Et = mag_Et .* (cos(phase_Et) + 1j*sin(phase_Et));
    Ep = mag_Ep .* (cos(phase_Ep) + 1j*sin(phase_Ep));
    
else
    error('Unknown file extension');
end

end


function [theta, phi, Et, Ep, ff_info] = import_wavestudio_ff(fname)
%
% Read coordinates and E-fields from Starlab data (wavestudio format).
% ------------------------------------------------------------------------
% INPUT  fname  : string, full filename
%
% OUTPUT theta  : (:,1), theta values in radians (repeated)
%        phi    : (:,1), phi values in radians
%        Et     : (:,1,:), E-field theta component
%        Ep     : (:,1,:), E-field phi component
%        ff_info: struct, frequencies and other info
% ------------------------------------------------------------------------
% 12.09.2023 Albert Salmi, Department of Electronics and Nanoengineering,
%                          Aalto University School of Electrical
%                          Engineering
% ------------------------------------------------------------------------
%

% Check file type and init ff_info
[path,name,~] = fileparts(fname);
ff_info.filename = name;
ff_info.filepath = path;
ff_info.frequency = 0;
ff_info.element = 0;  
    
% Read data
[Az, El, freq, Et_re, Et_im, Ep_re, Ep_im] = readvars(fname, 'NumHeaderLines', 2);

% Wrap Azimuth from 0->pi to 0->2*pi
Az = Az - pi*(El < 0);
Az = mod(Az, 2*pi);
% Wrap elevation to 0->pi
El = abs(El);

data = [freq, Az, El, Et_re, Et_im, Ep_re, Ep_im];
new = data(data(:,3) == 0,:); % rows where elevation = 0
new(:,2) = new(:,2) + pi;
data = [data; new];

data = sortrows(data, [1,2,3]);

freq_u = unique(data(:,1));
phi = data(data(:,1)==data(1,1),2);
theta = data(data(:,1)==data(1,1),3);

Ep_all = data(:,4) .* (cos(data(:,6)) + 1j*sin(data(:,6)));
Et_all = data(:,5) .* (cos(data(:,7)) + 1j*sin(data(:,7)));

Ep = reshape(Ep_all, [], 1, length(freq_u));
Et = reshape(Et_all, [], 1, length(freq_u));
% 
ff_info.frequency = freq_u;

end