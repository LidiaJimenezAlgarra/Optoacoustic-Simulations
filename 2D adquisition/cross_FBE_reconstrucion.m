%% CROSS RECONSTRUCTION USING FBE METHOD OF Precision assessment of label-free psoriasis biomarkers with ultra-broadband optoacoustic mesoscopy
% run the code after the cross aqsuisition to avoid possible parameter errors.
% creates a kwave grid with the parameters of the grid of the cross
% it divides the dataset in two frequency bands and reconstructes them
% individually with the Kwave function kspacePlaneRecon.
% After it, the code computes the absolute value of the
% Hilbert transform to obtain the envelope and later create a MIP.
% Finally the projections are equalized and fused in each plane

%% LOAD DATA AND PARAMETERS
close all
%load("") %load the raw data
l=500; %l, t length 
dx=10e-6; %m
dy=10e-6; %m
vs=1500; %m/s
p0=1; %u.au
Rs=10e-6; %m diameter

pdetX=(-2e-3:dx:2e-3);  % number of grid points in the x (row) direction
pdetY=(-1e-3:dy:1e-3);  % number of grid points in the y (column) direction
Nx=length(pdetX);
Ny=length(pdetY);
t=linspace(0,1.0667e-06,l); 
Fs=1/(t(2)-t(1));

%% GRID
kgrid.Nt=size(S_f,2);
Ny=201;Nx=401; dx=10e-6;dy=dx;
kgrid = kWaveGrid(Nx, dx, Ny, dy);
kgrid.dt=1/Fs;
kgrid.Nt=500;

%% LOW FREQUENCY RECONSTRUCTION
S=reshape(S_f,80601,500);

Slow=timeFiltS(double(S), [10 30]*1e6, 1/Fs); %filter the frequency band
figure(1); imagesc(Slow'); colormap('gray');title('Low')
figure(2); imagesc(S'); colormap('gray');title('Normal')
sensor_data_rs = reshape(Slow, kgrid.Ny, kgrid.Nx, kgrid.Nt);
sensor_data_rs=sensor_data_rs(1:end,1:end,:);
Rlow = kspacePlaneRecon(sensor_data_rs, kgrid.dx,kgrid.dy, kgrid.dt, 1500,...
    'DataOrder', 'yzt');

%% HIGH FREQUENCY RECONSTRUCTION
Shigh=timeFiltS(S, [30 120]*1e6, 1/Fs); 
figure(3); imagesc(Shigh'); colormap('gray');title('High')

sensor_data_rs = reshape(Shigh, kgrid.Ny, kgrid.Nx, kgrid.Nt);
sensor_data_rs=sensor_data_rs(1:end,1:end,:);
Rhigh = kspacePlaneRecon(sensor_data_rs, kgrid.dy,kgrid.dx, kgrid.dt, 1500,...
    'DataOrder', 'yzt');

%% IMAGE PROCESSING
prueba=Rlow;
for i=1:201
    pruebalow(:,i,:)=abs(hilbert(Rlow(:,i,:)));
    pruebahigh(:,i,:)=abs(hilbert(Rhigh(:,i,:)));
end
for i=1:3 
  projlow=squeeze(max(pruebalow,[],i));
  projhigh=squeeze(max(pruebahigh,[],i));

  [alphaval] = alphacalc(projlow,projhigh);
  if i==1
  fus1= imfuse(projlow,alphaval*projhigh,'falsecolor','Scaling','joint','ColorChannels',[1 2 0]);
  fus1 = imadjust(fus1,[0.05 0.05 0; 0.15 0.15 1],[]); % R1 G1 B1, R2 G2 B2. R1=G1, G1=G2
  figure(1); 
  imagesc(pdetX*1e3,pdetY*1e3,fus1);title('XY')
  elseif i==2
  fus2= imfuse(projlow,alphaval*projhigh,'falsecolor','Scaling','joint','ColorChannels',[1 2 0]);
  fus2 = imadjust(fus2,[0.05 0.05 0; 0.15 0.15 1],[]); % R1 G1 B1, R2 G2 B2. R1=G1, G1=G2
  figure(2); 
  imagesc(pdetX*1e3,t*1e3*vs,fus2);title('XZ')
  else
  fus3= imfuse(projlow,alphaval*projhigh,'falsecolor','Scaling','joint','ColorChannels',[1 2 0]);
  fus3 = imadjust(fus3,[0.05 0.05 0; 0.15 0.15 1],[]); % R1 G1 B1, R2 G2 B2. R1=G1, G1=G2
  figure(3);
  imagesc(pdetY*1e3,t*1e3*vs,fus3);title('YZ');
  end
end

%% SAVING OF RAW DATA
% saveFolderData = ''; %SAVE FOLDER 
% fileName = datestr(now, 'yyyymmddHHMMSS');

% fileName2   = [ fileName '_SAVE_NAME.mat'];
% save([saveFolderData fileName2], 'VARIABLE_NAME');

% fileName2   = [ fileName '_SAVE_NAME.png'];
% saveas(figure(3), [saveFolderData fileName2]);
% 
% fileName2   = [ fileName '_SAVE_NAME.png'];
% saveas(figure(4), [saveFolderData fileName2]);
