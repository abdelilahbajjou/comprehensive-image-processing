function varargout = Matlab_GUI(varargin)
% MATLAB_GUI MATLAB code for Matlab_GUI.fig
%      MATLAB_GUI, by itself, creates a new MATLAB_GUI or raises the existing
%      singleton*.
%
%      H = MATLAB_GUI returns the handle to a new MATLAB_GUI or the handle to
%      the existing singleton*.
%
%      MATLAB_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MATLAB_GUI.M with the given input arguments.
%
%      MATLAB_GUI('Property','Value',...) creates a new MATLAB_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Matlab_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Matlab_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Matlab_GUI

% Last Modified by GUIDE v2.5 16-Feb-2025 18:20:07

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Matlab_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @Matlab_GUI_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% --- Executes just before Matlab_GUI is made visible.
function Matlab_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Matlab_GUI (see VARARGIN)

% Choose default command line output for Matlab_GUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% This sets up the initial plot - only do when we are invisible
% so window can get raised using Matlab_GUI.
if strcmp(get(hObject,'Visible'),'off')
    plot(rand(5));
end

% UIWAIT makes Matlab_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Matlab_GUI_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.entree);
cla;

popup_sel_index = get(handles.popupmenu1, 'Value');
switch popup_sel_index
    case 1
        plot(rand(5));
    case 2
        plot(sin(1:0.01:25.99));
    case 3
        bar(1:.5:10);
    case 4
        plot(membrane);
    case 5
        surf(peaks);
end


% --------------------------------------------------------------------
function Fichier_Callback(hObject, eventdata, handles)
% hObject    handle to Fichier (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Ouvrir_Callback(hObject, eventdata, handles)
% hObject    handle to Ouvrir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[file, path] = uigetfile('*.*', 'Select an Image');
% Read and display the image
imagePath = fullfile(path, file);
handles.i = imread(imagePath);

% Update handles structure
handles.image_entree = handles.i;
handles.output = hObject;
guidata(hObject, handles);

% Display the image in the 'entree' axes only
axes(handles.entree); % Set focus to 'entree' axes
imshow(handles.image_entree); % Display the image in 'entree'

% Save the updated handles structure
guidata(hObject, handles);



% --------------------------------------------------------------------
function Enregistrer_Callback(hObject, eventdata, handles)
% hObject    handle to Enregistrer (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
image = handles.image_sortie;
[file,path] = uiputfile('*.png','Enregistrer Votre Image ...');
imwrite(image, sprintf('%s',path,file),'png');

% --------------------------------------------------------------------
function Quitter_Callback(hObject, eventdata, handles)
% hObject    handle to Quitter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
delete(handles.figure1)


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1


% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
     set(hObject,'BackgroundColor','white');
end

set(hObject, 'String', {'plot(rand(5))', 'plot(sin(1:0.01:25))', 'bar(1:.5:10)', 'plot(membrane)', 'surf(peaks)'});


% --------------------------------------------------------------------
function Traitement_Callback(hObject, eventdata, handles)
% hObject    handle to Traitement (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function FiltragePassePas_Callback(hObject, eventdata, handles)
% hObject    handle to FiltragePassePas (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function FiltragePasseHaut_Callback(hObject, eventdata, handles)
% hObject    handle to FiltragePasseHaut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function FiltrageFrequentielle_Callback(hObject, eventdata, handles)
% hObject    handle to FiltrageFrequentielle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Bruit_Callback(hObject, eventdata, handles)
% hObject    handle to Bruit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function PoivreEtSel_Callback(hObject, eventdata, handles)
% hObject    handle to PoivreEtSel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    % Check if there is an input image
    if ~isfield(handles, 'image_entree') || isempty(handles.image_entree)
        errordlg('No image loaded.', 'Error');
        return;
    end

    % Get the input image
    img = handles.image_entree;

    % Convert the image to grayscale if it's RGB
    if size(img, 3) == 3
        img = rgb2gray(img);
    end

    % Get the size of the image
    [n, m] = size(img);

    % Define the noise density (percentage of pixels to corrupt)
    noise_density = 0.02; % 2% of the pixels will be noisy

    % Generate salt-and-pepper noise
    noisy_image = img; % Copy the original image
    num_noisy_pixels = round(noise_density * n * m); % Number of noisy pixels

    % Add salt (white pixels)
    salt_indices = randperm(n * m, num_noisy_pixels); % Random indices
    noisy_image(salt_indices) = 255; % Set to white

    % Add pepper (black pixels)
    pepper_indices = randperm(n * m, num_noisy_pixels); % Random indices
    noisy_image(pepper_indices) = 0; % Set to black

    % Display the noisy image in the "sortie" axes
    axes(handles.sortie);
    imshow(noisy_image);

    % Store the noisy image in handles
    handles.image_sortie = noisy_image;

    % Update the handles structure
    guidata(hObject, handles);


% --------------------------------------------------------------------
function Gaussien_Callback(hObject, eventdata, handles)
% hObject    handle to Gaussien (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
img = handles.image_entree;
image_traitee = imnoise(img, 'gaussian');
axes(handles.sortie)
imshow(image_traitee)
handles.output = hObject;
guidata(hObject,handles);

% --------------------------------------------------------------------
function FFT_Callback(hObject, eventdata, handles)
% hObject    handle to FFT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Assuming handles.image_entree contains the input image
img = handles.image_entree;

% Convert the image to grayscale if it is RGB
if size(img, 3) == 3
    img = rgb2gray(img);
end

% Compute the 2D FFT of the image
fft_img = fft2(double(img));

% Shift the zero-frequency component to the center of the spectrum
fft_shifted = fftshift(fft_img);

% Compute the magnitude spectrum (log scale for better visualization)
magnitude_spectrum = log(1 + abs(fft_shifted));

% Display the magnitude spectrum
axes(handles.sortie);
imshow(magnitude_spectrum, []); % Use [] to scale the display automatically
title('Magnitude Spectrum (FFT)');

% Store the FFT result in handles for later use
handles.fft_result = fft_shifted;
guidata(hObject, handles);

% --------------------------------------------------------------------
function FiltrePasseBasFrequentiel_Callback(hObject, eventdata, handles)
% hObject    handle to FiltrePasseBasFrequentiel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Assuming handles.image_entree contains the input image
img = handles.image_entree;

% Convert the image to grayscale if it is RGB
if size(img, 3) == 3
    img = rgb2gray(img);
end

% Compute the 2D FFT of the image
fft_img = fft2(double(img));

% Shift the zero-frequency component to the center of the spectrum
fft_shifted = fftshift(fft_img);

% Get the size of the image
[M, N] = size(img);

% Create a low-pass filter (e.g., a Gaussian or ideal low-pass filter)
% Example: Gaussian Low-Pass Filter
D0 = 30; % Cutoff frequency (adjust as needed)
[u, v] = meshgrid(-N/2:N/2-1, -M/2:M/2-1);
D = sqrt(u.^2 + v.^2); % Distance from the center
H = exp(-(D.^2) / (2 * (D0^2))); % Gaussian low-pass filter

% Apply the filter to the shifted FFT
filtered_fft = fft_shifted .* H;

% Shift the zero-frequency component back to the original position
filtered_fft_original = ifftshift(filtered_fft);

% Compute the inverse FFT to get the filtered image
filtered_img = ifft2(filtered_fft_original);

% Take the real part (to remove any small imaginary components due to numerical errors)
filtered_img = real(filtered_img);

% Display the filtered image
axes(handles.sortie);
imshow(uint8(filtered_img), []);
title('Low-Pass Filtered Image');

% Store the filtered image in handles for later use
handles.filtered_image = filtered_img;
guidata(hObject, handles);

% --------------------------------------------------------------------
function FiltrePasseHautFrequentiel_Callback(hObject, eventdata, handles)
% hObject    handle to FiltrePasseHautFrequentiel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Assuming handles.image_entree contains the input image
img = handles.image_entree;

% Convert the image to grayscale if it is RGB
if size(img, 3) == 3
    img = rgb2gray(img);
end

% Compute the 2D FFT of the image
fft_img = fft2(double(img));

% Shift the zero-frequency component to the center of the spectrum
fft_shifted = fftshift(fft_img);

% Get the size of the image
[M, N] = size(img);

% Create a high-pass filter (e.g., a Gaussian or ideal high-pass filter)
% Example: Gaussian High-Pass Filter
D0 = 30; % Cutoff frequency (adjust as needed)
[u, v] = meshgrid(-N/2:N/2-1, -M/2:M/2-1);
D = sqrt(u.^2 + v.^2); % Distance from the center
H = 1 - exp(-(D.^2) / (2 * (D0^2))); % Gaussian high-pass filter

% Apply the filter to the shifted FFT
filtered_fft = fft_shifted .* H;

% Shift the zero-frequency component back to the original position
filtered_fft_original = ifftshift(filtered_fft);

% Compute the inverse FFT to get the filtered image
filtered_img = ifft2(filtered_fft_original);

% Take the real part (to remove any small imaginary components due to numerical errors)
filtered_img = real(filtered_img);

% Display the filtered image
axes(handles.sortie);
imshow(uint8(filtered_img), []);
title('High-Pass Filtered Image');

% Store the filtered image in handles for later use
handles.filtered_image = filtered_img;
guidata(hObject, handles);

% --------------------------------------------------------------------
function Laplacien_Callback(hObject, eventdata, handles)
% hObject    handle to Laplacien (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    % Check if there is an input image
    if ~isfield(handles, 'image_entree') || isempty(handles.image_entree)
        errordlg('No image loaded.', 'Error');
        return;
    end

    % Get the input image
    img = handles.image_entree;

    % Convert the image to grayscale if it's RGB
    if size(img, 3) == 3
        img = rgb2gray(img);
    end

    % Convert the image to double for processing
    image_double = double(img);

    % Get the size of the image
    [n, m] = size(image_double);

    % Initialize the output matrix
    S = zeros(n, m);

    % Define the Laplacian kernel
    Laplacien = [-1, -1, -1; -1, 8, -1; -1, -1, -1];

    % Apply the Laplacian filter
    for i = 2:n-1
        for j = 2:m-1
            % Extract the 3x3 neighborhood
            Vecteur = image_double((i-1:i+1), (j-1:j+1));
            
            % Apply the Laplacian kernel
            produit = Vecteur .* Laplacien;
            
            % Sum the result and store it
            S(i, j) = sum(produit(:));
        end
    end

    % Thresholding
    S(abs(S) < 100) = 0; % Set values below 100 to 0

    % Convert the result back to uint8
    image_traitee = uint8(S);

    % Clear previous image from the "sortie" axes
    cla(handles.sortie); % Clear current axes

    % Display the new processed image in "sortie" axes
    axes(handles.sortie);
    imshow(image_traitee);

    % Store the processed image in handles
    handles.image_sortie = image_traitee;

    % Update the handles structure
    guidata(hObject, handles);



% --------------------------------------------------------------------
function Gradient_Callback(hObject, eventdata, handles)
% hObject    handle to Gradient (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    % Check if there is an input image
    if ~isfield(handles, 'image_entree') || isempty(handles.image_entree)
        errordlg('No image loaded.', 'Error');
        return;
    end

    % Get the input image
    img = handles.image_entree;

    % Convert the image to grayscale if it's RGB
    if size(img, 3) == 3
        img = rgb2gray(img);
    end

    % Convert the image to double for processing
    image_double = double(img);

    % Get the size of the image
    [n, m] = size(image_double);

    % Initialize the output matrices
    S_horizontal = zeros(n, m); % Horizontal gradient
    S_vertical = zeros(n, m);   % Vertical gradient
    S = zeros(n, m);            % Gradient magnitude

    % Define the horizontal and vertical Sobel-like filters
    composition_hor = [0, 0, 0; -1, 0, 1; 0, 0, 0]; % Horizontal filter
    composition_ver = [0, -1, 0; 0, 0, 0; 0, 1, 0]; % Vertical filter

    % Apply the filters to compute horizontal and vertical gradients
    for i = 2:n-1
        for j = 2:m-1
            % Extract the 3x3 neighborhood
            neighborhood = image_double(i-1:i+1, j-1:j+1);
            
            % Compute the horizontal gradient
            S_horizontal(i, j) = sum(sum(neighborhood .* composition_hor)); % Use sum(sum(...)) for MATLAB 2013
            
            % Compute the vertical gradient
            S_vertical(i, j) = sum(sum(neighborhood .* composition_ver)); % Use sum(sum(...)) for MATLAB 2013
        end
    end

    % Compute the gradient magnitude
    S = sqrt(S_horizontal.^2 + S_vertical.^2);

    % Thresholding: Set values below 85 to 0
    S(S < 85) = 0;

    % Normalize the gradient magnitude to the range [0, 255]
    S = S / max(S(:)) * 255;

    % Convert the result to uint8
    image_traitee = uint8(S);

    % Clear previous image from the "sortie" axes
    cla(handles.sortie); % Clear current axes

    % Display the new processed image in "sortie" axes
    axes(handles.sortie);
    imshow(image_traitee);

    % Store the processed image in handles
    handles.image_sortie = image_traitee;

    % Update the handles structure
    guidata(hObject, handles);
% --------------------------------------------------------------------
function Prewitt_Callback(hObject, eventdata, handles)
% hObject    handle to Prewitt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    % Check if there is an input image
    if ~isfield(handles, 'image_entree') || isempty(handles.image_entree)
        errordlg('No image loaded.', 'Error');
        return;
    end

    % Get the input image
    img = handles.image_entree;

    % Convert the image to grayscale if it's RGB
    if size(img, 3) == 3
        img = rgb2gray(img);
    end

    % Convert the image to double for processing
    image_double = double(img);

    % Get the size of the image
    [n, m] = size(image_double);

    % Initialize the output matrices
    S_horizontal = zeros(n, m); % Horizontal gradient
    S_vertical = zeros(n, m);   % Vertical gradient
    S = zeros(n, m);            % Gradient magnitude

    % Define the Prewitt kernels
    composition_hor = [-1, 0, 1; -1, 0, 1; -1, 0, 1]; % Horizontal kernel
    composition_ver = [-1, -1, -1; 0, 0, 0; 1, 1, 1]; % Vertical kernel

    % Apply the Prewitt operator
    for i = 2:(n-1)
        for j = 2:(m-1)
            % Extract the 3x3 neighborhood
            neighborhood = image_double(i-1:i+1, j-1:j+1);
            
            % Compute the horizontal gradient
            S_horizontal(i, j) = sum(sum(neighborhood .* composition_hor));
            
            % Compute the vertical gradient
            S_vertical(i, j) = sum(sum(neighborhood .* composition_ver));
        end
    end

    % Compute the gradient magnitude
    S = sqrt(S_horizontal.^2 + S_vertical.^2);

    % Normalize the gradient magnitude to the range [0, 255]
    S = S / max(S(:)) * 255;

    % Thresholding: Set values below 100 to 0
    S(S < 100) = 0;

    % Convert the result to uint8
    image_traitee = uint8(S);

    % Clear previous image from the "sortie" axes
    cla(handles.sortie); % Clear current axes

    % Display the new processed image in "sortie" axes
    axes(handles.sortie);
    imshow(image_traitee);

    % Store the processed image in handles
    handles.image_sortie = image_traitee;

    % Update the handles structure
    guidata(hObject, handles);

% --------------------------------------------------------------------
function Roberts_Callback(hObject, eventdata, handles)
% hObject    handle to Roberts (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    % Check if there is an input image
    if ~isfield(handles, 'image_entree') || isempty(handles.image_entree)
        errordlg('No image loaded.', 'Error');
        return;
    end

    % Get the input image
    img = handles.image_entree;

    % Convert the image to grayscale if it's RGB
    if size(img, 3) == 3
        img = rgb2gray(img);
    end

    % Convert the image to double for processing
    image_double = double(img);

    % Get the size of the image
    [n, m] = size(image_double);

    % Initialize the output matrices
    S_horizontal = zeros(n, m); % Horizontal gradient
    S_vertical = zeros(n, m);   % Vertical gradient
    S = zeros(n, m);            % Gradient magnitude

    % Define the Roberts Cross kernels
    composition_hor = [-1, 0; 0, 1]; % Horizontal kernel
    composition_ver = [0, -1; 1, 0]; % Vertical kernel

    % Apply the Roberts Cross operator
    for i = 1:(n-1)
        for j = 1:(m-1)
            % Extract the 2x2 neighborhood
            neighborhood = image_double(i:i+1, j:j+1);
            
            % Compute the horizontal gradient
            S_horizontal(i, j) = sum(sum(neighborhood .* composition_hor));
            
            % Compute the vertical gradient
            S_vertical(i, j) = sum(sum(neighborhood .* composition_ver));
        end
    end

    % Compute the gradient magnitude
    S = sqrt(S_horizontal.^2 + S_vertical.^2);

    % Normalize the gradient magnitude to the range [0, 255]
    S = S / max(S(:)) * 255;

    % Thresholding: Set values below 100 to 0
    S(S < 100) = 0;

    % Convert the result to uint8
    image_traitee = uint8(S);

    % Clear previous image from the "sortie" axes
    cla(handles.sortie); % Clear current axes

    % Display the new processed image in "sortie" axes
    axes(handles.sortie);
    imshow(image_traitee);

    % Store the processed image in handles
    handles.image_sortie = image_traitee;

    % Update the handles structure
    guidata(hObject, handles);

% --------------------------------------------------------------------
function Lineaire_Callback(hObject, eventdata, handles)
% hObject    handle to Lineaire (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function NonLineaire_Callback(hObject, eventdata, handles)
% hObject    handle to NonLineaire (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function NiveauDeGris_Callback(hObject, eventdata, handles)
% hObject    handle to NiveauDeGris (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
image=handles.image_entree;
l = length(size(image));
if l==3
    image_traitee=rgb2gray(image);
elseif l==2
   image_traitee=image;
end
axes(handles.sortie);
imshow(image_traitee);
handles.output = hObject;
guidata(hObject, handles);

% --------------------------------------------------------------------
function Contraste_Callback(hObject, eventdata, handles)
% hObject    handle to Contraste (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
img=handles.image_entree;
 l = length(size(img));
        if l==3
            image = rgb2gray(img);
        elseif l==2
            image = img;
        end
axes(handles.entree);
imshow(image);
image = double(image);
[l,c]=size(image);
image = double(image);
%mi=min(image,[],'all');
%ma=max(image,[],'all');
mi = min(image(:));  % Minimum value of the image
ma = max(image(:)); 
T=255/(ma-mi);
p=image;
for i=1:l
    for j=1:c
       p(i,j) = T*(image(i,j)-mi); 
    end
end
handles.image_sortie=uint8(p); 
axes(handles.sortie);
imshow(handles.image_sortie);
handles.output = hObject;
guidata(hObject, handles);

% --------------------------------------------------------------------
function Luminosite_Callback(hObject, eventdata, handles)
% hObject    handle to Luminosite (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
image=handles.image_entree;
[n,m]=size(image);
image = double(image);
p=image;
for i=1:n
    for j=1:m
        pixel=image(i,j)+50;
         if(pixel>255)
            pixel=255;
         else if (pixel<0)
                pixel=0;
              end 
          end
       p(i,j)=pixel;    
    end
end  
image_traitee=uint8(p); 
axes(handles.sortie);
imshow(image_traitee);
handles.image_sortie=image_traitee;
handles.output = hObject;
guidata(hObject, handles);

% --------------------------------------------------------------------
function Negatif_Callback(hObject, eventdata, handles)
% hObject    handle to Negatif (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
image=handles.image_entree;
[n,m]=size(image);
image_double = double(image);
p=image_double;
for i=1:n
   for j=1:m
     p(i,j)=-image_double(i,j)+255;
    end
 end 
image_traitee=uint8(p); 
axes(handles.sortie);
imshow(image_traitee);
handles.output = hObject;
guidata(hObject, handles);

% --------------------------------------------------------------------
function Histogramme_Callback(hObject, eventdata, handles)
% hObject    handle to Histogramme (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
img = handles.image_entree;
        l = length(size(img));
        if l==3
            U = rgb2gray(img);
        elseif l==2
            U = img;
        end
axes(handles.entree);
imshow(U);

[n,m]=size(U);
t=double(U);
vecteur=1:256;
p=0;
for k=0:255 
    for i=1:n
        for j=1:m
            if t(i,j)==k 
               p=p+1;
            end
        end
    end
    vecteur(k+1)=p;
    p=0;
end
axes(handles.sortie);
plot(vecteur);
handles.output = hObject;
guidata(hObject, handles);
% --------------------------------------------------------------------
function Binarisation_Callback(hObject, eventdata, handles)
% hObject    handle to Binarisation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
image = handles.image_entree;

% Convert image to double before calculations
image = double(image);

% Flatten the image to a 1D vector
image_flattened = image(:);

% Calculating m0, m1, m2, m3 using the flattened image
m0 = 1;
m1 = mean(image_flattened);         % Mean of the image
m2 = mean(image_flattened.^2);      % Mean of the squared image
m3 = mean(image_flattened.^3); 

% Calculating C1 and C0
C1 = (m3 - (m1 * m2)) / (m2 - m1);
C0 = (-m2 - (C1 * m1)) / m0;

% Calculating z1 and z2
z1 = (-C1 - sqrt(C1^2 - 4 * C0)) / 2;
z2 = (-C1 + sqrt(C1^2 - 4 * C0)) / 2;

% Calculating the threshold
%seuil=(z1 + z2) / 2
seuil = 120;

% Binarize the image based on the threshold
bin = (image > seuil) * 255;

% Display the binarized image
axes(handles.sortie);
handles.image_sortie = bin;
imshow(handles.image_sortie);

% Update the handles structure
handles.output = hObject;
guidata(hObject, handles);



% --------------------------------------------------------------------
function Median_Callback(hObject, eventdata, handles)
% hObject    handle to Median (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Median 
img = handles.image_entree;
[n,m] = size(img);
image_double = double(img);
S = image_double;
for i = 2 : n-1
    for j = 2 : m-1
      fenetre = image_double(i-1:i+1,j-1:j+1);
      U = [fenetre(1,:) fenetre(2,:) fenetre(3,:)];
      sort(U);
      resultat = median(U);
      S(i,j) = resultat ;
    end 
end
image_traitee = uint8(S);
axes(handles.sortie);
imshow(image_traitee);
handles.image_sortie = image_traitee;
handles.output = hObject;
guidata(hObject, handles);

% --------------------------------------------------------------------
function Moyenneur3_3_Callback(hObject, eventdata, handles)
% hObject    handle to Moyenneur3_3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
img = handles.image_entree;
[n,m] = size(img);
image_double = double(img);
S = image_double;
M = (1/9)*[1, 1, 1 ; 1, 1, 1 ; 1, 1, 1 ];
for i = 2 : n-1
    for j = 2 : m-1
      fenetre = image_double(i-1:i+1,j-1:j+1);
      produit = fenetre.*M;
      S(i,j) = sum(produit(:));
    end 
end
image_traitee = uint8(S);
axes(handles.sortie);
imshow(image_traitee);
handles.image_sortie = image_traitee;
handles.output = hObject;
guidata(hObject, handles);

% --------------------------------------------------------------------
function Moyenneur5_5_Callback(hObject, eventdata, handles)
% hObject    handle to Moyenneur5_5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
img = handles.image_entree;
[n,m] = size(img);
image_double = double(img);
S = image_double;
M = (1/25)*[1, 1, 1, 1, 1 ; 1, 1, 1, 1, 1 ; 1, 1, 1, 1, 1 ; 1, 1, 1, 1, 1 ; 1, 1, 1, 1, 1 ];
for i = 3 : n-2
    for j = 3 : m-2
      fenetre = image_double(i-2:i+2,j-2:j+2);
      produit = fenetre.*M;
      S(i,j) = sum(produit(:));
    end 
end
handles.image_sortie = uint8(S);
axes(handles.sortie);
imshow(handles.image_sortie);
handles.output = hObject;
guidata(hObject, handles);

% --------------------------------------------------------------------
function Gaussien3_3_Callback(hObject, eventdata, handles)
% hObject    handle to Gaussien3_3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
img = handles.image_entree;
[n,m] = size(img);
image_double = double(img);
S = image_double;
G = (1/16)*[1, 2, 1 ;2, 4, 2 ; 1, 2, 1];
for i = 2 : n-1
    for j = 2 : m-1
      fenetre = image_double(i-1:i+1,j-1:j+1);
      produit = fenetre.*G;
      S(i,j) = sum(produit(:));
    end 
end
image_traitee = uint8(S);
axes(handles.sortie);
imshow(image_traitee);
handles.image_sortie = image_traitee;
handles.output = hObject;
guidata(hObject, handles);

% --------------------------------------------------------------------
function Gaussien5_5_Callback(hObject, eventdata, handles)
% hObject    handle to Gaussien5_5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
img = handles.image_entree;
[n,m] = size(img);
image_double = double(img);
S = image_double;
G = (1/256)*[1, 4, 6, 4, 1 ; 4, 16, 24, 16, 4 ; 6, 24, 36, 24, 6 ; 4, 16, 24, 16, 4 ; 1, 4, 6, 4, 1];
for i = 3 : n-2
    for j = 3 : m-2
      fenetre = image_double(i-2:i+2,j-2:j+2);
      produit = fenetre.*G;
      S(i,j) = sum(produit(:));
    end 
end
image_traitee = uint8(S);
axes(handles.sortie);
imshow(image_traitee);
handles.image_sortie = image_traitee ;
handles.output = hObject;
guidata(hObject, handles);

% --------------------------------------------------------------------
function Conique_Callback(hObject, eventdata, handles)
% hObject    handle to Conique (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
img = handles.image_entree;
[n,m] = size(img);
image_double = double(img);
S = image_double;
C = (1/25)*[0, 0, 1, 0, 0 ; 0, 2, 2, 2, 0 ; 1, 2, 5, 2, 1 ; 0, 2, 2, 2, 0 ; 0, 0, 1, 0, 0];
for i = 3 : n-2
    for j = 3 : m-2
      fenetre = image_double(i-2:i+2,j-2:j+2);
      produit = fenetre.*C;
      S(i,j) = sum(produit(:));
    end 
end
image_traitee = uint8(S);
axes(handles.sortie);
imshow(image_traitee);
handles.image_sortie = image_traitee ;
handles.output = hObject;
guidata(hObject, handles);

% --------------------------------------------------------------------
function Pyramidal_Callback(hObject, eventdata, handles)
% hObject    handle to Pyramidal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
img = handles.image_entree;
[n,m] = size(img);
image_double = double(img);
S = image_double;
P = (1/81)*[ 1,2,3,2,1 ; 2,4, 6, 4, 2 ; 3, 6, 9, 6, 3 ; 2, 4, 6, 4, 2 ; 1, 2, 3, 2, 1 ] ;
for i = 3 : n-2
    for j = 3 : m-2
      fenetre = image_double(i-2:i+2,j-2:j+2);
      produit = fenetre.*P;
      S(i,j) = sum(produit(:));
    end 
end
image_traitee = uint8(S);
axes(handles.sortie);
imshow(image_traitee);
handles.image_sortie = image_traitee ;
handles.output = hObject;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function entree_CreateFcn(hObject, eventdata, handles)
% hObject    handle to entree (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
set(handles.entree, 'XTick', [], 'YTick', []);


% --- Executes during object creation, after setting all properties.
function sortie_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sortie (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --------------------------------------------------------------------
function Sobel_Callback(hObject, eventdata, handles)
% hObject    handle to Sobel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    % Check if there is an input image
    if ~isfield(handles, 'image_entree') || isempty(handles.image_entree)
        errordlg('No image loaded.', 'Error');
        return;
    end

    % Get the input image
    img = handles.image_entree;

    % Convert the image to grayscale if it's RGB
    if size(img, 3) == 3
        img = rgb2gray(img);
    end

    % Convert the image to double for processing
    image_double = double(img);

    % Get the size of the image
    [n, m] = size(image_double);

    % Initialize the output matrices
    S_horizontal = zeros(n, m); % Horizontal gradient
    S_vertical = zeros(n, m);   % Vertical gradient
    S = zeros(n, m);            % Gradient magnitude

    % Define the Sobel kernels
    composition_hor = [-1, 0, 1; -2, 0, 2; -1, 0, 1]; % Horizontal kernel
    composition_ver = [-1, -2, -1; 0, 0, 0; 1, 2, 1]; % Vertical kernel

    % Apply the Sobel operator
    for i = 2:(n-1)
        for j = 2:(m-1)
            % Extract the 3x3 neighborhood
            neighborhood = image_double(i-1:i+1, j-1:j+1);
            
            % Compute the horizontal gradient
            S_horizontal(i, j) = sum(sum(neighborhood .* composition_hor));
            
            % Compute the vertical gradient
            S_vertical(i, j) = sum(sum(neighborhood .* composition_ver));
        end
    end

    % Compute the gradient magnitude
    S = sqrt(S_horizontal.^2 + S_vertical.^2);

    % Normalize the gradient magnitude to the range [0, 255]
    S = S / max(S(:)) * 255;

    % Thresholding: Set values below 100 to 0
    S(S < 100) = 0;

    % Convert the result to uint8
    image_traitee = uint8(S);

    % Clear previous image from the "sortie" axes
    cla(handles.sortie); % Clear current axes

    % Display the new processed image in "sortie" axes
    axes(handles.sortie);
    imshow(image_traitee);

    % Store the processed image in handles
    handles.image_sortie = image_traitee;

    % Update the handles structure
    guidata(hObject, handles);


% --------------------------------------------------------------------
function Compresser_Callback(hObject, eventdata, handles)
% hObject    handle to Compresser (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Save the image in PNG format (lossless compresser)
    % Get the current image from handles
    img = handles.image_entree;

    % Ask the user to choose the compresser type
    choice = questdlg('Choose compression type:', ...
                      'Compression', ...
                      'Lossless (PNG)', 'Lossy (JPEG)', 'Quantization', 'Lossless (PNG)');

    % Perform compresser based on user choice
    switch choice
        case 'Lossless (PNG)'
            % Save the image in PNG format (lossless compresser)
            compressed_image_path = 'compressed_image.png';
            imwrite(img, compressed_image_path, 'png');
            compressed_image = imread(compressed_image_path);  % Read back the compressed image
            disp('Image saved with lossless compression (PNG).');

        case 'Lossy (JPEG)'
            % Save the image in JPEG format (lossy compresser)
            compressed_image_path = 'compressed_image.jpg';
            imwrite(img, compressed_image_path, 'jpg', 'Quality', 90);  % Adjust quality (0-100)
            compressed_image = imread(compressed_image_path);  % Read back the compressed image
            disp('Image saved with lossy compression (JPEG).');

        case 'Quantization'
            % Custom compresser by reducing intensity levels
            levels = 32;  % Reduce to 32 intensity levels
            quantized_image = round(double(img) / (256 / levels)) * (256 / levels);
            compressed_image = uint8(quantized_image);  % Convert back to uint8
            disp('Image compressed by quantization.');
    end

    % Display the compressed image in the "Image Traitée" axes
    axes(handles.sortie);  % Use the axes handle for "Image Traitée"
    imshow(compressed_image);
    title('Compressed Image');

    % Save the compressed image in handles for further processing
    handles.image_sortie = compressed_image;

    % Update handles structure
    guidata(hObject, handles);


% --------------------------------------------------------------------
function IsolerLesObjets_Callback(hObject, eventdata, handles)
% hObject    handle to IsolerLesObjets (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    % Get the current image
    img = handles.image_entree;

    % Convert to grayscale if necessary
    if size(img, 3) == 3
        img = rgb2gray(img);
    end

    % Ensure the image is in double format for processing
    img = double(img);

    % Compute the gradient magnitude (to highlight edges)
    hy = fspecial('sobel');  % Sobel filter for vertical edges
    hx = hy';                % Sobel filter for horizontal edges
    gy = imfilter(img, hy, 'replicate');
    gx = imfilter(img, hx, 'replicate');
    gradient_magnitude = sqrt(gx.^2 + gy.^2);

    % Apply the watershed transformation on the gradient magnitude
    watershed_labels = watershed(gradient_magnitude);

    % Colorize the regions for visualization
    segmented_image = label2rgb(watershed_labels, 'jet', 'w', 'shuffle');

    % Display the segmented image in the "Image Traitée" axes
    axes(handles.sortie);
    imshow(segmented_image);

    % Save the segmented image in handles for further processing
    handles.image_sortie = segmented_image;

    % Update handles structure
    guidata(hObject, handles);


    

% --------------------------------------------------------------------
function TriangleDePascal_Callback(hObject, eventdata, handles)
% hObject    handle to TriangleDePascal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --- Callback function for Triangle de Pascal
    % hObject    handle to TriangleDePascal (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % --- Callback function for Triangle de Pascal
    % Get the current image
    img = handles.image_entree;

    % Convert to grayscale if necessary
    if size(img, 3) == 3
        img = rgb2gray(img);
    end

    % Step 1: Generate Pascal's Triangle (5x5 filter)
    n = 5;  % Size of the filter (5x5)
    pascal_triangle = zeros(n, n);
    for i = 1:n
        for j = 1:i
            if j == 1 || j == i
                pascal_triangle(i, j) = 1;
            else
                pascal_triangle(i, j) = pascal_triangle(i-1, j-1) + pascal_triangle(i-1, j);
            end
        end
    end

    % Step 2: Create the smoothing filter
    pascal_row = pascal_triangle(n, :);  % Get the nth row
    pascal_row = pascal_row / sum(pascal_row);  % Normalize the row
    smoothing_filter = pascal_row' * pascal_row;  % Create 2D filter

    % Step 3: Apply the smoothing filter to the image
    smoothed_image = imfilter(double(img), smoothing_filter, 'conv', 'replicate');

    % Convert back to uint8 for display
    smoothed_image = uint8(smoothed_image);

    % Display the smoothed image in the "Image Traitée" axes
    axes(handles.sortie);
    imshow(smoothed_image);

    % Save the smoothed image in handles for further processing
    handles.image_sortie = smoothed_image;

    % Update handles structure
    guidata(hObject, handles);


% --------------------------------------------------------------------
function Kirsch_Callback(hObject, eventdata, handles)
% hObject    handle to Kirsch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    % hObject    handle to Kirsch (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % --- Callback function for Kirsch edge detection
    % Get the current image
    img = handles.image_entree;

    % Convert to grayscale if necessary
    if size(img, 3) == 3
        img = rgb2gray(img);
    end

    % Convert the image to double for processing
    img = double(img);

    % Define the 8 Kirsch kernels (3x3)
    kirsch_kernels = {
        [ 5,  5,  5; -3,  0, -3; -3, -3, -3],  % North
        [ 5,  5, -3;  5,  0, -3; -3, -3, -3],  % Northeast
        [ 5, -3, -3;  5,  0, -3;  5, -3, -3],  % East
        [-3, -3, -3;  5,  0, -3;  5,  5, -3],  % Southeast
        [-3, -3, -3; -3,  0, -3;  5,  5,  5],  % South
        [-3, -3, -3; -3,  0,  5; -3,  5,  5],  % Southwest
        [-3, -3,  5; -3,  0,  5; -3, -3,  5],  % West
        [-3,  5,  5; -3,  0,  5; -3, -3, -3]   % Northwest
    };

    % Initialize the edge magnitude image
    edge_magnitude = zeros(size(img));

    % Apply each Kirsch kernel and compute the maximum response
    for i = 1:8
        kernel = kirsch_kernels{i};
        filtered_image = imfilter(img, kernel, 'conv', 'replicate');
        edge_magnitude = max(edge_magnitude, filtered_image);
    end

    % Normalize the edge magnitude to the range [0, 255]
    edge_magnitude = uint8(255 * mat2gray(edge_magnitude));

    % Display the edge magnitude image in the "Image Traitée" axes
    axes(handles.sortie);
    imshow(edge_magnitude);
    % Save the edge magnitude image in handles for further processing
    handles.image_sortie = edge_magnitude;

    % Update handles structure
    guidata(hObject, handles);


% --------------------------------------------------------------------
function Marr_Hildreth_Callback(hObject, eventdata, handles)
% hObject    handle to Marr_Hildreth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    % Get the current image
    img = handles.image_entree;

    % Convert to grayscale if necessary
    if size(img, 3) == 3
        img = rgb2gray(img);
    end

    % Convert the image to double for processing
    img = double(img);

    % Step 1: Apply Gaussian smoothing
    sigma = 2;  % Standard deviation of the Gaussian filter
    gaussian_filter = fspecial('gaussian', [5 5], sigma);  % Create a 5x5 Gaussian filter
    smoothed_image = imfilter(img, gaussian_filter, 'conv', 'replicate');

    % Step 2: Apply the Laplacian operator
    laplacian_filter = fspecial('laplacian', 0);  % Create a Laplacian filter
    laplacian_image = imfilter(smoothed_image, laplacian_filter, 'conv', 'replicate');

    % Step 3: Detect zero-crossings
    zero_crossings = edge(smoothed_image, 'zerocross', 0, laplacian_filter);

    % Convert the zero-crossings image to uint8 for display
    zero_crossings = uint8(255 * zero_crossings);

    % Display the zero-crossings image in the "Image Traitée" axes
    axes(handles.sortie);
    imshow(zero_crossings);
    title('Edge Detection (Marr-Hildreth)');

    % Save the zero-crossings image in handles for further processing
    handles.image_sortie = zero_crossings;

    % Update handles structure
    guidata(hObject, handles);


% --------------------------------------------------------------------
function Canny_Callback(hObject, eventdata, handles)
% hObject    handle to Canny (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    % Get the current image
    img = handles.image_entree;

    % Convert to grayscale if necessary
    if size(img, 3) == 3
        img = rgb2gray(img);
    end

    % Step 1: Apply Gaussian smoothing (alternative to imgaussfilt)
    sigma = 1.5;  % Standard deviation of the Gaussian filter
    filter_size = 2 * ceil(3 * sigma) + 1;  % Determine filter size based on sigma
    gaussian_filter = fspecial('gaussian', [filter_size filter_size], sigma);  % Create Gaussian filter
    smoothed_image = imfilter(double(img), gaussian_filter, 'conv', 'replicate');  % Apply filter

    % Step 2: Compute gradient magnitude and direction using Sobel filters
    [Gx, Gy] = imgradientxy(smoothed_image, 'sobel');  % Gradient in x and y directions
    [gradient_magnitude, gradient_direction] = imgradient(Gx, Gy);  % Gradient magnitude and direction

    % Step 3: Use MATLAB's built-in Canny edge detector
    edge_image = edge(smoothed_image, 'canny');  % Apply Canny edge detection

    % Convert the edge image to uint8 for display
    edge_image = uint8(255 * edge_image);

    % Display the edge-detected image in the "Image Traitée" axes
    axes(handles.sortie);
    imshow(edge_image);
    title('Edge Detection (Canny)');

    % Save the edge-detected image in handles for further processing
    handles.image_sortie = edge_image;

    % Update handles structure
    guidata(hObject, handles);


% --------------------------------------------------------------------
function FiltrePasseBasIdeal_Callback(hObject, eventdata, handles)
% hObject    handle to FiltrePasseBasIdeal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --------------------------------------------------------------------
img = handles.image_entree;

% Convert the image to grayscale if it is RGB
if size(img, 3) == 3
    img = rgb2gray(img);
end

% Compute the 2D FFT of the image
fft_img = fft2(double(img));

% Shift the zero-frequency component to the center of the spectrum
fft_shifted = fftshift(fft_img);

% Get the size of the image
[M, N] = size(img);

% Create an ideal low-pass filter
D0 = 30; % Cutoff frequency (adjust as needed)
[u, v] = meshgrid(-N/2:N/2-1, -M/2:M/2-1);
D = sqrt(u.^2 + v.^2); % Distance from the center
H = double(D <= D0); % Ideal low-pass filter (1 inside the cutoff, 0 outside)

% Apply the filter to the shifted FFT
filtered_fft = fft_shifted .* H;

% Shift the zero-frequency component back to the original position
filtered_fft_original = ifftshift(filtered_fft);

% Compute the inverse FFT to get the filtered image
filtered_img = ifft2(filtered_fft_original);

% Take the real part (to remove any small imaginary components due to numerical errors)
filtered_img = real(filtered_img);

% Display the filtered image
axes(handles.sortie);
imshow(uint8(filtered_img), []);

% Store the filtered image in handles for later use
handles.filtered_image = filtered_img;
guidata(hObject, handles);


% --------------------------------------------------------------------
function FiltrePasseBasDeButterwoth_Callback(hObject, eventdata, handles)
% hObject    handle to FiltrePasseBasDeButterwoth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --------------------------------------------------------------------
img = handles.image_entree;

% Convert the image to grayscale if it is RGB
if size(img, 3) == 3
    img = rgb2gray(img);
end

% Compute the 2D FFT of the image
fft_img = fft2(double(img));

% Shift the zero-frequency component to the center of the spectrum
fft_shifted = fftshift(fft_img);

% Get the size of the image
[M, N] = size(img);

% Create a Butterworth low-pass filter
D0 = 30; % Cutoff frequency (adjust as needed)
n = 2; % Order of the Butterworth filter (adjust as needed)
[u, v] = meshgrid(-N/2:N/2-1, -M/2:M/2-1);
D = sqrt(u.^2 + v.^2); % Distance from the center
H = 1 ./ (1 + (D ./ D0).^(2*n)); % Butterworth low-pass filter

% Apply the filter to the shifted FFT
filtered_fft = fft_shifted .* H;

% Shift the zero-frequency component back to the original position
filtered_fft_original = ifftshift(filtered_fft);

% Compute the inverse FFT to get the filtered image
filtered_img = ifft2(filtered_fft_original);

% Take the real part (to remove any small imaginary components due to numerical errors)
filtered_img = real(filtered_img);

% Display the filtered image
axes(handles.sortie);
imshow(uint8(filtered_img), []);

% Store the filtered image in handles for later use
handles.filtered_image = filtered_img;
guidata(hObject, handles);


% --------------------------------------------------------------------
function FiltrePasseHautIdeal_Callback(hObject, eventdata, handles)
% hObject    handle to FiltrePasseHautIdeal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Assuming handles.image_entree contains the input image
img = handles.image_entree;

% Convert the image to grayscale if it is RGB
if size(img, 3) == 3
    img = rgb2gray(img);
end

% Compute the 2D FFT of the image
fft_img = fft2(double(img));

% Shift the zero-frequency component to the center of the spectrum
fft_shifted = fftshift(fft_img);

% Get the size of the image
[M, N] = size(img);

% Create an ideal high-pass filter
D0 = 30; % Cutoff frequency (adjust as needed)
[u, v] = meshgrid(-N/2:N/2-1, -M/2:M/2-1);
D = sqrt(u.^2 + v.^2); % Distance from the center
H = double(D > D0); % Ideal high-pass filter (1 outside the cutoff, 0 inside)

% Apply the filter to the shifted FFT
filtered_fft = fft_shifted .* H;

% Shift the zero-frequency component back to the original position
filtered_fft_original = ifftshift(filtered_fft);

% Compute the inverse FFT to get the filtered image
filtered_img = ifft2(filtered_fft_original);

% Take the real part (to remove any small imaginary components due to numerical errors)
filtered_img = real(filtered_img);

% Display the filtered image
axes(handles.sortie);
imshow(uint8(filtered_img), []);

% Store the filtered image in handles for later use
handles.filtered_image = filtered_img;
guidata(hObject, handles);

% --------------------------------------------------------------------
function FiltrePasseHautDeButterworth_Callback(hObject, eventdata, handles)
% hObject    handle to FiltrePasseHautButterworth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Assuming handles.image_entree contains the input image
img = handles.image_entree;

% Convert the image to grayscale if it is RGB
if size(img, 3) == 3
    img = rgb2gray(img);
end

% Compute the 2D FFT of the image
fft_img = fft2(double(img));

% Shift the zero-frequency component to the center of the spectrum
fft_shifted = fftshift(fft_img);

% Get the size of the image
[M, N] = size(img);

% Create a Butterworth high-pass filter
D0 = 30; % Cutoff frequency (adjust as needed)
n = 2; % Order of the Butterworth filter (adjust as needed)
[u, v] = meshgrid(-N/2:N/2-1, -M/2:M/2-1);
D = sqrt(u.^2 + v.^2); % Distance from the center
H = 1 ./ (1 + (D0 ./ D).^(2*n)); % Butterworth high-pass filter

% Apply the filter to the shifted FFT
filtered_fft = fft_shifted .* H;

% Shift the zero-frequency component back to the original position
filtered_fft_original = ifftshift(filtered_fft);

% Compute the inverse FFT to get the filtered image
filtered_img = ifft2(filtered_fft_original);

% Take the real part (to remove any small imaginary components due to numerical errors)
filtered_img = real(filtered_img);

% Display the filtered image
axes(handles.sortie);
imshow(uint8(filtered_img), []);

% Store the filtered image in handles for later use
handles.filtered_image = filtered_img;
guidata(hObject, handles);


% --------------------------------------------------------------------
function FiltrePasseBandeIdeal_Callback(hObject, eventdata, handles)
% hObject    handle to FiltrePasseBandeIdeal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Assuming handles.image_entree contains the input image
img = handles.image_entree;

% Convert the image to grayscale if it is RGB
if size(img, 3) == 3
    img = rgb2gray(img);
end

% Compute the 2D FFT of the image
fft_img = fft2(double(img));

% Shift the zero-frequency component to the center of the spectrum
fft_shifted = fftshift(fft_img);

% Get the size of the image
[M, N] = size(img);

% Define the band-pass filter parameters
D0_low = 20; % Lower cutoff frequency (adjust as needed)
D0_high = 60; % Upper cutoff frequency (adjust as needed)

% Create a meshgrid for the frequency domain
[u, v] = meshgrid(-N/2:N/2-1, -M/2:M/2-1);
D = sqrt(u.^2 + v.^2); % Distance from the center

% Ideal band-pass filter
H = double(D >= D0_low & D <= D0_high); % Pass frequencies within the band

% Apply the filter to the shifted FFT
filtered_fft = fft_shifted .* H;

% Shift the zero-frequency component back to the original position
filtered_fft_original = ifftshift(filtered_fft);

% Compute the inverse FFT to get the filtered image
filtered_img = ifft2(filtered_fft_original);

% Take the real part (to remove any small imaginary components due to numerical errors)
filtered_img = real(filtered_img);

% Normalize the image to the range [0, 255]
filtered_img = uint8(255 * mat2gray(filtered_img));

% Display the filtered image
axes(handles.sortie);
imshow(filtered_img, []);
% Store the filtered image in handles for later use
handles.filtered_image = filtered_img;
guidata(hObject, handles);


% --------------------------------------------------------------------
function FiltrePasseBandedeButterworth_Callback(hObject, eventdata, handles)
% hObject    handle to FiltrePasseBandeButterworth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get input image
img = handles.image_entree;

% Convert to grayscale if RGB
if size(img, 3) == 3
    img = rgb2gray(img);
end

% Compute 2D FFT
fft_img = fft2(double(img));
fft_shifted = fftshift(fft_img);

% Get image dimensions
[M, N] = size(img);

% Define filter parameters
D0_low = 20;  % Lower cutoff frequency
D0_high = 60; % Upper cutoff frequency
n = 2;        % Filter order

% Create frequency domain meshgrid
[u, v] = meshgrid(-N/2:N/2-1, -M/2:M/2-1);
D = sqrt(u.^2 + v.^2);

% Create Butterworth filters according to course formula
H_lowpass = 1 ./ (1 + (D ./ D0_high).^(2*n));    % Low-pass with higher cutoff
H_highpass = 1 ./ (1 + (D0_low ./ D).^(2*n));    % High-pass with lower cutoff

% Combine to create band-pass filter
H = H_lowpass - H_highpass;

% Apply filter
filtered_fft = fft_shifted .* H;

% Inverse shift and compute IFFT
filtered_fft_original = ifftshift(filtered_fft);
filtered_img = ifft2(filtered_fft_original);

% Convert to real image and normalize
filtered_img = real(filtered_img);
filtered_img = uint8(255 * mat2gray(filtered_img));

% Display result
axes(handles.sortie);
imshow(filtered_img, []);

% Store filtered image
handles.filtered_image = filtered_img;
guidata(hObject, handles);


% --------------------------------------------------------------------
function FiltreSpectraleLocale_Callback(hObject, eventdata, handles)
% hObject    handle to FiltreSpectreLocale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Assuming handles.image_entree contains the input image
img = handles.image_entree;

% Convert the image to grayscale if it is RGB
if size(img, 3) == 3
    img = rgb2gray(img);
end

% Compute the 2D FFT of the image
fft_img = fft2(double(img));

% Shift the zero-frequency component to the center of the spectrum
fft_shifted = fftshift(fft_img);

% Get the size of the image
[M, N] = size(img);

% Define the local region for filtering
% Example: A rectangular region around the center
x_center = M/2; % Center row
y_center = N/2; % Center column
width = 50; % Width of the local region (adjust as needed)
height = 50; % Height of the local region (adjust as needed)

% Create a mask for the local region
mask = zeros(M, N);
mask(x_center - height/2:x_center + height/2, y_center - width/2:y_center + width/2) = 1;

% Apply the mask to the shifted FFT
local_fft = fft_shifted .* mask;

% Shift the zero-frequency component back to the original position
local_fft_original = ifftshift(local_fft);

% Compute the inverse FFT to get the filtered image
local_img = ifft2(local_fft_original);

% Take the real part (to remove any small imaginary components due to numerical errors)
local_img = real(local_img);

% Normalize the image to the range [0, 255]
local_img = uint8(255 * mat2gray(local_img));

% Display the filtered image
axes(handles.sortie);
imshow(local_img, []);

% Store the filtered image in handles for later use
handles.filtered_image = local_img;
guidata(hObject, handles);


% --------------------------------------------------------------------
function FiltreARejectiondeBande_Callback(hObject, eventdata, handles)
% hObject    handle to FiltreRejectionBande (voir GCBO)
% eventdata  réservé - à définir dans une future version de MATLAB
% handles    structure avec handles et données utilisateur (voir GUIDATA)

% Récupérer l'image d'entrée
img = handles.image_entree;

% Convertir en niveaux de gris si l'image est en couleur
if size(img, 3) == 3
    img = rgb2gray(img);
end

% Appliquer la transformée de Fourier 2D
fft_img = fft2(double(img));

% Déplacer les basses fréquences au centre du spectre
fft_shifted = fftshift(fft_img);

% Taille de l'image
[M, N] = size(img);

% Définition des paramètres du filtre
D0 = 40;  % Seuil de coupure centrale (ajuster selon besoin)
W = 20;   % Largeur de la bande rejetée (ajuster selon besoin)
D1 = D0 - W/2; % Borne inférieure
D2 = D0 + W/2; % Borne supérieure

% Création de la grille fréquentielle
[u, v] = meshgrid(-N/2:N/2-1, -M/2:M/2-1);
D = sqrt(u.^2 + v.^2); % Distance par rapport au centre

% Construction du masque de réjection de bande
H = ones(M, N);  % Initialise avec 1 (laisser passer toutes les fréquences)
H(D1 <= D & D <= D2) = 0; % Réjection des fréquences dans la bande [D1, D2]

% Appliquer le filtre sur le spectre fréquentiel
filtered_fft = fft_shifted .* H;

% Déplacer les basses fréquences à leur position d'origine
filtered_fft_original = ifftshift(filtered_fft);

% Appliquer la transformée de Fourier inverse
filtered_img = ifft2(filtered_fft_original);

% Prendre la partie réelle de l'image filtrée
filtered_img = real(filtered_img);

% Normaliser l'image dans la plage [0, 255]
filtered_img = uint8(255 * mat2gray(filtered_img));

% Afficher l'image filtrée
axes(handles.sortie);
imshow(filtered_img, []);

% Sauvegarder l'image filtrée dans handles
handles.filtered_image = filtered_img;
guidata(hObject, handles);


% --------------------------------------------------------------------
function FiltreHomomorphique_Callback(hObject, eventdata, handles)
% hObject    handle to FiltreHomomorphique (voir GCBO)
% eventdata  réservé - à définir dans une future version de MATLAB
% handles    structure avec handles et données utilisateur (voir GUIDATA)

% Récupérer l'image d'entrée
img = handles.image_entree;

% Convertir en niveaux de gris si l'image est en couleur
if size(img, 3) == 3
    img = rgb2gray(img);
end

% Convertir l'image en double pour éviter les erreurs de calcul
img = double(img) + 1; % Éviter log(0)

% Appliquer la transformation logarithmique
log_img = log(img);

% Appliquer la transformée de Fourier 2D
fft_img = fft2(log_img);

% Déplacer les basses fréquences au centre du spectre
fft_shifted = fftshift(fft_img);

% Obtenir la taille de l'image
[M, N] = size(img);

% Paramètres du filtre homomorphique
D0 = 50;   % Fréquence de coupure
n = 2;     % Ordre du filtre
gammaH = 2; % Amplification des hautes fréquences
gammaL = 0.5; % Atténuation des basses fréquences

% Création de la grille fréquentielle
[u, v] = meshgrid(-N/2:N/2-1, -M/2:M/2-1);
D = sqrt(u.^2 + v.^2); % Distance par rapport au centre

% Construction du filtre passe-haut de Butterworth
H = (gammaH - gammaL) * (1 - exp(-(D.^2) / (2 * D0^2))) + gammaL;

% Appliquer le filtre sur le spectre fréquentiel
filtered_fft = fft_shifted .* H;

% Déplacer les basses fréquences à leur position d'origine
filtered_fft_original = ifftshift(filtered_fft);

% Appliquer la transformée de Fourier inverse
filtered_log_img = ifft2(filtered_fft_original);

% Exponentiation pour retrouver l'image
filtered_img = exp(real(filtered_log_img));

% Normalisation de l'image
filtered_img = uint8(255 * mat2gray(filtered_img));

% Affichage du résultat
axes(handles.sortie);
imshow(filtered_img, []);

% Sauvegarde de l'image filtrée dans handles
handles.filtered_image = filtered_img;
guidata(hObject, handles);
