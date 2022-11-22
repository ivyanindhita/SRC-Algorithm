function varargout = FaceRecognitionTool_TA(varargin)
% FACERECOGNITIONTOOL_TA MATLAB code for FaceRecognitionTool_TA.fig
%      FACERECOGNITIONTOOL_TA, by itself, creates a new FACERECOGNITIONTOOL_TA or raises the existing
%      singleton*.
%
%      H = FACERECOGNITIONTOOL_TA returns the handle to a new FACERECOGNITIONTOOL_TA or the handle to
%      the existing singleton*.
%
%      FACERECOGNITIONTOOL_TA('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FACERECOGNITIONTOOL_TA.M with the given input arguments.
%
%      FACERECOGNITIONTOOL_TA('Property','Value',...) creates a new FACERECOGNITIONTOOL_TA or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before FaceRecognitionTool_TA_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to FaceRecognitionTool_TA_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help FaceRecognitionTool_TA

% Last Modified by GUIDE v2.5 09-Mar-2022 10:55:03

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @FaceRecognitionTool_TA_OpeningFcn, ...
                   'gui_OutputFcn',  @FaceRecognitionTool_TA_OutputFcn, ...
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


% --- Executes just before FaceRecognitionTool_TA is made visible.
function FaceRecognitionTool_TA_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to FaceRecognitionTool_TA (see VARARGIN)
set(handles.text18,'visible','off');
set(handles.edit12,'visible','off');
set(handles.text16,'visible','off');
set(handles.text17,'visible','off');
set(handles.uipanel10,'visible','off');
set(handles.uipanel9,'visible','off');
axes(handles.axes4)
cla
axes(handles.axes3)
cla
% Choose default command line output for FaceRecognitionTool_TA
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes FaceRecognitionTool_TA wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = FaceRecognitionTool_TA_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function edit8_Callback(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit8 as text
%        str2double(get(hObject,'String')) returns contents of edit8 as a double


% --- Executes during object creation, after setting all properties.
function edit8_CreateFcn(hObject, eventdata, ~)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles) 
%-----------------------------------------------------------------------
%                           TRAIN PUSHBUTTON                           %
%-----------------------------------------------------------------------

% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global A m1 n1 Mphi Nphi Phi1 Phi2 Phi3 Phi4 No_Files_In_Class_Folder Class_Count Training_Set_Folder  Extract_Choice

if ismac
    Database_Path = '/TA/FaceRecogTool_nia/FaceDatabase2/Database2'; %database asli 
    %Database_Path = '/TA/FaceRecogTool_nia/FaceDatabase/Database1'; %modifdataset
elseif ispc
    Database_Path = 'D:\TA\FaceRecogTool_nia\FaceDatabase2\Database2'; %database asli
    %Database_Path = 'D:\TA\FaceRecogTool_nia\FaceDatabase\Database1'; %modifdataset
end
Listing_Training_Set_Folder = [uigetdir(Database_Path)];
% --- fixed error on cancel
if Listing_Training_Set_Folder == 0
    return
end

Extract_Choice = str2num(handles.edit13.String)
m1 = 8; %Merubah ukuran dimensi citra
n1 = 6;


%--Set The visible of
set(handles.text18,'visible','off');
set(handles.edit12,'visible','off');
set(handles.text16,'visible','off');
set(handles.text17,'visible','off');
set(handles.uipanel10,'visible','off');
set(handles.uipanel9,'visible','off');
set(handles.axes4,'visible','off');
axes(handles.axes4)
cla
set(handles.uipanel9,'visible','off');
axes(handles.axes3)
cla

% CS Preparation
        Mphi = m1*n1;
        Nphi = 112*92;
        
        %--------------------------------------------NOTES------------------------------------------------------------------------------%
        %(rand = random integers) *rand = decimals starts at 0
        %sign =  returns an array Y the same size as X, where each element of Y is, 1 = lebih dari 0, 0 = 0, -1 = kurang dari 0
        %ones = Create an array of all ones
        %Gaussian / normal = randn (Terdapat titik yang maksimum)
        %Uniform = rand (Rata atau setara)
        
        %-------------------------------------------------------------------------------------------------------------------------------%
        %Cara mengecek Gaussian dan Uniform, 1. X1 = Phi1(1,:) 2. lihat grafik > hist(X1,100) (begitupun sebaliknya untuk phi yang lain)
        %----------- KODING PENURUNAN DIMENSI DENGAN 4 TEKNIK -------------
        
        %--CS1
        Phi1 = (sign(randn(Mphi,Nphi))+ones(Mphi,Nphi))/2; %Random Uniform Binary
        
        %--CS2
        Phi2 = randn(Mphi,Nphi,'double'); %Random Gaussian Integers
        
        %--CS3
        Phi3 = rand(Mphi,Nphi); %Random Uniform Integers
        
        
        %--CS4
        MaxPhi4 = 256;
        %Phi4 = randi(MaxPhi4,Mphi,Nphi);
        TempPhi4 = randi(MaxPhi4,Mphi,Nphi);
        Phi4 = TempPhi4/MaxPhi4;
        
                     
if ispc
    Training_Set_Folder = [Listing_Training_Set_Folder '\'];
elseif ismac
    Training_Set_Folder = [Listing_Training_Set_Folder '/'];
end
TS_Vector = dir(Training_Set_Folder);
No_Folders_In_Training_Set_Folder = length(TS_Vector);
File_Count = 1;
Class_Count = 1;
h = waitbar(0,'Reading Test Images,Please wait...');
tic
for k = 3:No_Folders_In_Training_Set_Folder
    waitbar(k/(No_Folders_In_Training_Set_Folder-2))
    if ispc
        Class_Folder = [Training_Set_Folder '\' TS_Vector(k).name,'\'];
    elseif ismac
        Class_Folder = [Training_Set_Folder '/' TS_Vector(k).name,'/'];
    end
    CF_Tensor = dir(Class_Folder);
    No_Files_In_Class_Folder(Class_Count) = length(CF_Tensor)-2;
    strr = sprintf('Reading Test Images...!, # of Classes = %d, Now Reading %d ',No_Folders_In_Training_Set_Folder-2,Class_Count);
    
    drawnow;
    for p = 3:No_Files_In_Class_Folder(Class_Count)+2
        Tmp_Image_Path = Class_Folder; %image training dataset (6-10)
        Tmp_Image_Name = CF_Tensor(p).name;
        Tmp_Image_Path_Name = [Tmp_Image_Path,Tmp_Image_Name];
        if strcmp(Tmp_Image_Name,'Thumbs.db')
            break
        end
        
        %----Capture PIC (TRAINING) 
        CapPic = imread(Tmp_Image_Path_Name); %dataset dibaca 112 x 92 piksel
        CapPic = double(CapPic); %dataset dalam bentuk desimal
        if Extract_Choice >= 3 
        CapPic = CapPic/norm(CapPic(:)); % normalized image (agar komponen desimal ,menjadi 0)
        CapPic = CapPic - mean(CapPic(:)); % zero mean image
        else
        end
        
        %----Show Image Original Training Samples Pic
        
        axes(handles.axes3) 
        imshow(Tmp_Image_Path_Name);  %image training dataset (6-10)
        realimgsize = size(CapPic); %ukuran asli citra 
        title(['Trained Image Pixels: ' strrep(num2str(realimgsize),"",'x')],'Color','white','FontSize',12); %mengambil dimensi pada keluaran "112x92"
        drawnow; %menampilkan tulisan trained image 112x92
        set(handles.edit8,'string',['Class Trained: ' TS_Vector(k).name ' | ' 'Image Trained: ' Tmp_Image_Name]); %memunculkan class trained (kotak sebelah kiri train)
        
        % RGB to Gray conversion
        if length(size(CapPic))==3
            Tmp_Image = rgb2gray(CapPic); %merubah dataset rgb menjadi greyscale
        else
            Tmp_Image = CapPic; %dataset real sudah greyscale
        end
        
        % Feature Extraction pada CapPic
        if Extract_Choice == 1
        % ---1.Down Sampled
        Tmp_Image_Down_Sampled = double(imresize(Tmp_Image,[m1 n1])); 
        Tmp_Image_Extraction = Tmp_Image_Down_Sampled(:); % Choice 1
        %menurunkan dimensi citra dengan downscale (imresize merupakan kodingan utk ds) sesuai dengan faktor reduksi + scanning
        %-------------------------------------------------------------------------------------------------------------
        
        %-------------------------------------------------NOTES-----------------------------------------------------------------------%
        %Cara mengecek teknik scanning (lakukan debug, pada command window 1. ketik CapPic, 2. ketik Tmp_Image(:,1); (mengambil row 1)%
        %-----------------------------------------------------------------------------------------------------------------------------%
        elseif  Extract_Choice == 2 %Random Uniform Binary
        % ---3. Compressive Sampling
        % --- CS1
        Xphi = Tmp_Image(:); %scanning column to row scanning (kebawah)
        Tmp_Image_CS1 = Phi1*Xphi; 
        Tmp_Image_Extraction = Tmp_Image_CS1; % Choice 3a
        
        %sizetmpimagecs1 = size(Tmp_Image_CS1) 
        %---------------------------------------------------------------------------------------%
        elseif Extract_Choice == 3 %Random Gaussian
        % ---CS2
        Xphi = Tmp_Image(:); %scanning column to row
        Tmp_Image_CS2 = Phi2*Xphi; 
        Tmp_Image_Extraction = Tmp_Image_CS2; % Choice 3b
        %---------------------------------------------------------------------------------------%
        elseif Extract_Choice == 4 %Random Uniform Integers
        %---CS3 
        Xphi = Tmp_Image(:); %scanning column to row
        Tmp_Image_CS3 = Phi3*Xphi;
        Tmp_Image_Extraction = Tmp_Image_CS3; % Choice 3c
        else
        %---------------------------------------------------------------------------------------%
        %---CS4
        Xphi = Tmp_Image(:);
        Tmp_Image_CS4 = Phi4*Xphi;
        Tmp_Image_Extraction = Tmp_Image_CS4; % Choice 3d
        end        
        % Image Matrix Choice
        Image_Data_Matrix(:,File_Count) = Tmp_Image_Extraction; %count dataset selanjutnya, kembali ke awal line 197 !
        %sizeimagedatmat = size(Image_Data_Matrix)
        File_Count = File_Count+1;
    end
    Class_Count = Class_Count+1;

end
close(h)
A = Image_Data_Matrix; %A berisi dataset yang telah diturunkan dengan extract choice yang digunakan
A = A/(diag(sqrt(diag(A'*A))));
timeElapsed = toc %waktu
%profile report
set(handles.text18,'visible','on');
set(handles.text17,'visible','on');
handles.text18.String = [num2str(timeElapsed) ' secs'];
axes(handles.axes3)
title('');
cla
set(handles.edit8,'string','');

% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
%----------------------------------------------------------------------------------------%
%                                  SINGLE TESTING
%----------------------------------------------------------------------------------------% 

% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global A m1 n1 Phi1 Phi2 Phi3 No_Files_In_Class_Folder Class_Count Training_Set_Folder Total_Test Extract_Choice
set(handles.text18,'visible','on');
set(handles.edit11,'visible','on');
set(handles.uipanel9,'visible','on');


% [Test_File Test_File_Path] = uigetfile('*.jpg;*.pgm;*.png;*.tif','Select a Test Image');
% --- where database live
if ismac
    Database_Path = '/TA/FaceRecogTool_nia/FaceDatabase/Database1'; %dataset oklusi
    %Database_Path = '/TA/FaceRecogTool_nia/FaceDatabase2/Database2'; %realdataset
    Database_Path = [Database_Path '/'];
elseif ispc
    Database_Path = 'D:\TA\FaceRecogTool_nia\FaceDatabase\Database1'; %dataset oklusi
    %Database_Path = 'D:\TA\FaceRecogTool_nia\FaceDatabase2\Database2'; %real dataset
    Database_Path = [Database_Path '\'];
end
[Test_File Test_File_Path] = uigetfile([Database_Path '*.jpg;*.pgm;*.png;*.tif'],'Select a Test Image');
if Test_File == 0
    return
end
%profile on
tic 

% --- file to test
test_image_path = [Test_File_Path Test_File];
% --- folder class identification begin
str_test_image_path = test_image_path;
[testfilepath] = fileparts(str_test_image_path);
if ispc
    idx = strfind(testfilepath,'\');
elseif ismac
    idx = strfind(testfilepath,'/');
end
testfoldername = testfilepath(idx(end)+1:end)
% --- folder class identification end
axes(handles.axes3)
cla
axes(handles.axes4)
cla
axes(handles.axes3)
imshow(test_image_path);
title('Test Image','Color','red','FontSize',15);
drawnow;
Test_File = [Test_File_Path Test_File];

% ---- Capture Pic
CapPic = imread(Test_File);
CapPic = double(CapPic);
CapPic = CapPic/norm(CapPic(:)); % normalized image
CapPic = CapPic - mean(CapPic(:)); % zero mean image

if length(size(CapPic))==3
    Test_Image = rgb2gray(CapPic);
else
    Test_Image = CapPic;
end
                % FEATURE EXTRACTION
                if Extract_Choice == 1
                %---1 Down Sampled
                Test_Image_Down_Sampled = imresize(Test_Image,[m1 n1]); % --- down sampling
                Test_Image_Extraction = double(Test_Image_Down_Sampled(:)); %choice 1
                %elseif Extract_Choice == 2       
                elseif Extract_Choice == 2  %Random Uniform Binary
                %---3 Compressive Sampling
                Xphi = Test_Image(:);
                %-----CS1 Choice 3a           
                Test_Image_CS1 = Phi1*Xphi; 
                Test_Image_Extraction = Test_Image_CS1; %choice 3a
                elseif Extract_Choice == 3 %Random Gaussian
                %-----CS2 Choice 3b %Random Uniform binary
                Xphi = Test_Image(:);
                Test_Image_CS2 = Phi2*Xphi;
                Test_Image_Extraction = Test_Image_CS2; %choice 3b
                else
                Xphi = Test_Image(:);
                Test_Image_CS3 = Phi3*Xphi;
                Test_Image_Extraction = Test_Image_CS3; %choice 3c   
                end
                
y = Test_Image_Extraction; %berisi dataset yang telah diturunkan dimensinya
n = size(A,2);
f=ones(2*n,1); %Create an array of all ones
Aeq=[A -A]; %nilai A diganti dengan -A
lb=zeros(2*n,1);
% options = optimoptions('linprog','Display','none');
% x1 = linprog(f,[],[],Aeq,y,lb,[],options);

x1 = linprog(f,[],[],Aeq,y,lb,[],[],[]); %lasso, tmpt y juga
x1 = x1(1:n)-x1(n+1:2*n);
nn = No_Files_In_Class_Folder
nn = cumsum(nn);
tmp_var = 0;
k1 = Class_Count-1;
Total_Test = 0;
for i = 1:k1
    i;
    delta_xi = zeros(length(x1),1);
    if i == 1
        delta_xi(1:nn(i)) = x1(1:nn(i));
    else
        tmp_var = tmp_var + nn(i-1);
        begs = nn(i-1)+1;
        ends = nn(i);
        delta_xi(begs:ends) = x1(begs:ends);
    end
    tmp(i) = norm(y-A*delta_xi,2);
    tmp1(i) = norm(delta_xi,1)/norm(x1,1);
    % figure(1);
    % histogram(tmp(i),i);
    Total_Test = Total_Test + i;
end
Sparse_Conc_Index = (k1*max(tmp1)-1)/(k1-1); % --- Sparsity Concentration Index (SCI)
clss = find(tmp==min(tmp)) % --- train image class detection with minimum reconstruction error (kelas)
cccc = dir([Training_Set_Folder]);
% cccc(clss+2).name
if ispc
    Which_Folder = dir([Training_Set_Folder,cccc(clss+2).name,'\']);
elseif ismac
    Which_Folder = dir([Training_Set_Folder,cccc(clss+2).name,'/']);
end
Which_Image = randsample(3:length(Which_Folder),1);
if ispc
    Image_Path = [Training_Set_Folder,cccc(clss+2).name,'\',Which_Folder(Which_Image).name];
elseif ismac
    Image_Path = [Training_Set_Folder,cccc(clss+2).name,'/',Which_Folder(Which_Image).name];
end
Class_Image = (Image_Path); % --- show image file
%x = imread(Class_Image) (melihat bit greyscale)
% --- folder class identification begin
str_train_image_path = Image_Path;
[trainfilepath] = fileparts(str_train_image_path);
if ispc
    idx = strfind(trainfilepath,'\');
elseif ismac
    idx = strfind(trainfilepath,'/');
end
trainfoldername = trainfilepath(idx(end)+1:end)
% --- folder class identification end
axes(handles.axes4);
imshow(Class_Image)
title('Detected Image','Color','black','FontSize',15)
% --- compare detected image class versus original source class
if trainfoldername == testfoldername
   disp('Match!')
   handles.text12.ForegroundColor = 'blue';
   handles.text12.String = 'Match';
else
   disp('Failed!')
   handles.text12.ForegroundColor = 'red';
   handles.text12.String = 'Failed';
end
set(handles.edit8,'string',['Class Test: ' testfoldername ' | ' 'Class Detected: ' trainfoldername]);
timeElapsed = toc
%profile report
handles.text18.String = [num2str(timeElapsed) ' secs'];
handles.edit11.String = num2str(Total_Test);



% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles) 
%----------------------------------------------------------------------------------------%
%                                  MULTI TESTING
%----------------------------------------------------------------------------------------% 
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global A m1 n1 Nphi Mphi Phi1 Phi2 Phi3 Phi4 No_Files_In_Class_Folder Class_Count Training_Set_Folder Extract_Choice
set(handles.text18,'visible','off');
set(handles.edit12,'visible','off');
set(handles.text16,'visible','off');
set(handles.text17,'visible','off');
set(handles.uipanel10,'visible','off');
set(handles.uipanel9,'visible','off');

% FullPath=uigetdir('');
if ismac
   % Database_Path = '/TA/FaceRecogTool_nia/FaceDatabase/Database1'; %modif dataset oklusi
   Database_Path = '/TA/FaceRecogTool_nia/FaceDatabase2/Database2'; %real dataset
elseif ispc
    %Database_Path = 'D:\TA\FaceRecogTool_nia\FaceDatabase\Database1'; %modif dataset oklusi
  Database_Path = 'D:\TA\FaceRecogTool_nia\FaceDatabase2\Database2'; %real dataset
end
FullPath = [uigetdir(Database_Path)];
% --- fixed error on cancel
if FullPath == 0
    return
end
%profile on
tic
IsTrue=0;
IsFalse=0;
TotImg=0;
TestFiles=dir(FullPath);
set(handles.uipanel9,'visible','on');
for k=1:length(TestFiles)
   % k
  if ~strcmp(TestFiles(k,1).name(1),'.')
       if ispc
            Imgfiles=dir([FullPath '\' TestFiles(k).name])
        elseif ismac
            Imgfiles=dir([FullPath '/' TestFiles(k).name]);
        end
        for m=1:length(Imgfiles)
        %    m
            if ~strcmp(Imgfiles(m,1).name(1),'.')
                axes(handles.axes3)
                if ispc
                    Test_File = [FullPath '\' TestFiles(k,1).name '\' Imgfiles(m,1).name]
                elseif ismac
                    Test_File = [FullPath '/' TestFiles(k,1).name '/' Imgfiles(m,1).name];
                end
                % --- file to test
                test_image_path = [Test_File];
                % --- folder class identification begin
                str_test_image_path = test_image_path;
                [testfilepath] = fileparts(str_test_image_path);
                if ispc
                    idx = strfind(testfilepath,'\')
                elseif ismac
                    idx = strfind(testfilepath,'/');
                end
                testfoldername = testfilepath(idx(end)+1:end);
                
                % --- folder class identification end
                if ispc
                    set(handles.edit8,'string',[TestFiles(k,1).name '\' Imgfiles(m,1).name])
                elseif ismac
                    set(handles.edit8,'string',[TestFiles(k,1).name '/' Imgfiles(m,1).name]);
                end
                % show image original
                imshow(Test_File);
                title('Test Image','Color','White','FontSize',12);

                drawnow;
               
                %----Capture Pic
                CapPic = imread(Test_File);
                CapPic = double(CapPic);
                if Extract_Choice >= 3
                CapPic = CapPic/norm(CapPic(:)); % normalized image
                CapPic = CapPic - mean(CapPic(:)); % zero mean image
                else
                end
                
                if length(size(CapPic))==3
                    Test_Image = rgb2gray(CapPic);
                else
                    Test_Image = CapPic;
                end
                % FEATURE EXTRACTION
                if Extract_Choice == 1
                %---1 Down Sampled
                Test_Image_Down_Sampled = imresize(Test_Image,[m1 n1]); % --- down sampling
                Test_Image_Extraction = double(Test_Image_Down_Sampled(:)); %choice 1
                      
                elseif Extract_Choice == 2 
                %---3 Compressive Sampling
                Xphi = Test_Image(:);
                %-----CS1 Choice 3a           
                Test_Image_CS1 = Phi1*Xphi; %dikali dengan test image 
                Test_Image_Extraction = Test_Image_CS1; %choice 3a
                elseif Extract_Choice == 3
                %-----CS2 Choice 3b
                Xphi = Test_Image(:);
                Test_Image_CS2 = Phi2*Xphi;
                Test_Image_Extraction = Test_Image_CS2; %choice 3b
                elseif Extract_Choice == 4
                %-----CS3 Choice 3c
                Xphi = Test_Image(:);
                Test_Image_CS3 = Phi3*Xphi;
                Test_Image_Extraction = Test_Image_CS3; %choice 3c 
                else
                %-----CS4 Choice 3d
                Xphi = Test_Image(:);
                Test_Image_CS4 = Phi4*Xphi;
                Test_Image_Extraction = Test_Image_CS4; %choice 3d 
                end
               
                % IMAGE EXTRACTION
                y = Test_Image_Extraction;
                n = size(A,2);
                f = ones(2*n,1);
                Aeq = [A -A];
                lb = zeros(2*n,1);
                x1 = linprog(f,[],[],Aeq,y,lb,[],[],[]); %LASSO 
                x1 = x1(1:n)-x1(n+1:2*n);
                nn = No_Files_In_Class_Folder
                nn = cumsum(nn);
                tmp_var = 0;
                k1 = Class_Count-1;
                for i = 1:k1
                    delta_xi = zeros(length(x1),1);
                    if i == 1
                        delta_xi(1:nn(i)) = x1(1:nn(i));
                    else
                        tmp_var = tmp_var + nn(i-1);
                        begs = nn(i-1)+1;
                        ends = nn(i);
                        delta_xi(begs:ends) = x1(begs:ends);
                    end
                    tmp(i) = norm(y-A*delta_xi,2); %residuals
                    tmp1(i) = norm(delta_xi,1)/norm(x1,1);
                end
                TotImg=TotImg+1; %proses mencari kelas dari x
                Sparse_Conc_Index(TotImg) = (k1*max(tmp1)-1)/(k1-1);
                clss = find(tmp==min(tmp));
                % figure,plot(tmp)
                ssttrr = sprintf('The Test Image Corresponds to Class: %d',clss);
                cccc = dir([Training_Set_Folder]);
                if ispc
                    Which_Folder = dir([Training_Set_Folder,cccc(clss+2).name,'\'])
                elseif ismac
                    Which_Folder = dir([Training_Set_Folder,cccc(clss+2).name,'/']);
                end
                details(Which_Folder)
                length(Which_Folder)
                Which_Image = randsample(3:length(Which_Folder),1)
                if ispc
                    Image_Path = [Training_Set_Folder,cccc(clss+2).name,'\',Which_Folder(Which_Image).name]
                elseif ismac
                    Image_Path = [Training_Set_Folder,cccc(clss+2).name,'/',Which_Folder(Which_Image).name];
                end
                Class_Image = (Image_Path);
                
                % --- folder class identification begin
                str_train_image_path = Image_Path;
                [trainfilepath] = fileparts(str_train_image_path);
                if ispc
                    idx = strfind(trainfilepath,'\');
                elseif ismac
                    idx = strfind(trainfilepath,'/');
                end
                trainfoldername = trainfilepath(idx(end)+1:end);
                % --- folder class identification end
                % show detected image
                axes(handles.axes4);
                imshow(Class_Image);
                
                title('Detected Image','Color','White','FontSize',12);
                %pause(1)

                % --- compare detected image class versus original source class
                if trainfoldername == testfoldername
                   disp('Match!');
                   %purplecol = handles.edit2.ForegroundColor
                   handles.text12.ForegroundColor = handles.edit12.ForegroundColor;
                   handles.text12.String = 'Match';
                   IsTrue=IsTrue+1;
                else
                   disp('Failed!');
                   handles.text12.ForegroundColor = 'red';
                   handles.text12.String = 'Failed';
                   IsFalse=IsFalse+1;
                end
                                
                %axes(handles.axes4)
                %
                eta = (IsTrue/TotImg)*100
                set(handles.edit12,'visible','on');
                set(handles.text16,'visible','on');
                set(handles.edit12,'String',[num2str(eta) '%']);
                drawnow;
                timeElapsed = toc;
                
                set(handles.text18,'visible','on');
                set(handles.text17,'visible','on');
                handles.text18.String = [num2str(timeElapsed) ' secs'];
                set(handles.uipanel10,'visible','on');
                handles.edit9.String = num2str(IsTrue);
                handles.edit10.String = num2str(IsFalse);
                handles.edit11.String = num2str(TotImg);

            end
        end
    end
end
%profile report
%if eta >= 96 && Extract_Choice == 3
 %   writematrix(Phi2, 'Phi2Test.xlsx');
  %  writematrix(Phi2, 'Phi2.xlsx');
%writematrix(Phi2,'Phi2.txt');
%writematrix(Phi2,'Phi2_tab.txt','Delimiter','tab')
%writematrix(Phi2,'Phi2.xls')
%else
 %   writematrix(Phi2, 'Phi2Test.xlsx');
%end
%tempphi = readmatrix('Phi2.txt');
%tempphi = readmatrix('Phi2.xlsx');
%writematrix(Phi2,'Phi_tab.txt','Delimiter','tab')
%writematrix(Phi2,'Phi.xls')
%}
%{
eta = (IsTrue/TotImg)*100;
set(handles.edit2,'visible','on');
set(handles.text4,'visible','on');
set(handles.edit2,'String',[num2str(eta) '%']);
drawnow;
timeElapsed = toc
set(handles.text3,'visible','on');
set(handles.text6,'visible','on');
handles.text3.String = [num2str(timeElapsed) ' secs'];
set(handles.uipanel7,'visible','on');
handles.edit4.String = num2str(IsTrue);
handles.edit5.String = num2str(IsFalse);
handles.edit3.String = num2str(TotImg);
%}


function edit9_Callback(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit9 as text
%        str2double(get(hObject,'String')) returns contents of edit9 as a double


% --- Executes during object creation, after setting all properties.
function edit9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit10_Callback(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit10 as text
%        str2double(get(hObject,'String')) returns contents of edit10 as a double


% --- Executes during object creation, after setting all properties.
function edit10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit11_Callback(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit11 as text
%        str2double(get(hObject,'String')) returns contents of edit11 as a double


% --- Executes during object creation, after setting all properties.
function edit11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit12_Callback(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit12 as text
%        str2double(get(hObject,'String')) returns contents of edit12 as a double


% --- Executes during object creation, after setting all properties.
function edit12_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit13_Callback(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit13 as text
%        str2double(get(hObject,'String')) returns contents of edit13 as a double


% --- Executes during object creation, after setting all properties.
function edit13_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
