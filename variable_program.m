
function varargout = variable_program(varargin)
% VARIABLE_PROGRAM MATLAB code for variable_program.fig
%      VARIABLE_PROGRAM, by itself, creates a new VARIABLE_PROGRAM or raises the existing
%      singleton*.
%
%      H = VARIABLE_PROGRAM returns the handle to a new VARIABLE_PROGRAM or the handle to
%      the existing singleton*.
%
%      VARIABLE_PROGRAM('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in VARIABLE_PROGRAM.M with the given input arguments.
%
%      VARIABLE_PROGRAM('Property','Value',...) creates a new VARIABLE_PROGRAM or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before variable_program_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to variable_program_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help variable_program

% Last Modified by GUIDE v2.5 14-Jun-2019 21:01:26

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @variable_program_OpeningFcn, ...
                   'gui_OutputFcn',  @variable_program_OutputFcn, ...
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



% --- Executes just before variable_program is made visible.
function variable_program_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to variable_program (see VARARGIN)

% Choose default command line output for variable_program
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes variable_program wait for user response (see UIRESUME)
% uiwait(handles.figure1);
landareas = shaperead('landareas.shp','UseGeoCoords',true);
[n,k]=size(landareas);
total = 1;
lon=zeros(1,20000);
lat=zeros(1,20000);
for j=1:n
    [a,b]=size(landareas(j).Lon);
    for m=1:b
        lon(1,total) = [landareas(j).Lon(1,m)] ;
        lat(1,total) = [landareas(j).Lat(1,m)] ;
        total=total +1;
    end
end

axesm ('variable_ortho', 'Frame', 'on', 'Grid', 'on','origin',[40 270 0],'rngz1',pi/3);   
geoshow(lat,lon);
%若需要填充，则直接使用下面一句，效率会降低
geoshow(landareas,'FaceColor',[0.5 1 0.5],'EdgeColor',[.6 .6 .6]);
tissot;
mdistort;

% --- Outputs from this function are returned to the command line.
function varargout = variable_program_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes during object creation, after setting all properties.
function axes1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes1


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% 确定时调用投影
cla reset
landareas = shaperead('landareas.shp','UseGeoCoords',true);
[n,k]=size(landareas);
total = 1;
lon=zeros(1,20000);
lat=zeros(1,20000);
for j=1:n
    [a,b]=size(landareas(j).Lon);
    for m=1:b
        lon(1,total) = [landareas(j).Lon(1,m)] ;
        lat(1,total) = [landareas(j).Lat(1,m)] ;
        total=total +1;
    end
end
% 重置坐标轴，防止图像重叠

global phi00 lambda00;
axesm ('variable_ortho', 'Frame', 'on', 'Grid', 'on','origin',[phi00 lambda00 0],'rngz1',[str2double(get(handles.edit1,'string'))],'zoom_factor',[str2double(get(handles.edit2,'string'))]);

geoshow(lat,lon);
%若需要填充，则直接使用下面一句，效率会降低
geoshow(landareas,'FaceColor',[0.5 1 0.5],'EdgeColor',[.6 .6 .6]);
tissot;
mdistort;
%测试参数
%gcm

% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global phi00 lambda00;
%由于放大镜的原理，改变了外围的正算公式，导致投影反解公式在放大镜外围失效，从而无法inputm选点，此bug待解决，不过不影响点击放大镜范围内的点的正常使用
%拾取经纬度坐标
[lat,lon]=inputm(1);
lon;

if ~isempty(lat) && ~isempty(lon) 
    phi00=lat;
    lambda00=lon;
else
    error(message('坐标反算失败'))
end

%set(handles.pushbutton3,'value',[phi00,lambda00]);

% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over pushbutton3.
function pushbutton3_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
