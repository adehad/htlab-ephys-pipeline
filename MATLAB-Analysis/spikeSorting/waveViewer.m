classdef waveViewer < handle
    % waveViewer is a handle object so its properties get passed by
    % reference instead of by value.
    properties (GetAccess = 'public', SetAccess = 'public')
        GUI;
    end
    
    methods(Access = public)
        %% Constructor
        function guiObj = waveViewer()
            %Initialize the GUI
            GUI_main(guiObj);
        end
    end
end


%% Graphical User Interface
function GUI_main(obj)
%Initializes the GUI elements in the .GUI property of the object 'obj'
% passed into the function


backgroundCol = [0.96 0.96 0.96];
obj.GUI.backgroundCol = backgroundCol; % store into obj.GUI
%% Create main figure window
obj.GUI.fig_main = figure( ...
    'units',             'pixels',...
    'position',          [100 100 1000 500],...
    'menubar',           'none',...
    'numbertitle',       'off',...
    'resize',            'on',...
    'name',              'Wave Viewer',...
    'Color',            backgroundCol,...
    'CloseRequestFcn',   {@myClose,obj});

%% Waveform Axis
obj.GUI.ax_Waveform = axes(...
        'Parent',           obj.GUI.fig_main,...
        'units',            'pixels',...
        'Position',         [0 0 450 450]+[65+500 50 -50 -50]);
title(obj.GUI.ax_Waveform, 'Selected Waveforms - Samples')
xlabel(obj.GUI.ax_Waveform, 'Samples' )
ylabel(obj.GUI.ax_Waveform, 'Amplitude')
    % obj.GUI.ax_Waveform.PlotBoxAspectRatio = [1 0.868632707774799 0.868632707774799];
    
%% Tab Group Initialising
obj.GUI.tabGroup_data = uitabgroup(...
    'Parent',           obj.GUI.fig_main,...
    'units',            'pixels', ...
    'Position',         [0 0 500 500]);
obj.GUI.tabSettings    = uitab('Parent', obj.GUI.tabGroup_data, 'Title', 'Settings ');

%% Tab: Settings
    %% Load Data Buttom
obj.GUI.btn_Load = uicontrol(...
    'style',        'pushbutton', ...
    'units',        'pixels', ...
    'parent',       obj.GUI.tabSettings, ...
    'position',     [20 460 60 22], ...
    'Callback',     {@buttonPush_Load,obj},...
    'String',       'Load Data');

    %% Unit Selection List Box
obj.GUI.listBox_unitList = uicontrol(...
    'style',            'listbox', ...
    'units',            'pixels', ...
    'parent',           obj.GUI.tabSettings, ...
    'position',         [380 280 100 200], ...
    'ForegroundColor',  [0 0 0], ...
    'String',           {});

    %% Unit Filter List Box
 obj.GUI.listBox_unitTypes = uicontrol(...
    'style',        'listbox', ...
    'units',        'pixels', ...
    'parent',       obj.GUI.tabSettings, ...
    'position',     [260 420 100 60], ...
    'String',       'Good|MUA|Unsorted|Noise', ...
    'Value',        [1,2,3,4], ... % default all are selected
    'Max',          4, ... % max number you can select
    'Visible',      1); 

    %% Plot Selected Unit Buttom
obj.GUI.btn_Plot = uicontrol(...
    'style',        'pushbutton', ...
    'units',        'pixels', ...
    'parent',       obj.GUI.tabSettings, ...
    'position',     [260 378 100 22], ...
    'Callback',     {@wavePlotter,obj,obj.GUI.ax_Waveform},...
    'String',       'Plot', ...
    'Visible',      0);

%% Export Data Buttons
obj.GUI.btn_Export = uicontrol(...
    'style',        'pushbutton', ...
    'units',        'pixels', ...
    'parent',       obj.GUI.tabSettings, ...
    'position',     [260 336 100 22], ...
    'Callback',     {@buttonPush_Export,obj},...
    'String',       'Export Data', ...
    'Visible',      0);


%% Create Input Fields & Labels
% Edit Fields - where the user can enter data
% Labels - Just labels
% Display Fields - where the current value is displayed
% -----------------------------------------------------------HIP--------%
% ---Labels
% obj.GUI.editFieldLabel_RHIP = uicontrol(...
%     'style',            'edit', ...
%     'units',            'pixels', ...
%     'parent',           obj.GUI.tabSettings, ...
%     'position',         [164 74 25 22], ...
%     'ForegroundColor',  [0 0 1], ...
%     'fontweight',       'bold', ...
%     'String',           'R', ...
%     'Enable',           'inactive', ...
%     'BackgroundColor',  backgroundCol);
    % custom java to remove border color
%     jEdit = findjobj(obj.GUI.editFieldLabel_RHIP);
    % jEdit.Border = [];

%% Time Range Control Boxes
% Time Range Field and Label
% obj.GUI.editField_timeRangeLL = uicontrol(...
%     'style',            'edit', ...
%     'units',            'pixels', ...
%     'parent',           obj.GUI.tabSettings, ...
%     'position',         [116 150 25 22], ...
%     'ForegroundColor',  [0 0 0], ...
%     'String',           '2', ...
%     'Callback',         {@obj.isNumber,-Inf,Inf,'timeRangeLL'});
% 
% obj.GUI.editField_timeRangeUL = uicontrol(...
%     'style',            'edit', ...
%     'units',            'pixels', ...
%     'parent',           obj.GUI.tabSettings, ...
%     'position',         [140 150 32 22], ...
%     'ForegroundColor',  [0 0 0], ...
%     'String',           '3', ...
%     'Callback',         {@obj.isNumber,-Inf,Inf,'timeRangeUL'});
% 
%     % Label
% obj.GUI.editFieldLabel_timeRange = uicontrol(...
%     'style',            'pushbutton', ...
%     'units',            'pixels', ...
%     'parent',           obj.GUI.tabSettings, ...
%     'position',         [14 150 98 22], ...
%     'ForegroundColor',  [0 0 0], ...
%     'String',           'Time Range (s)', ...
%     'BackgroundColor',  backgroundCol-0.01, ...
%     'Callback',         @obj.timeRangeSet);
%     % custom java to remove border color
%     % jEdit = findjobj(obj.GUI.editFieldLabel_timeRange);
%     % jEdit.Border = [];
%% Plot Handle Creation
% -------------------------------------------------------------- Angle ---
% HIP
% obj.GUI.hipAngL_animatedLine = animatedline(obj.GUI.ax_Waveform);
% obj.GUI.hipAngR_animatedLine = animatedline(obj.GUI.ax_Waveform);
% obj.GUI.hipAngL_animatedLine.Color = 'r';
% obj.GUI.hipAngR_animatedLine.Color = 'b';


%% Make Scalable
% Turns all 'units' properties to relative - so stuff can scale somewhat

for fieldN = fieldnames(obj.GUI)'
    try % if setting is already not normalized, set to normalized
        if ~strcmp(obj.GUI.(fieldN{1}).Units,'normalized') 
            obj.GUI.(fieldN{1}).Units = 'normalized';         
        end
    catch        % Do nothing if does not have Units property
%         obj.GUI.(fieldN{1}).Units
    end
end

%% Confirm Completion of GUI
display('GUI Setup done')
end


%% READ ME
% When calling a callback function it automatically sends handle for
% the source object which called the callback, and the event data,
% when we are calling a callback from another object i.e. the 'obj' it will
% additionally send the handle for that object.

% !!!! If we are editing any properties of the 'obj' class !!!!
% The below are callback templates that we can copy paste into the 'obj' class
% method section and change the input arguments in the obj class to: 
% (obj,src,eventdata, ... etc)
% Then in this GUI_main code we also need to change the callback to
% {@<Function Name>, ...}


%% Func: Number Verifier
function isNumber(src,eventdata,minVal,maxVal,var) % src = the handle to the source of the call 
  str = get(src, 'String'); % Correct way to get data 'edit' fields
  if isnan(str2double(str))     % is not a number?
      set(src, 'String','0');  % reset to 0
      warndlg('Input must be numerical'); % warn dialog
  end
  if (str2double(str))<maxVal && (str2double(str))>minVal  
      % do stuff if valid
      if size(var,1) == 1
        obj.(var) = str2double(str);
      else
        obj.(var(1))(str2double(var(2))) = str2double(str);
      end
  else
      set(src, 'String','0');  % reset to 0
      warndlg(strcat({'Input must be between '},num2str(minVal), ...
          {' & '}, num2str(maxVal),'.')); % warn dialog
  end
end

%% Func: Load Button Press Handler
function buttonPush_Load(src,eventdata,obj)
   [fileName, filePath] = uigetfile('*.mat','Select _extraced experiment .mat file:');

   load([filePath, fileName],'s');
   obj.GUI.s = s; % store loaded data into GUI Object - maybe store in new waveViewer property?
   % tried these but stuff was passed by value everytime
    %    evalin('base',['load(''',filePath, fileName,''',''s'')']);
    %     assignin('caller',obj.s,s);
    %     evalin('caller',['load(''',filePath, fileName,''',''s'')']);
    %     setappdata(obj,'s',s);
   
   % Update the Unit List Box : <UnitNumber> (<unitType>, <numberOfSpikes>)
   for ii=1:length(obj.GUI.s.clusters)
       clusterType = char(obj.GUI.s.cluster_groups{ii});
       clusterType = clusterType(1);
       newlistBox{ii} = string(strcat( ...
                           obj.GUI.s.clusters{ii}," (",...
                           clusterType,",", ...
                           num2str(size(obj.GUI.s.(['unit_',obj.GUI.s.clusters{ii}]),1)), ...
                           ")" ));
   end
   
   obj.GUI.listBox_unitList.String = newlistBox;
   obj.GUI.listBox_unitList.Max = length(obj.GUI.listBox_unitList.String);
   
   % Make stuff visible
   obj.GUI.btn_Plot.Visible = 1;
   obj.GUI.btn_Export.Visible = 1;

end

%% Func: Wave Plotter
function wavePlotter(src, eventdata,obj,axHandle)
    if length(obj.GUI.listBox_unitList.Value)>1
        hold(axHandle,'on');
    end
    
    % default "color order" of matlab
    colorScheme = [0,0.44700,0.74100;0.85000,0.32500,0.098000;0.92900,0.69400,0.12500;0.49400,0.18400,0.55600;0.46600,0.67400,0.18800;0.30100,0.74500,0.93300;0.63500,0.078000,0.18400];
    
   for ii=1:length(obj.GUI.listBox_unitList.Value)
        clusterNum = char(obj.GUI.listBox_unitList.String( ...
                            obj.GUI.listBox_unitList.Value(ii)) );
        clusterNum = clusterNum(1:(find(clusterNum==' ')-1));
        
        wavesToPlot = (obj.GUI.s.(['waves_',clusterNum])-mean(obj.GUI.s.(['waves_',clusterNum]),2));
        
        if ii>length(colorScheme)
            
        end
        plot( axHandle, ...
                wavesToPlot', ...
                'Color',colorScheme(ii-length(colorScheme)*floor(ii/length(colorScheme)),:) ...
             )
    end
    xlim([1 size(obj.GUI.s.(['waves_',clusterNum]),2)]);
    hold(axHandle,'off');
    drawnow
end

%% Func: Export Data Verifier
function buttonPush_Export(src,eventdata,obj)
    % Plots into a new figure so you can Save As or whatever
    figHan = figure;
    axHan = axes('Parent',figHan);
    wavePlotter(src,eventdata,obj,axHan);
end

%% Func: Toggle Visibility (Graphs)
function toggleGraphVisibility(src, eventdata, objH_graph,objH_graph2)
    
    if logical( (src.Value) )
       objH_graph.Visible = 1;
       objH_graph2.Visible = 1;
    else
        objH_graph.Visible = 0; 
        objH_graph2.Visible = 0;
    end

end

%% Func: Force Set of Time Range
function timeRangeSet(src,eventdata)

    axesList = { obj.GUI.ax_Waveform, ...
                 obj.GUI.figure_Angle_KNEE, ...
                 obj.GUI.figure_Angle_ANKLE, ...
                 ... % Velocity Plots
                 obj.GUI.figure_Velocity_HIP, ...
                 obj.GUI.figure_Velocity_KNEE, ...
                 obj.GUI.figure_Velocity_ANKLE, ...
        };

%             axesListEMG = {obj.GUI.figure_EMG}; % might have a different time axis

    for k=1:length(axesList)
        eval( strcat( "set(", (axesList{k}), ",'XLim',([obj.timeData(obj.inc)-obj.timeRangeLL, obj.timeData(obj.inc)+obj.timeRangeUL])" ) );
    end

%             for k=1:length(axesListEMG)
%                 eval( strcat( "set(", (axesListEMG{k}), ",'XLim',([obj.timeData(obj.inc)-obj.timeRangeLL, obj.timeData(obj.inc)+obj.timeRangeUL])" ) );
%             end

end
%% Func: Main Figure Close Destructor
function myClose(src,eventdata,obj)
    delete(src); %close the window
end
