function gui_phaseplot(hObject,ParaTask) %#ok<INUSD>
% user interface for exploration of complex functions (internal use) 

% Usage: gui_phaseplot(hObject,ParaTask)
% internal use only, call CFEGUI instead

% with the assistence of Frank Martin 

% Part of the Complex Function Explorer (former Phase Plot) package 
% Version 1.1, February 1, 2014
% Copyright(c) 2012-2014 by Elias Wegert (elias@wegert.com, www.wegert.com)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% size of the menu and the top of a figure
width  = 350;
height = 480+30;
height_top_of_window = 60;

scrsz = get(0, 'ScreenSize');

XPosMain = (scrsz(3)-width)/2;
YPosMain = scrsz(4)-height-height_top_of_window;
YPosGui  = YPosMain;
XPosScm  = XPosMain+width+15;
XPosFct  = XPosMain-480-15;

LBnd = 15;

%%%%%%%%%%%% initiation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin <= 0
    
    %%%% main window user interface %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    HndlGui = findobj('Tag', 'PhasePlot');
    
    if isempty(HndlGui)
        HndlGui = figure('Name', 'Complex Function Explorer', ...
            'NumberTitle', 'off', ...
            'Position', [XPosMain, YPosGui, width, height], ...
            'Resize', 'off', ...
            'Tag', 'PhasePlot', ...
            'CloseRequestFcn', @gui_phaseplot);
    else
        figure(HndlGui)
        return
    end
    
	color = get(HndlGui,'Color');
        
    %%%%% window for phase plot of function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    hfig_phaseplot = findobj('Tag', 'fig_phaseplot');
    
    if isempty(hfig_phaseplot)
        figure('Name', 'Visualization of the function on its domain', ...
            'NumberTitle', 'off', ...
            'Position', [XPosFct, YPosMain, 480, 480], ...
            'Tag', 'fig_phaseplot');
    end
    
    %%%%% window for phase plot of the identity %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    hfig_identity = findobj('Tag','fig_identity');
    
    if isempty(hfig_identity)
        figure('Name', 'The color scheme of the image plane', ...
            'NumberTitle', 'off', ...
            'Position', [XPosScm, YPosMain, 480, 480], ...
            'Tag', 'fig_identity');
    end
    
    %%%%% building the user interface %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    figure(HndlGui)
    
    % axes for description of user interface
    
    axes('Parent', HndlGui, ...
        'Units', get(HndlGui, 'Units'), ...
        'Position', [0, 0, width, height], ...
        'XTick', [], 'YTick', [], ...
        'XColor', color, ...
        'YColor', color, ...
        'Box', 'on', ...
        'Color', color);
    
	% headline
    uicontrol(HndlGui, ...
        'Style', 'text', ...
        'Position', [LBnd, 0.93*height, width-2*LBnd, 0.04*height], ...
        'BackgroundColor', color, ...
        'String', 'VISUALIZATION OF COMPLEX FUNCTION', ...
        'FontSize', 12, ...
        'FontName', 'Helvetica');
    
    % function f(z)
    LHght = 0.85*height; 
    
    uicontrol(HndlGui, ...
        'Style', 'text', ...
        'Position', [LBnd, LHght, 50, 0.04*height], ...
        'String', 'f(z) =', ...
        'HorizontalAlignment','left', ...    
        'BackgroundColor', color, ...
        'FontName', 'Helvetica', ...
        'FontSize', 12);
    
    uicontrol(HndlGui, ...
        'Position', [LBnd+45, LHght, width-2*LBnd-45, 0.04*height], ...
        'Style', 'edit', ...
        'Tag', 'function', ...
        'String', 'MyFunction(z)', ...%(z^5-1)/(z^10+0.1)', ... 
        'BackgroundColor', 'w', ...
        'FontSize', 10, ...
        'FontName', 'Helvetica', ...
        'UserData',[], ...            
        'CallBack', ...
              'set(findobj(''Tag'',''function''),''UserData'',[]);');
        
    % specify domain of function
    LHght = 0.77*height; 
    
    HndlRect=uicontrol(HndlGui, ...
        'Style', 'text', ...
        'Tag','DomainRect', ...
        'Position', [LBnd, LHght, 200, 0.04*height], ...
        'HorizontalAlignment','left', ...    
        'String', 'defined in rectangle', ...
        'BackgroundColor', color, ...
        'FontName', 'Helvetica', ...
        'FontSize', 12, ...
        'CallBack', 'set(findobj(''tag'',''xmin''),''string'',''-2'');'); %#ok<NASGU>
    
    uicontrol(HndlGui, ...
        'Style', 'pushbutton', ...
        'Tag','DomainRect', ...
        'Position', [width-LBnd-153, LHght, 150, 0.04*height], ...
        'HorizontalAlignment','center', ...    
        'String', '(reset to default values)', ...
        'BackgroundColor', color, ...
        'FontName', 'Helvetica', ...
        'FontSize', 10, ...
        'CallBack', ...
        ['set(findobj(''tag'',''xmin''),''string'',''-2''); ',...
        'set(findobj(''tag'',''xmax''),''string'',''2''); ',...
        'set(findobj(''tag'',''ymin''),''string'',''-2''); ',...
        'set(findobj(''tag'',''ymax''),''string'',''2''); ']);
    
    uicontrol(HndlGui, ...
        'Style', 'text', ...
        'Tag','DomainSphere', ...
        'visible','off', ...
        'Position', [LBnd, LHght, 330, 0.04*height], ...
        'HorizontalAlignment','left', ...    
        'String', 'defined on Riemann sphere', ...
        'BackgroundColor', color, ...
        'FontName', 'Helvetica', ...
        'FontSize', 12);
       
    % control of x-range
    LHght = 0.7*height; 
    
    uicontrol(HndlGui, ...
        'Style', 'edit', ...
        'Tag', 'xmin', ...
        'Position', [LBnd, LHght, 40, 0.04*height], ...
        'HorizontalAlignment','center', ...    
        'String', '-2', ... 
        'UserData','-2',...
        'BackgroundColor', 'w', ...
        'FontSize', 10, ...
        'FontName', 'Helvetica', ...
        'CallBack','set(findobj(''Tag'',''function''),''UserData'',[]);');
              
    uicontrol(HndlGui, ...
        'Style', 'text', ...
        'tag','txtx', ...   
        'Position', [LBnd+43, LHght, 100, 0.04*height], ...
        'HorizontalAlignment','left', ...    
        'String', '< Re(z) <', ... 
        'BackgroundColor', color, ...
        'FontSize', 12, ...
        'FontName', 'Helvetica');
          
    uicontrol(HndlGui, ...
        'Style', 'edit', ...
        'Tag', 'xmax', ...
        'Position', [0.5*width-50, LHght, 40, 0.04*height], ...
        'HorizontalAlignment','center', ...    
        'String', '2', ... 
        'UserData','2',...
        'BackgroundColor', 'w', ...
        'FontSize', 10, ...
        'FontName', 'Helvetica', ...
        'CallBack', ...
              'set(findobj(''Tag'',''function''),''UserData'',[]);');
               
    % control of y-range
    LHght = 0.7*height; 
    
    uicontrol(HndlGui, ...
        'Style', 'edit', ...
        'Tag', 'ymin', ...
        'Position', [LBnd-3+width/2, LHght, 40, 0.04*height], ...
        'HorizontalAlignment','center', ...    
        'String', '-2', ... 
        'UserData','-2',...
        'BackgroundColor', 'w', ...
        'FontSize', 10, ...
        'FontName', 'Helvetica', ...
        'CallBack', ...
              'set(findobj(''Tag'',''function''),''UserData'',[]);');
       
    uicontrol(HndlGui, ...
        'Style', 'text', ...
        'tag','txty', ...
        'Position', [width-122, LHght, 100, 0.04*height], ...
        'HorizontalAlignment','left', ...    
        'String', '< Im(z) <', ... 
        'BackgroundColor', color, ...
        'FontSize', 12, ...
        'FontName', 'Helvetica');
    
    uicontrol(HndlGui, ...
        'Style', 'edit', ...
        'Tag', 'ymax', ...
        'Position', [width-57, LHght, 40, 0.04*height], ...
        'HorizontalAlignment','center', ...    
        'String', '2', ... 
        'UserData','2',...
        'BackgroundColor', 'w', ...
        'FontSize', 10, ...
        'FontName', 'Helvetica', ...
        'CallBack', ...
              'set(findobj(''Tag'',''function''),''UserData'',[]);');
        
   %%%% control of discretization (grid) 
    
    LHght = 0.6*height; 
    
    uicontrol(HndlGui, ...
        'Style', 'text', ...
        'tag','txtdis', ...
        'Position', [LBnd, LHght, 100, 0.04*height], ...
        'HorizontalAlignment','left', ...    
        'String', 'discretized at', ...
        'BackgroundColor', color, ...
        'FontName', 'Helvetica', ...
        'FontSize', 12);
    
    uicontrol(HndlGui, ...
        'Style', 'edit', ...
        'Tag', 'xres', ...
        'Position', [LBnd+111, LHght, 50, 0.04 * height], ...
        'String', '1000', ... 
        'BackgroundColor', 'w', ...
        'FontSize', 10, ...
        'FontName', 'Helvetica', ...
        'CallBack', ...
              ['if get(findobj(''tag'',''DomainType''),''Value'')==3, ', ...
              'XRES=get(findobj(''tag'',''xres''),''String''); ', ...
              'set(findobj(''tag'',''yres''),''String'',XRES); end, ', ...
              'set(findobj(''Tag'',''function''),''UserData'',[]);']);
        
    uicontrol(HndlGui, ...
        'Style', 'text', ...
        'tag','txttimes', ...
        'Position', [LBnd+163, LHght, 10, 0.04*height], ...
        'HorizontalAlignment','left', ...    
        'String', 'x', ...
        'BackgroundColor', color, ...
        'FontName', 'Helvetica', ...
        'FontSize', 12);
    
    uicontrol(HndlGui, ...
        'Style', 'edit', ...
        'Tag', 'yres', ...
        'Position', [LBnd+176, LHght, 50, 0.04 * height], ...
        'String', '1000', ... 
        'BackgroundColor', 'w', ...
        'FontSize', 10, ...
        'FontName', 'Helvetica', ...
        'CallBack', ...
         ['if get(findobj(''tag'',''DomainType''),''Value'')==3, ', ...
              'YRES=get(findobj(''tag'',''yres''),''String''); ', ...
              'set(findobj(''tag'',''xres''),''String'',YRES); end, ', ...
              'set(findobj(''Tag'',''function''),''UserData'',[]);']);
              
    uicontrol(HndlGui, ...
        'Style', 'text', ...
        'tag','TxtgridR', ...
        'Position', [LBnd+232, LHght, 100, 0.04*height], ...
        'HorizontalAlignment','left', ...    
        'String', 'grid points', ...
        'BackgroundColor', color, ...
        'FontName', 'Helvetica', ...
        'FontSize', 12);
    
     %%%% specification of color scheme %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
     LHght = 0.5*height; 
    
     uicontrol(HndlGui, ...
        'Style', 'text', ...
        'Position', [LBnd, LHght, 150, 0.04*height], ...
        'HorizontalAlignment','left', ...    
        'String', 'with color scheme', ...
        'BackgroundColor', color, ...
        'FontName', 'Helvetica', ...
        'FontSize', 12);
    
    XColCont = width/2-22;
    
    uicontrol(HndlGui, 'Style', 'popupmenu', ...
        'Tag', 'colorscheme', ...
        'Position', [XColCont, LHght, 182, 0.04 * height], ...
        'Value', 2, ...
        'String', {' plain phase plot',' phase & modulus',...
          ' phase + conformal grid',' standard domain coloring', ...
          ' enhanced domain coloring', ' polar chessboard ', ...
          ' cartesian chessboard', ' alternating b&w phase', ...
          ' alternating b&w modulus', ...
          ' b&w stripes (real part)', ' b&w stripes (imag part)'}, ... 
        'BackgroundColor', 'w', ...
        'FontSize', 10, ...
        'FontName', 'Helvetica', ...
        'UserData', 0, ...
        'Callback', 'gui_phaseplot(gcbo)');
    
    % set phase resolution
    
    LHght = 0.41*height;
    
    uicontrol(HndlGui, 'Style', 'slider', ...
        'Tag', 'PhaseResolution', ...
        'Position', [XColCont, LHght, 182, 0.04 * height], ...
        'BackgroundColor', color, ...
        'Min', 1, 'Max', 50, ...
        'Value', 20, ...
        'SliderStep', [1 5] / 49, ...
        'Callback', 'gui_phaseplot(gcbo)');
    
    LHght = 0.37*height;
    
    uicontrol(HndlGui, 'Style', 'text', ...
        'Tag', 'TxtPhaseResolution', ...
        'Position', [XColCont, LHght, 182, 0.04 * height], ...
        'BackgroundColor', color, ...
        'String','coarse       resolution       fine', ...
        'HorizontalAlignment','center', ...    
        'FontName', 'Helvetica', ...
        'FontSize', 10);
      
    
    % modify color scheme according to NIST coloring
    
    LHght = 0.32*height;
    
    uicontrol(HndlGui, 'Style', 'checkbox', ...
        'Tag', 'NISTcoloring', ...
        'Position', [XColCont, LHght, 182, 0.04 * height], ...
        'Value', 0, ...
        'String', 'NIST standard coloring', ... 
        'BackgroundColor', color, ...
        'FontSize', 12, ...
        'FontName', 'Helvetica', ...
        'Callback', 'gui_phaseplot(gcbo)');
     
    
    % display color scheme as icon
    
    ax_col = axes('Parent', HndlGui, ...
        'Units', get(HndlGui, 'Units'), ...
        'Position', [40, 0.33*height, 70, 70], ...
        'XTick', [], 'YTick', [], ...
        'XColor', color, ...
        'YColor', color, ...
        'Box', 'on', ...
        'Color', color);
    
        zz = zdomain(-1-1i,1+1i,120,120);
        PhasePlot(zz, zz, 'm');
    
    set(HndlGui, 'UserData', ax_col)
    
    % choose if domain is depicted in plane or on sphere %%%%%%%%%%%%%%%%%%
    
    LHght = 0.25*height;
    XSurfCont = 40+70;
    
    uicontrol(HndlGui, ...
        'Style', 'text', ...
        'Position', [LBnd, LHght, 150, 0.04*height], ...
        'HorizontalAlignment','left', ...    
        'String', 'plotted on', ...
        'BackgroundColor', color, ...
        'FontName', 'Helvetica', ...
        'FontSize', 12);
    
    uicontrol(HndlGui, 'Style', 'popupmenu', ...
        'Tag', 'DomainType', ...
        'Position', [XSurfCont, LHght, width-XSurfCont-LBnd, 0.04 * height], ...
        'Value', 1, ...
        'String', {' domain in complex plane',' domain on Riemann sphere',...
        ' entire Riemann sphere'}, ... 
        'BackgroundColor', 'w', ...
        'FontSize', 10, ...
        'FontName', 'Helvetica', ...
        'Callback', ...
            ['set(findobj(''Tag'',''function''),''UserData'',[]);', ...
            'gui_phaseplot(gcbo)']);
    
    % choose type of surface  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    LHght = 0.20*height;
    
    uicontrol(HndlGui, ...
        'Style', 'text', ...
        'Position', [LBnd, LHght, 150, 0.04*height], ...
        'HorizontalAlignment','left', ...    
        'String', 'as surface', ...
        'BackgroundColor', color, ...
        'FontName', 'Helvetica', ...
        'FontSize', 12);
  
   % the next uicontrol depends on the surface to be plotted on
      
     uicontrol(HndlGui, 'Style', 'popupmenu', ...
        'Tag', 'SurfaceTypeP', ...
        'Position', [XSurfCont, LHght, width-XSurfCont-LBnd, 0.04 * height], ...
        'Value', 1, ...
        'String', {' flat phase plot',' analytic landscape', ...
        ' logarithmic analytic landscape',' compressed analytic landscape'}, ... 
        'BackgroundColor', 'w', ...
        'FontSize', 10, ...
        'FontName', 'Helvetica', ...
        'Callback', 'gui_phaseplot(gcbo)',...
        'Visible','on');
          
    uicontrol(HndlGui, 'Style', 'popupmenu', ...
        'Tag', 'SurfaceTypeS', ...
        'Position', [XSurfCont, LHght, width-XSurfCont-LBnd, 0.04 * height], ...
        'Value', 1, ...
        'String', {' flat phase plot',' compressed analytic landscape'}, ... 
        'BackgroundColor', 'w', ...
        'FontSize', 10, ...
        'FontName', 'Helvetica', ...
        'Callback', 'gui_phaseplot(gcbo)',...
        'Visible','off');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
           
    % allow scaling of deformation for surface above sphere 
    LHght = 0.148*height;

    uicontrol(HndlGui, ...
        'Tag', 'TextDeformation', ...
        'Style', 'text', ...
        'Position', [LBnd, LHght, 150, 0.04*height], ...
        'String', 'deformation', ...
        'HorizontalAlignment','left', ...    
        'BackgroundColor', color, ...
        'FontName', 'Helvetica', ...
        'Visible','off', ...
        'FontSize', 12);
      
    uicontrol(HndlGui, 'Style', 'slider', ...
        'Tag', 'SliderDeformation', ...
        'Position', [XSurfCont, LHght, width-XSurfCont-LBnd, 0.04 * height],...
        'BackgroundColor', color, ...
        'Min', 0, 'Max', .95, ...
        'Value', 0.45, ...
        'Visible','off', ...
        'SliderStep', [1 5]/60, ...
        'Callback', 'gui_phaseplot(gcbo)');
    
    LHght = 0.11*height;
    
    uicontrol(HndlGui, 'Style', 'text', ...
        'Tag', 'TextDeformation2', ...
        'Position', [XSurfCont, LHght, width-XSurfCont-LBnd, 0.04 * height], ...
        'BackgroundColor', color, ...
        'String','small                                         large', ...
        'HorizontalAlignment','center', ...    
        'FontName', 'Helvetica', ...
        'Visible','off', ...
        'FontSize', 10);
     
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    
    LHght = 0.148*height;

    uicontrol(HndlGui, ...
        'Tag', 'HeightTxt1', ...
        'Style', 'text', ...
        'Position', [LBnd, LHght, 150, 0.04*height], ...
        'String', 'shown for ', ...
        'HorizontalAlignment','left', ...    
        'BackgroundColor', color, ...
        'FontName', 'Helvetica', ...
        'Visible','off', ...
        'FontSize', 12);
      
    uicontrol(HndlGui, 'Style', 'edit', ...
        'Tag', 'HeightLowBnd', ...
        'Position', [XSurfCont, LHght, 70, 0.04 * height],...
        'BackgroundColor','white', ...
        'String', '0', ...
        'HorizontalAlignment','center', ...    
        'FontName', 'Helvetica', ...
        'Visible','off', ...
        'FontSize', 10);
   
    uicontrol(HndlGui, ...
        'Tag', 'HeightTxt2', ...
        'Style', 'text', ...
        'Position', [XSurfCont+70, LHght, 84, 0.04*height], ...
        'String', ' < | f(z) | < ', ...
        'HorizontalAlignment','center', ...    
        'BackgroundColor', color, ...
        'FontName', 'Helvetica', ...
        'Visible','off', ...
        'FontSize', 12);
   
    uicontrol(HndlGui, 'Style', 'edit', ...
        'Tag', 'HeightUpBnd', ...
        'Position', [XSurfCont+155, LHght, 70, 0.04 * height],...
        'BackgroundColor', 'white', ...
        'String', '5.0', ...
        'HorizontalAlignment','center', ...    
        'FontName', 'Helvetica', ...
        'Visible','off', ...
        'FontSize', 10);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    % control button for depicting phase plot of the function
    uicontrol(HndlGui, 'Style', 'pushbutton',...
        'Tag', 'showphaseplot',...
        'Position', [LBnd, 0.05 * height, 80, 0.04 * height], ...
        'String', 'Show',...
        'FontSize', 12, ...
        'Callback', 'gui_phaseplot(gcbo)', ...
        'FontName', 'Helvetica');
    
    % control button for printing phase plot of the function in png format
    uicontrol(HndlGui, 'Style', 'pushbutton',...
        'Tag', 'printphaseplot',...
        'Position', [LBnd+81, 0.05 * height, 80, 0.04 * height], ...
        'String', 'Print',...
        'FontSize', 12, ...
        'Callback', 'gui_phaseplot(gcbo)', ...
        'FontName', 'Helvetica');
   
    % control button for saving interface data
    uicontrol(HndlGui, 'Style', 'pushbutton',...
        'Tag', 'SaveGUIData',...
        'Position', [LBnd+162, 0.05 * height, 80, 0.04 * height], ...
        'String', 'Save',...
        'FontSize', 12, ...
        'Callback', 'gui_phaseplot(gcbo)', ...
        'FontName', 'Helvetica');
   
    % control button for loading interface data
    uicontrol(HndlGui, 'Style', 'pushbutton',...
        'Tag', 'LoadGUIData',...
        'Position', [LBnd+243, 0.05 * height, 80, 0.04 * height], ...
        'String', 'Load',...
        'FontSize', 12, ...
        'Callback', 'gui_phaseplot(gcbo)', ...
        'FontName', 'Helvetica');
   
     set(HndlGui,'Position',[XPosMain,YPosMain,width,height]);   
    
	return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


tag = get(hObject,'Tag');
color_schemes = {'p', 'm', 'c', 'd', 'e', 'u', 'v', 'a', 'b', 'x', 'y'};

% react on change of color scheme

if strcmp(tag,'colorscheme')
   
    % main window
    HndlGui = findobj('Tag', 'PhasePlot');
    
    % axes to draw color scheme
    ax_col = get(HndlGui, 'UserData');
    
    % draw icon  
    val = get(hObject, 'Value');
    col = color_schemes{val};
    
    % modify if checkbox for NIST coloring is selected
    if get(findobj('Tag', 'NISTcoloring'), 'Value') >= 1
        col = cat(2, col, 'n');
    end
    
    % phase resolution
    pres = round(get(findobj('Tag', 'PhaseResolution'), 'Value'));
    
    % draw icon
    zz = zdomain(-1-1i,1+1i,120,120);
    axes(ax_col)
    PhasePlot(zz, zz, col, pres);

    % show or hide phaseresolution depending on color scheme
    
    CS = get(findobj('Tag','colorscheme'), 'Value');
       
    if CS==1 || CS==4
      set(findobj('Tag','PhaseResolution'),'Visible','off');  
      set(findobj('Tag','TxtPhaseResolution'),'Visible','off');  
    else
      set(findobj('Tag','PhaseResolution'),'Visible','on');
      set(findobj('Tag','TxtPhaseResolution'),'Visible','on');  
    end
    
    set(findobj('tag','colorscheme'),'UserData',0);
    
elseif strcmp(tag, 'NISTcoloring')
    
    % redraw color scheme icon according to NIST coloring or not
    gui_phaseplot(findobj('Tag', 'colorscheme'))
     
%%%%%%%%%% react on change of phase resolution %%%%%%%%%%%%%%%%%%%%%%%%%%%%

elseif strcmp(tag, 'PhaseResolution')
    
    % update the value
    
    val = get(hObject, 'Value');
    set(findobj('Tag', 'ValuePhaseResolution'), 'String', num2str(val))
    
    % redraw the color icon in user interface 
    gui_phaseplot(findobj('Tag', 'colorscheme'))

%%%%%%%%%% react on change of domain type %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

elseif strcmp(tag, 'DomainType')
    
    % determine type of domain
    hobj = findobj('Tag', 'DomainType');
    DomType = get(hobj, 'Value');
  
    % save actual values if DomType is changed to entire sphere
    XMIN=get(findobj('tag','xmin'),'string');
    
    if (DomType==3) && ~strcmp(XMIN,'-Inf')
        set(findobj('tag','xmin'),'UserData', ...
              get(findobj('tag','xmin'),'string'));
        set(findobj('tag','xmax'),'UserData', ...
              get(findobj('tag','xmax'),'string'));
        set(findobj('tag','ymin'),'UserData', ...
              get(findobj('tag','ymin'),'string'));
        set(findobj('tag','ymax'),'UserData', ...
              get(findobj('tag','ymax'),'string'));
    end
    
    % reload saved values if DomType changed from entire sphere to rectangle
    if DomType<3 && strcmp(XMIN,'-Inf')
        set(findobj('tag','xmin'),'string', ...
              get(findobj('tag','xmin'),'UserData'));
        set(findobj('tag','xmax'),'string', ...
              get(findobj('tag','xmax'),'UserData'));
        set(findobj('tag','ymin'),'string', ...
              get(findobj('tag','ymin'),'UserData'));
        set(findobj('tag','ymax'),'string', ...
              get(findobj('tag','ymax'),'UserData'));
    end
    
    STP = findobj('Tag','SurfaceTypeP');
    STS = findobj('Tag','SurfaceTypeS');
    SDS = findobj('Tag','SliderDeformation');
    TDS = findobj('Tag','TextDeformation');
    TDS2 = findobj('Tag','TextDeformation2');
    HSL = findobj('Tag','HeightLowBnd');
    HSU = findobj('Tag','HeightUpBnd');
    HT1 = findobj('Tag','HeightTxt1');
    HT2 = findobj('Tag','HeightTxt2');
      
    if DomType==1 % plane
        
        if strcmp(get(STP,'Visible'),'off')
            set(STP,'Visible','on');
            set(STS,'Visible','off');
            STPval = get(STP,'Value');
  
            if STPval>=2, 
                STSval=2; 
            else
                STSval=1;
                if STPval==1;
                  set(HSL,'Visible','off');
                  set(HSU,'Visible','off');
                  set(HT1,'Visible','off');
                  set(HT2,'Visible','off');
                else
                  set(HSL,'Visible','on');
                  set(HSU,'Visible','on');
                  set(HT1,'Visible','on');
                  set(HT2,'Visible','on');
                end
            end
            
            set(STS,'Value',STSval);
            set(SDS,'Visible','off');
            set(TDS,'Visible','off');
            set(TDS2,'Visible','off');
       end
        
    elseif DomType>=2 % sphere etc
        
        set(HSL,'Visible','off');
        set(HSU,'Visible','off');
        set(HT1,'Visible','off');
        set(HT2,'Visible','off');
            
        if strcmp(get(STS,'Visible'),'off')
            set(STS,'Visible','on');
            set(STP,'Visible','off');
            STSval = get(STS,'Value');
            set(STP,'Value',STSval);
            if get(STP,'Value')==2
              set(SDS,'Visible','on');
              set(TDS,'Visible','on');
              set(TDS2,'Visible','on');
            end
        end
        
        set(findobj('Tag','SurfaceTypeP'),'Visible','off')
        set(findobj('Tag','SurfaceTypeS'),'Visible','on')
      
    end

    % adapt grid generation
    
    if DomType==3 % entire sphere
    
        set(findobj('tag','xmin'),'string','-Inf');
        set(findobj('tag','xmax'),'string','+Inf');
        set(findobj('tag','ymin'),'string','-Inf');
        set(findobj('tag','ymax'),'string','+Inf');
        
        set(findobj('tag','DomainSphere'),'Visible','on');
        
    else %rectangle
        
        set(findobj('tag','DomainSphere'),'Visible','off');
        
    end
    
%%%%%%%%% react on change of surface type %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

elseif strcmp(tag, 'SurfaceTypeP') % surface when domain is flat

    SDS = findobj('Tag', 'SliderDeformation');
    TDS = findobj('Tag', 'TextDeformation');
    TDS2 = findobj('Tag', 'TextDeformation2');
        
    set(SDS,'Visible','off');
    set(TDS,'Visible','off');
    set(TDS2,'Visible','off');
   
    % hide or show cutoff options for analytic landscape 
    
    HSL = findobj('Tag','HeightLowBnd');
    HSU = findobj('Tag','HeightUpBnd');
    HT1 = findobj('Tag','HeightTxt1');
    HT2 = findobj('Tag','HeightTxt2');
      
    STP = get(findobj('Tag','SurfaceTypeP'),'Value');
    
    if STP==1 || STP ==4 
      set(HSL,'Visible','off');
      set(HSU,'Visible','off');
      set(HT1,'Visible','off');
      set(HT2,'Visible','off');
    elseif STP==2 % ana land
      set(HSL,'Visible','on');
      set(HSU,'Visible','on');
      set(HT1,'Visible','on');
      set(HT2,'Visible','on');
      set(HT2,'String',' < | f(z) | < ');
      set(HSL,'String','0');
      set(HSU,'String','5.0');
    elseif STP==3 % log ana land
      set(HSL,'Visible','on');
      set(HSU,'Visible','on');
      set(HT1,'Visible','on');
      set(HT2,'Visible','on');
      set(HT2,'String','< log |f(z)| <');
      set(HSL,'String','-5.0');
      set(HSU,'String','5.0');
    end
    
elseif strcmp(tag, 'SurfaceTypeS') % surface when domain is spherical
   
     HSL = findobj('Tag','HeightLowBnd');
     HSU = findobj('Tag','HeightUpBnd');
     HT1 = findobj('Tag','HeightTxt1');
     HT2 = findobj('Tag','HeightTxt2');
     
     set(HSL,'Visible','off');
     set(HSU,'Visible','off');
     set(HT1,'Visible','off');
     set(HT2,'Visible','off');
       
     SDS = findobj('Tag', 'SliderDeformation');
     TDS = findobj('Tag', 'TextDeformation');
     TDS2 = findobj('Tag', 'TextDeformation2');
      
    
    if get(findobj('Tag', 'SurfaceTypeS'),'Value')==2;
    
      set(SDS,'Visible','on');
      set(TDS,'Visible','on');
      set(TDS2,'Visible','on');
      
    else
    
      set(SDS,'Visible','off');
      set(TDS,'Visible','off');
      set(TDS2,'Visible','off');
      
     
    end
    
%%%%%%%%% show phase plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

elseif strcmp(tag, 'showphaseplot')
 
   % show phase plot
   hfig_phaseplot = findobj('Tag', 'fig_phaseplot');
     
    % check if values f(z) must be computed (updated)
    
    % read data from user interface
    xmin = eval(get(findobj('Tag', 'xmin'), 'String'));
    xmax = eval(get(findobj('Tag', 'xmax'), 'String'));
    xres = eval(get(findobj('Tag', 'xres'), 'String'));
    ymin = eval(get(findobj('Tag', 'ymin'), 'String'));
    ymax = eval(get(findobj('Tag', 'ymax'), 'String'));
    yres = eval(get(findobj('Tag', 'yres'), 'String'));
        
    if isempty(get(findobj('Tag','function'),'UserData'))
         
      % argument z of function f(z) to be depicted
      if get(findobj('Tag', 'DomainType'), 'Value') == 3
        % phase plot on entire sphere
        z = zplanePP(xres); 
      else
        % phase plot on rectangle
        z = zdomain(xmin+1i*ymin, xmax+1i*ymax, xres, yres);
      end
      
      % allow function definitions depending on x=real(z) and y=imag(z)
      x = real(z); %#ok<NASGU>
      y = imag(z); %#ok<NASGU>
    
      % compute function values
    
      % read definition of f as string from interface
      fct = get(findobj('Tag', 'function'), 'String');
    
      % modify string to allow simpler syntax for array operations
      fct = strrep(fct,'.*','*');  
      fct = strrep(fct,'./','/');  
      fct = strrep(fct,'.^','^');  
    
      % save the simplified symbolic representation
      fsym = fct;
      set(hfig_phaseplot, 'UserData', fsym)
         
      % change everything to array operations
      fct = strrep(fct,'*','.*');
      fct = strrep(fct,'/','./');
      fct = strrep(fct,'^','.^');
    
      % evaluate function values
      fct = eval(fct);
      
      if size(fct)==1, fct=fct*ones(size(z)); end
      
      % save values of z and f(z) in the figure
      set(findobj('Tag', 'function'), 'UserData',{z, fct});
         
   else
       
      valofzandf = get(findobj('Tag','function'),'UserData');
      z = valofzandf{1};  fct =valofzandf{2};
      
   end 
   
      % color scheme
      hobj = findobj('Tag', 'colorscheme');
      val  = get(hobj, 'Value');
      col  = color_schemes{val};

      % modify color scheme for NIST coloring 
      if get(findobj('Tag', 'NISTcoloring'), 'Value') >= 1
         col = cat(2,col,'n');
      end
    
      % determine type of domain
      hobj = findobj('Tag', 'DomainType');
      DomType = get(hobj, 'Value');
    
      % determine type of surface
      STP=findobj('Tag','SurfaceTypeP');
      STS=findobj('Tag','SurfaceTypeS');
    
      if strcmp(get(STP,'Visible'),'on')
       SurfType = get(STP, 'Value');
      elseif strcmp(get(STS,'Visible'),'on')
        SurfType = get(STS, 'Value');
      end
    
      % phase resolution
      pres = round(get(findobj('Tag', 'PhaseResolution'), 'Value'));
      
      %%%%%%%%%%%%%%%%%%%%%% create figure %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      if isempty(hfig_phaseplot)
       
        hfig_phaseplot = figure('Name', ...
            'Visualization of the function on its domain', ...
            'NumberTitle', 'off', ...
            'Position', [XPosFct, YPosMain, 480, 480], ...
            'Tag', 'fig_phaseplot');
        
      end
    
      figure(hfig_phaseplot)
    
      % load symbolic representation of function
      fsym = get(hfig_phaseplot, 'UserData');
          
      % determine type of representation (surface) 
    
      if DomType >= 2
        % phase plot on sphere
        
        if SurfType>=2 
          %compressed analytic landscape on sphere
            
          deform = get(findobj('Tag','SliderDeformation'),'Value');
          PPOnSphere(z, fct, col, pres,[],deform);
                      
        else
          PPOnSphere(z, fct, col, pres,[],0);
        
        end
        
          TitFct=title(['f(z)=',fsym,' on Riemann sphere']);
          set(TitFct,'Interpreter','none','FontSize',10)
         
          view(50,20)
          
          lghtdull
        
       elseif DomType==1 % phase plot on flat domain
           
         if SurfType==1  % phase plot
            azel = [0,90];    
            bestlight = '';
            PhasePlot(z, fct, col, pres);
            
         elseif SurfType==2  % analytic landscape
            azel = [40,20];
            bestlight = 'lghtdull';
             
            hmin = str2double(get(findobj('Tag','HeightLowBnd'),'String'));
            hmax = str2double(get(findobj('Tag','HeightUpBnd'),'String'));
            
            if hmax<0, hmax = abs(hmax); end
            if hmin<0, hmin=0; end
            if hmin>hmax, hmaxs=hmax; hmax=hmin; hmin=hmaxs; end
            if hmin==hmax, hmax=hmin+1; end
            
            set(findobj('Tag','HeightLowBnd'),'String',num2str(hmin));
            set(findobj('Tag','HeightUpBnd'),'String',num2str(hmax));
            
            fcthght=abs(fct); 
            small = (fcthght<hmin);
            large = (fcthght>hmax);
            fcthght(large)=hmax;
            fcthght(small)=hmin;
            hmin = min(min(fcthght)); hmax = max(max(fcthght));
            fcthght = (1/(hmax-hmin))*(fcthght-hmin);
           
            PhasePlot(z, fct, col, pres, [], fcthght);
                     
            % simpler call with internal computation of height
            %PhasePlot(z, fct, [col,'a'], pres);
        
         elseif SurfType==3  % logarithmic analytic landscape
            azel = [40,20];
            bestlight = 'lghtdull';
         
            hmin = str2double(get(findobj('Tag','HeightLowBnd'),'String'));
            hmax = str2double(get(findobj('Tag','HeightUpBnd'),'String'));
            
            if hmin>hmax, hmaxs=hmax; hmax=hmin; hmin=hmaxs; end
            if hmin==hmax, hmax=hmin+1; end
            
            set(findobj('Tag','HeightLowBnd'),'String',num2str(hmin));
            set(findobj('Tag','HeightUpBnd'),'String',num2str(hmax));
            
            absf = abs(fct);
            fcthght = log(absf);
            small = (fcthght<hmin);
            large = (fcthght>hmax);
            fcthght(small) = hmin;
            fcthght(large) = hmax;
            hmin = min(min(fcthght)); hmax = max(max(fcthght));
            fcthght = (1/(hmax-hmin))*(fcthght-hmin);
                 
            PhasePlot(z, fct, col, pres, [], fcthght);
        
            % simpler call with internal computation of height
            %PhasePlot(z, fct, [col,'l'], pres);
        
         elseif SurfType==4  % compressed analytic landscape
            azel = [40,20];
            bestlight = 'lghtdull';
            
            fcthght=(2/pi)*atan(abs(fct)); %hmin=0; hmax=1;
            
            PhasePlot(z, fct, col, pres, [], fcthght);
            
            % simpler call with internal computation of height
            %PhasePlot(z, fct, [col,'c'], pres);
           
         end
         
         view(azel);
       
         TitFct=title(['f(z)=',fsym,' in ', ...
           num2str(xmin),'<Re z<',num2str(xmax),', ', ...
           num2str(ymin),'<Im z<',num2str(ymax)]);
         
         set(TitFct,'Interpreter','none','FontSize',10)
         
         eval(bestlight)
       
      end
 
    % redraw phase plot of identity (color scheme) in separate window
    
    hfig_identity = findobj('Tag', 'fig_identity');
    
    if isempty(hfig_identity)
        hfig_identity = ...
            figure('Name', 'The color scheme of the image plane', ...
            'NumberTitle', 'off', ...
            'Position', [(scrsz(3)+width)/2+15, YPosMain, 480, 480], ...
            'Tag', 'fig_identity');
    end
    
    if get(findobj('tag','colorscheme'),'UserData')==0
    
      XResi = min(xres,500); YResi = min(yres,500);
      z = zdomain(-1-1i,1+1i,XResi,YResi); 
    
      figure(hfig_identity)
      PhasePlot(z,z,col,pres);
      TitClm=title('color scheme in square -1<Re z<1, -1<Im z<1');
      set(TitClm,'Interpreter','none','FontSize',10)
 
      % activate phase plot of function and user interface
      figure(hfig_phaseplot)
      figure(findobj('tag','PhasePlot'))
      
      set(findobj('tag','colorscheme'),'UserData',1);

    end
    
%%%%% print phaseplot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
elseif strcmp(tag, 'printphaseplot')
    
    CurDir = cd;
    
    if isdir('Figures')==0; mkdir('Figures'), end
    
    if isunix
        SavePath=[CurDir,'/Figures/'];
    else
        SavePath=[CurDir,'\Figures\'];
    end
    
    % print phase plot to png file
    [file, pathname] = uiputfile([SavePath,'*.png'], ...
        'Export current visualization');
     
    if file == 0
        % no file was specified; do nothing
        return
    else
        file = [pathname, file];
        
        % make the figure of the phase plot the current figure
        hfig_phaseplot = findobj('Tag', 'fig_phaseplot');
        
        if isempty(hfig_phaseplot)
           return
        end
            
        figure(hfig_phaseplot)
        
        % send phase plot to png file
        print(gcf, '-r300', '-dpng', file);
        
        % for higher resolution activate next line 
        % print(gcf, '-r600', '-dpng', file);
      
    end

    figure(findobj('tag','PhasePlot'))
    
%%%%%%%% save phase plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% save all parameters from interface to recompute phase plot

elseif strcmp(tag, 'SaveGUIData')
    
    % read parameters from user interface
    
    FctStr    = get(findobj('Tag','function'),'String');
    XMin      = get(findobj('Tag','xmin'),'String');
    XMax      = get(findobj('Tag','xmax'),'String');
    YMin      = get(findobj('Tag','ymin'),'String');
    YMax      = get(findobj('Tag','ymax'),'String');
    XResStr   = get(findobj('Tag','xres'),'String');
    YResStr   = get(findobj('Tag','yres'),'String');
    ColScVal  = get(findobj('Tag','colorscheme'),'Value');
    PhsResVal = get(findobj('Tag','PhaseResolution'),'Value');
    NISTVal   = get(findobj('Tag','NISTcoloring'),'Value');
    DomType   = get(findobj('Tag','DomainType'),'Value');
    SurfTypeP = get(findobj('Tag','SurfaceTypeP'),'Value');
    SurfTypeS = get(findobj('Tag','SurfaceTypeS'),'Value');
    HghtScal  = get(findobj('Tag','SliderDeformation'),'Value');
    HghtMax   = get(findobj('Tag','HeightUpBnd'),'String');
    HghtMin   = get(findobj('Tag','HeightLowBnd'),'String');
    
    % build structure to save
    
    PPGuiData = struct('FctStr',FctStr, ...
    'XMin', XMin, 'XMax', XMax, 'YMin', YMin, 'YMax', YMax, ...
    'XResStr', XResStr, 'YResStr', YResStr, ...
    'ColScVal', ColScVal, 'PhsResVal', PhsResVal, 'NISTVal', NISTVal, ...
    'DomType', DomType, 'SurfTypeP', SurfTypeP, 'SurfTypeS', SurfTypeS, ...
    'HghtScal', HghtScal,'HghtMin', HghtMin,'HghtMax', HghtMax); %#ok<NASGU>
    
    % choose file to save structure
 
    if ~isdir('Examples')
        mkdir('Examples')
    end
    
    CurDir = cd;  
    
    if isunix
        SaveFile = [CurDir,'/Examples/Example.mat'];
    else
        SaveFile = [CurDir,'\Examples\Example.mat'];
    end
    
    uisave('PPGuiData',SaveFile);
    
    figure(findobj('tag','PhasePlot'))
        
%%%%%%%% load phase plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

elseif strcmp(tag, 'LoadGUIData')

    % load parameters from file to recompute phase plot
    
    if ~isdir('Examples')
        mkdir('Examples')
    end
        
    CurDir = cd;  
    
    if isunix
       LoadFile = [CurDir,'/Examples/Example.mat']; 
    else
       LoadFile = [CurDir,'\Examples\Example.mat'];
    end
    
    [ExampleFile,ExamplePath] = uigetfile(LoadFile);
    
    if ExampleFile==0
       return
    end
    
    % this might not work in former MATLAB versions
    PPGuiData = load([ExamplePath,ExampleFile]);
    
    PPGuiData = PPGuiData.PPGuiData;
           
    set(findobj('tag','function'),'string',PPGuiData.FctStr);
    set(findobj('tag','function'),'UserData',[]);
    set(findobj('tag','xmin'),'string',PPGuiData.XMin);
    set(findobj('tag','xmax'),'string',PPGuiData.XMax);
    set(findobj('tag','ymin'),'string',PPGuiData.YMin);
    set(findobj('tag','ymax'),'string',PPGuiData.YMax);
    set(findobj('Tag','xres'),'String',PPGuiData.XResStr);
    set(findobj('Tag','yres'),'String',PPGuiData.YResStr);
    set(findobj('Tag','colorscheme'),'Value',PPGuiData.ColScVal);
    set(findobj('Tag','PhaseResolution'),'Value',PPGuiData.PhsResVal);
    set(findobj('Tag','NISTcoloring'),'Value',PPGuiData.NISTVal);
    set(findobj('Tag','DomainType'),'Value',PPGuiData.DomType);
    set(findobj('Tag','SurfaceTypeP'),'Value',PPGuiData.SurfTypeP);
    set(findobj('Tag','SurfaceTypeS'),'Value',PPGuiData.SurfTypeS);
    set(findobj('Tag','SliderDeformation'),'Value',PPGuiData.HghtScal);
    set(findobj('Tag','HeightLowBnd'),'String',PPGuiData.HghtMin);
    set(findobj('Tag','HeightUpBnd'),'String',PPGuiData.HghtMax);
    set(findobj('Tag','colorscheme'),'UserData',0);
    
    PPDomType=PPGuiData.DomType;
    PPSurfTypeS=PPGuiData.SurfTypeS;
    PPSurfTypeP=PPGuiData.SurfTypeP;
    
    if PPDomType==1 %plane
       set(findobj('Tag', 'SurfaceTypeS'),'Visible','off');
       set(findobj('Tag', 'SurfaceTypeP'),'Visible','on');
    elseif PPDomType>=2 %sphere
       set(findobj('Tag', 'SurfaceTypeS'),'Visible','on');
       set(findobj('Tag', 'SurfaceTypeP'),'Visible','off');
    end
          
    if PPDomType~=1 && PPSurfTypeS==2 %sphere
       set(findobj('Tag', 'SliderDeformation'),'Visible','on');
       set(findobj('Tag', 'TextDeformation'),'Visible','on');
       set(findobj('Tag', 'TextDeformation2'),'Visible','on');
    else
       set(findobj('Tag', 'SliderDeformation'),'Visible','off');
       set(findobj('Tag', 'TextDeformation'),'Visible','off');
       set(findobj('Tag', 'TextDeformation2'),'Visible','off');
    end
    
    if PPDomType==1 && PPSurfTypeP==2 %analytic landscape
       set(findobj('Tag', 'HeightTxt1'),'Visible','on');
       set(findobj('Tag', 'HeightLowBnd'),'Visible','on');
       set(findobj('Tag', 'HeightTxt2'),'Visible','on');
       set(findobj('Tag', 'HeightTxt2'),'String',' < | f(z) | < ');
       set(findobj('Tag', 'HeightUpBnd'),'Visible','on');
    elseif PPDomType==1 && PPSurfTypeP==3 % log ana land
       set(findobj('Tag', 'HeightTxt1'),'Visible','on');
       set(findobj('Tag', 'HeightLowBnd'),'Visible','on');
       set(findobj('Tag', 'HeightTxt2'),'Visible','on');
       set(findobj('Tag', 'HeightTxt2'),'String','< log |f(z)| <');
       set(findobj('Tag', 'HeightUpBnd'),'Visible','on');
    else
       set(findobj('Tag', 'HeightTxt1'),'Visible','off');
       set(findobj('Tag', 'HeightLowBnd'),'Visible','off');
       set(findobj('Tag', 'HeightTxt2'),'Visible','off');
       set(findobj('Tag', 'HeightUpBnd'),'Visible','off');
    end
    
    gui_phaseplot(findobj('tag','DomainType'));
    gui_phaseplot(findobj('tag','colorscheme'));
    
    figure(findobj('tag','PhasePlot'))
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
elseif strcmp(tag, 'PhasePlot')
   
    % close the figure for displaying (if not already done)
    hfig_phaseplot = findobj('Tag', 'fig_phaseplot');

    if ~isempty(hfig_phaseplot)
        delete(hfig_phaseplot)
    end
    
    hfig_identity = findobj('Tag', 'fig_identity');
    
    if ~isempty(hfig_identity)
        delete(hfig_identity)
    end
    
    % close main window
    delete(hObject)

end

end

% internal function to control light, can be adapted

function lghtdull

    light
    material dull
 
    %set(gcf,'Renderer','zbuffer')
    %set(findobj(gca,'type','surface'),...
    %   'FaceLighting','phong',...
    %   'AmbientStrength',.7,'DiffuseStrength',.3,...
    %   'SpecularStrength',.5,'SpecularExponent',25,...
    %   'BackFaceLighting','unlit')

 
end