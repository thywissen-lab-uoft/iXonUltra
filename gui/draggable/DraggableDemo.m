function DraggableDemo(demotitle, mode)
%{
This function is a demo for the draggable.m function and should be distributed
along with it. It basically presents some of draggable.m"s features: users of
draggable.m are invited to read dragdemo"s source.

>> dragdemo(demotitle, mode);

runs the demo which title is contained in the string demotitle. Mode determines
whether to use axes() and figure() or uiaxes() and uifigure().


demotitle:

"crosshair" - A draggable crosshair is drawn; its movements are limited so that
    its center never leaves the axes. This is the default demo.

"dragtest" - Some graphical objects (such as rectangles, lines, patches and
    plots) are drawn, following various constraints. For more details, please
    read the comments in the source code.

"snapgrid" - An example on how to use draggable"s "motionfcn" argument so that a
    draggable object snaps to a grid.

"polymove" - A polygon with draggable vertices is drawn; draggable"s "motionfcn"
    argument is used to redraw the polygon each time a vertex is moved.

"sliders" - A set of vertical sliders are created; a single call to draggable is
    used for all sliders since they share the same parameters. A "motionfcn"
    argument is used to update the displayed slider values.


mode:

"axes" - Draw demo in axes() and figure() objects.

"uiaxes" - Draw demo in uiaxes() and uifigure() objects.

(C) Copyright 2004-2020
François Bouffard
fbouffard@gmail.com

(C) Copyright 2020
William Warriner
wwarriner@gmail.com

%}

AXES = "axes";
UIAXES = "uiaxes";

if nargin < 1
    demotitle = "crosshair";
end

if nargin < 2
    mode = AXES;
end

switch lower(mode)
    case AXES
        figure_fn = @figure;
        axes_fn = @axes;
    case UIAXES
        figure_fn = @uifigure;
        axes_fn = @uiaxes;
    otherwise
        error("Unknown mode.");
end

switch lower(demotitle)
    case "crosshair"
        % We first build the cross, using plot() with a NaN y-element so that it
        % is non-continuous. The cross spans 0.2 units:
        f = figure_fn();
        
        ax = axes_fn(f);
        axis(ax, "equal");
        set(ax, "DataAspectRatio", [1 1 1], "Xlim", [0 1], "YLim", [0 1]);
        title(ax, "CROSSHAIR demo for draggable.m");
        hold(ax, "on");
        
        p = plot(ax, [0.5 0.5 0.4 0.4 0.6], [0.4 0.6 NaN 0.5 0.5]);
        set(p, "Color", "k", "LineWidth", 2);
        
        % Creating a draggable graphics object, using no movement constraint and
        % limits corresponding to the axes limits, plus half the size of the
        % cross, so that its center always stays in the axes:
        dg = draggable(p);
        dg.xlim = [-0.1 1.1];
        dg.ylim = [-0.1 1.1];
    case "dragtest"
        % This demo is used in draggable development to test some of its
        % features. Users interested in modifying draggable.m should ensure that
        % it works with this demo in order to verify backward compatibility.
        
        % Figure creation; we set up an initial WindowButtonDownFcn in order to
        % test proper figure properties recovery after the object is dragged:
        f = figure_fn();
        wbd_fn = @(varargin)disp("Window Click!");
        set(f, "WindowButtonDownFcn", wbd_fn);
        
        ax = axes_fn(f);
        axis(ax, "equal");
        box(ax, "on");
        set(ax, "DataAspectRatio", [1 1 1], "XLim", [0.3 1], "YLim", [0 1] );
        title(ax, "DRAGTEST demo for draggable.m");
        hold(ax, "on");
        
        % The Green Rectangle tests default behavior
        greenrect = rectangle(ax, "Position", [0.7 0.1 0.1 0.2]);
        set(greenrect, "FaceColor", "g");
        draggable(greenrect);
        
        % The Red Square tests right mouse button only, it shouldn't be
        % draggable by left or middle mouse buttons.
        redsquare = rectangle(ax, "Position", [0.5 0.5 0.1 0.1]);
        set(redsquare, "FaceColor", "r");
        dr = draggable(redsquare);
        dr.button = draggable.RIGHT_BUTTON;
        
        % The Blue Rectangle will demonstrate a case in which the object cannot
        % be dragged past the axis limits on the right, left and bottom.
        % However, the object can be dragged past the axis limit on top;
        % furthermore, the object starts off-limits to the left, so that it can
        % be dragged inside the limits, but cannot be dragged back off-limits
        % afterwards.
        bluerect = rectangle(ax, "Position", [0.2 0.4 0.2 0.2]);
        set(bluerect, "FaceColor", "b");
        draggable(bluerect, "none", [0.3 1 0 inf]);
        
        % The Magenta Line will demonstrate a line object being dragged
        % with horizontal movement constraint with default parameters.
        magline = line(ax, [0.7 0.9], [0.7 0.9]);
        set(magline, "Color", "m", "LineWidth", 2);
        draggable(magline, "h");
        
        % The Cyan Cross will demonstrate a plot object being dragged, its
        % center always forced to be in the axes (as in the "crosshair" demo).
        cyancross = plot(ax, [0.7 0.7 0.6 0.6 0.8], [0.4 0.6 NaN 0.5 0.5]);
        set(cyancross, "Color", "c", "LineWidth", 2);
        draggable(cyancross, "n", [0.2 1.1 -0.1 1.1]);
        
        % The Yellow Triangle demonstrate a patch object being dragged with
        % vertical movement constraint.
        yellowtri = patch(ax, [0.4 0.6 0.6], [0.1 0.1 0.2], "y");
        draggable(yellowtri, "v", [0.1 0.7]);
        
        % The Gray Text demonstrate diagonal constraint operating on text
        % objects. We test a negative, gentle slope crossing the horizontal
        % range limits and a positive, large slope crossing the vertical range
        % limits but in which the object will be nonetheless stopped by the
        % leftmost limit.
        
        % first text
        x1 = 0.4;
        y1 = 0.9;
        m1 = -1/3;
        b1 = y1 - (m1 * x1);
        gray_1 = [0.33 0.33 0.33];
        x = ax.XLim;
        plot(ax, x, m1*x+b1, "LineStyle", "--", "Color", gray_1);
        graytext1 = text(ax, x1, y1, "This is a test");
        set(graytext1, "Color", gray_1, "EdgeColor", gray_1, ...
            "BackgroundColor", "w", "LineStyle", ":")
        draggable(graytext1, "d", m1);
        
        % second text
        x2 = 0.7;
        y2 = 0.35;
        m2 = 3;
        b2 = y2 - m2*x2;
        x = ax.XLim;
        plot(ax, x, m2*x+b2, "LineStyle", "--", "Color", gray_1);
        graytext2 = text(ax, x2, y2, "Another test");
        set(graytext2, "Color", gray_1, "EdgeColor", gray_1, ...
            "BackgroundColor", "w", "LineStyle", ":")
        draggable(graytext2, "d", m2);
    case "snapgrid"
        % The "motionfcn" argument or "on_move_callback" argument/property of
        % draggable is used here to set up a grid on which the object movement
        % is constrained.
        %
        % Furthermore, we use the function to display a "fleur" figure pointer
        % while the object moves, and the figure's WindowButtonMotionFcn to
        % display a standard "arrow" figure pointer while the mouse moves but
        % the object is not dragged.
        
        % First we set up the figure and axes.
        f = figure_fn();
        
        ax = axes_fn(f);
        axis(ax, "equal");
        box(ax, "on");
        grid(ax, "on");
        set(ax, "DataAspectRatio", [1 1 1], "XLim", [0 10], "YLim", [0 10], ...
            "XTick", 0:10, "YTick", 0:10);
        title(ax, "SNAPGRID demo for draggable.m");
        hold(ax, "on");
        
        % Now we create a cross which will snap on a grid with 1-unit spacing.
        % This is done by giving the handle to the move_cross function (see
        % below) as the "movefcn" argument of draggable.m.
        cross = plot(ax, [5 5 4.5 4.5 5.5], [4.5 5.5 NaN 5 5]);
        set(cross, "Color", "r", "LineWidth", 3);
        dr = draggable(cross, "n", [-0.5 10.5 -0.5 10.5]);
        dr.on_click_callback = @(g)cross_on_click(f, g);
        dr.on_move_callback = @cross_on_move;
        dr.on_release_callback = @(g)cross_on_release(f, g);
    case "polymove"
        % The "motionfcn" argument of draggable.m is used here to redraw the
        % polygon each time one of its vertex is moved.
        
        % Setting up the figure and axes.
        f = figure_fn();
        
        ax = axes_fn(f);
        axis(ax, "equal");
        box(ax, "on");
        set(ax, "DataAspectRatio", [1 1 1], "Xlim", [-2 2], "YLim", [-2 2]);
        title(ax, "POLYMOVE demo for draggable.m")
        hold(ax, "on");
        
        % Creating the polygon vertices
        v1 = plot(ax, -1, 0);
        v2 = plot(ax, 0, -1);
        v3 = plot(ax, 1, 0);
        v4 = plot(ax, 0, 1);
        vv = [v1 v2 v3 v4];
        % For visible vertices:
        set(vv, "Marker", "o", "MarkerSize", 10, "MarkerFaceColor", "b", ...
            "MarkerEdgeColor", "none");
        % For invisible vertices:
        % set(vv,"Marker","o","MarkerSize",10,"MarkerFaceColor","none", ...
        %     "MarkerEdgeColor","none");
        
        % Saving the vertex vector as application data in the current axes
        % (along with empty element p which will later hold the handle to the
        % polygon itself)
        setappdata(ax, "vv", vv);
        setappdata(ax, "p", []);
        
        % Calling draggable on each of the vertices, passing as an argument the
        % handle to the redraw_poly fucntion (see below)
        draggable(v1,@(g)redraw_poly(ax, g));
        draggable(v2,@(g)redraw_poly(ax, g));
        draggable(v3,@(g)redraw_poly(ax, g));
        draggable(v4,@(g)redraw_poly(ax, g));
        
        % Finally we draw the polygon itself using the redraw_poly function,
        % which can be found below
        redraw_poly(ax);
    case "sliders"
        % Creates five sliders with different draggable parameters. Slider
        % values are displayed and updated in real time. Parameters are
        % initialized uniformly and then changed using custom subsref and
        % subsasgn.
        
        % Slider numbers:
        %  1   2   3   4   5
        % Slider 1 is left mouse button.
        % Sliders 2 and 3 are middle mouse button.
        % Slider 4 is right mouse button.
        % Slider 5 is any mouse button.
        
        % Figure layout
        f = figure_fn();
        
        ax = axes_fn(f);
        axis(ax, "equal");
        box(ax, "on");
        set(ax, 'DataAspectRatio', [1 1 1], ...
            'Xlim', [0 6], 'Ylim', [-0.25 1.5], ...
            'XTick', 1:5, 'YTick', []);
        title(ax, "SLIDERS demo for draggable.m")
        hold(ax, "on");
        
        % Creating the objects
        w = 0.6; h = 0.1;
        draggables = [];
        for x = 1:5
            plot(ax, [x x], [0 1], ":k");
            value = text(ax, x, 1.2, "", "HorizontalAlignment", "center");
            slider = rectangle(ax, "Position", [x-w/2, 0.5-h/2, w, h], ...
                "Curvature", [0.05 0.05], ...
                "EdgeColor", [0 0 0.7], ...
                "FaceColor", [0 0 0.9]);
            % Storing the handle to the value text object for each slider
            setappdata(slider, "value_handle", value);
            % This will be also used as the motionfcn argument. See the function
            % definition later in the file.
            update_value(slider);
            
            % Calling draggable on sliders with a vertical constraint. We must
            % take into account the height of the sliders when setting range
            % limits. We use update_value as the "motionfcn" argument, and add
            % an "endfcn" argument which makes the slider and slider value blink
            % when the user stops dragging.
            draggables = [draggables draggable(slider, "v", [0-h/2 1+h/2], ...
                @update_value, "endfcn", @blink_value)]; %#ok<AGROW>
        end
        
        % subsasgn test
        draggables(1).FaceColor = [0.9 0.0 0.0];
        draggables(1).button = draggable.LEFT_BUTTON;
        [draggables(2:3).button] = deal(draggable.MIDDLE_BUTTON);
        draggables(4).button = draggable.RIGHT_BUTTON;
        
        % subsref test
        fprintf("draggables(1).FaceColor: [");
        fprintf("%0.1f ", draggables(1).FaceColor);
        fprintf("\b]" + newline);
        fprintf("draggables(2).FaceColor: [");
        fprintf("%0.1f ", draggables(2).FaceColor);
        fprintf("\b]" + newline);
        fprintf("draggables(1:4).button: [");
        fprintf("%s ", [draggables(1:4).button]);
        fprintf("\b\b]" + newline);
        fprintf("draggables.button: [");
        fprintf("%s, ", [draggables.button]);
        fprintf("\b\b]" + newline);
    otherwise
        disp(["Demo """ demotitle """ is not available."]);
end

end


function cross_on_click(figh, varargin)
%{
This function is passed as the "on_click_callback" property of a draggable in
the SNAPGRID demo.
%}

% We first set up the figure pointer to "fleur"
try
    figh.Pointer = "fleur";
catch
    disp( "uifigure.Pointer not supported before R2020a" );
end

end


function cross_on_move(g)
%{
This function is passed as the "on_move_callback" property of a draggable in the
SNAPGRID demo.
%}

% Then we retrieve the current cross position
cross_xdata = get(g, "XData");
cross_ydata = get(g, "YData");
cross_center = [cross_xdata(1) cross_ydata(5)];

% Computing the new position of the cross
new_position = round(cross_center);

% Updating the cross" XData and YData properties
delta = new_position - cross_center;
set(g, "XData", cross_xdata+delta(1), "YData", cross_ydata+delta(2));

end


function cross_on_release(figh, varargin)
%{
This function is passed as the "on_release_callback" property of a draggable in
the SNAPGRID demo.
%}

% Return figure pointer to "arrow"
try
    figh.Pointer = "arrow";
catch
    disp( "uifigure.Pointer not supported before R2020a" );
end

end


function redraw_poly(axh, varargin)
%{
This function is passed as the "motionfcn" argument to draggable.m in the
POLYMOVE demo. It recieves the handle to the object being dragged as its only
argument, but it is not actually used in this function.
%}

% Deleting the previous polygon
delete(getappdata(axh, "p"));

% Retrieving the vertex vector and corresponding xdata and ydata
vv = getappdata(axh, "vv");
xdata = cell2mat(get(vv, "xdata"));
ydata = cell2mat(get(vv, "ydata"));

% Plotting the new polygon and saving its handle as application data
p = plot(axh, [xdata.' xdata(1)],[ydata.' ydata(1)], "Color", "k");
setappdata(axh, "p", p);

% Putting the vertices on top of the polygon so that they are easier to drag (or
% else, the polygone line get in the way)
set(axh, "Children", [vv p]);

end


function update_value(slider)
%{
This function is passed as the "motionfcn" argument to draggable.m in the
SLIDERS demo. It recieves the handle to the slider being dragged as its only
argument and uses that to retrieve the handle to a text object which was stored
with the slider object.
%}

value_handle = getappdata(slider,"value_handle");
slider_position = get(slider,"Position");
value = slider_position(2) + slider_position(4)/2;
set(value_handle,"String",sprintf("%0.2f",value));

end


function blink_value(slider)
%{
This makes the displayed value blink; just to show usage of the "endfcn"
argument.
%}

value_handle = getappdata(slider,"value_handle");
initedgecolor = get(slider,"EdgeColor");
initfacecolor = get(slider,"FaceColor");
set(value_handle,"FontWeight","bold");
set(slider,"EdgeColor",[0 0 0.8],"FaceColor",slider.FaceColor/2);
pause(0.2)
set(value_handle,"FontWeight","normal");
set(slider,"EdgeColor",initedgecolor,"FaceColor",initfacecolor);

end
