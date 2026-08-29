function options = EELS_ZEST_ExportDialog(title_text, scopes, formats)
%EELS_ZEST_EXPORTDIALOG Choose export scope, format, and workspace copy.
if nargin < 1 || isempty(title_text), title_text = 'Export'; end
if nargin < 2 || isempty(scopes), scopes = {'All'}; end
if nargin < 3 || isempty(formats), formats = {'MAT'}; end
options = struct('cancelled', true, 'scope_index', 1, 'format', formats{1}, ...
    'assign_workspace', false);
fig = dialog('Name', title_text, 'Position', [420 310 430 245], ...
    'Color', 'w', 'WindowStyle', 'modal', 'CloseRequestFcn', @cancelDialog);
uicontrol(fig, 'Style', 'text', 'String', 'Range', 'Position', [35 180 80 25], ...
    'BackgroundColor', 'w', 'HorizontalAlignment', 'left');
scope_control = uicontrol(fig, 'Style', 'popupmenu', 'String', scopes, ...
    'Position', [120 180 270 28]);
uicontrol(fig, 'Style', 'text', 'String', 'Format', 'Position', [35 130 80 25], ...
    'BackgroundColor', 'w', 'HorizontalAlignment', 'left');
format_control = uicontrol(fig, 'Style', 'popupmenu', 'String', formats, ...
    'Position', [120 130 270 28]);
workspace_control = uicontrol(fig, 'Style', 'checkbox', ...
    'String', 'Also assign to base workspace', 'Value', 0, ...
    'Position', [120 85 240 25], 'BackgroundColor', 'w');
uicontrol(fig, 'Style', 'pushbutton', 'String', 'Export', 'Position', [210 28 85 32], ...
    'FontWeight', 'bold', 'Callback', @acceptDialog);
uicontrol(fig, 'Style', 'pushbutton', 'String', 'Cancel', 'Position', [305 28 85 32], ...
    'Callback', @cancelDialog);
uiwait(fig);
if isgraphics(fig), delete(fig); end

    function acceptDialog(~, ~)
        options.cancelled = false;
        options.scope_index = scope_control.Value;
        options.format = formats{format_control.Value};
        options.assign_workspace = logical(workspace_control.Value);
        uiresume(fig);
    end
    function cancelDialog(~, ~)
        options.cancelled = true;
        uiresume(fig);
    end
end
