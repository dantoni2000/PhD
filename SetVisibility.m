function SetVisibility(hObject,eventdata)
switch eventdata.Key
    case {'uparrow','up'}
        direction=1;
    case {'downarrow','down'}
        direction=-1;
    otherwise
        return
end
handles=guidata(hObject);
%set the old plot to invisible
set(handles.h(handles.visible_index),'Visible','off');
%find the new index to be set
handles.visible_index=mod(handles.visible_index+direction,length(handles.h));
if handles.visible_index==0,handles.visible_index=length(handles.h);end
%set the new plot to visible
set(handles.h(handles.visible_index),'Visible','on');
%optional: set the title to the index value
title(sprintf('selected plot: %d',handles.visible_index))
guidata(hObject,handles)
end