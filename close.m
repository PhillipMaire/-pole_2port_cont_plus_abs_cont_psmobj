% [] = close(obj)    Method that gets called when the protocol is closed
%
% This method deletes any and all SoloParamHandles that are owned by
% obj. It assumes that it has been given read-only access to a
% SoloParamHandle named myfig which holds a handle to the main protocol
% figure; this method deletes that figure.
%
% Add any other necessary figure closures or deletions to this code.
%

% CDB Feb 06

function [] = close(obj)

GetSoloFunctionArgs;

if ~isempty(whos('global','pole_motors'))
    global pole_motors
    if isa(pole_motors,'zaberTCD1000')
        close_and_cleanup(pole_motors)
    end
    clear global pole_motors
end

if ~isempty(whos('global','lickport_motors'))
    global lickport_motors
    if isa(lickport_motors,'zaberTCD1000')
        close_and_cleanup(lickport_motors)
    end
    clear global lickport_motors
end

if ~isempty(whos('global','state_machine_properties'))
    global state_machine_properties;
    Close(state_machine_properties.sm); % Close state machine on RTLinux box, to keep number of threads under control.
end

delete(value(myfig));

delete_sphandle('owner', ['^@' class(obj) '$']);



