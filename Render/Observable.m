%Class redefining the observer pattern.
%
%Matlab's built-in handle class already provides the addlistener and notify
%functions. The purpose of this class is to overcome a bug in the event
%notification of the handle class (in some cases the notification fails
%maybe due to interference with java components. E.g. happens when dragging
%a java jSlider too fast out of its limits).
%
%Usage : The class of the source object sending notifications should simply
%inherit from the Observable class instead of the handle class. The
%addlistener and notify functions follow the same syntax as in Matlab's
%handle class.
%
%Note: This class does not use matlab's built-in event system. Here, an
%event is simply defined by its name. The notify function directly calls
%the callbacks functions of the listeners without any safety mechanism
%(e.g. preventing infinite calls to the same callback).

classdef (Abstract) Observable < handle
    properties( SetAccess = private )
        callbacks={}
        eventStrings={}
        numListeners=0;
    end
    
    methods
        
        function addlistener(obj, eventStr, callback )
            if(isa(eventStr,'char') && isa(callback,'function_handle'))
                obj.callbacks{end+1} = callback;
                obj.eventStrings{end+1} = eventStr;
                obj.numListeners = obj.numListeners+1;
            else
                error('Wrong argument type : ''eventStr'' should be a string and ''callback'' should be function callback.')
            end
        end
        
        
        function notify(obj,eventStr)
            if(isa(eventStr,'char'))
                for i=1:obj.numListeners
                    if(strcmp(obj.eventStrings{i},eventStr))
                        event.EventName = eventStr;
                        obj.callbacks{i}(obj, event);
                    end
                end
            else
                error('Wrong argument type : ''eventStr'' should be a string.');
            end
        end
        
    end
    
end