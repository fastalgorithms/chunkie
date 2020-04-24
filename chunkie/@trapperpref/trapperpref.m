classdef trapperpref
%TRAPPERPREF class for storing preferences for trapper objects and
% sets defaults

    properties
        dim
        nptstor
        nptmax
    end
    
    methods
        function obj = trapperpref(pref)
            %TRAPPERPREF construct trapperpref from either a given
            % trapperpref or struct listing properties to set

            assert(nargin <= 1,'number of arguments must be 0 or 1');
            defpref = trapperpref.defaultprefs();
            if nargin < 1 || isempty(pref)
                pref = defpref;
            end
            assert(or(isempty(pref),or(isstruct(pref),isa(pref,'trapperpref'))), ...
            'bad input, must be empty, a struct, or a trapperpref object');
            if isa(pref,'trapperpref')
                % do nothing if input is a trapperpref
                obj = pref;
            else
                deffields = fieldnames(defpref);
                fields = fieldnames(pref);
                
                % assign defaults first, then overwrite
                for i = 1:length(deffields)
                    obj.(deffields{i}) = defpref.(deffields{i});
                end
                for i = 1:length(fields)
                    obj.(fields{i}) = pref.(fields{i});
                end
            end
        end
    end
    methods (Static)
        function defpref = defaultprefs()
            defpref = [];
            defpref.nptstor = 10000;
            defpref.nptmax = 1000000;
            defpref.dim = 2;
        end
    end
end
