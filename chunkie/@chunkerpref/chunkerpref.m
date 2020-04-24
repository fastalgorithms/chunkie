classdef chunkerpref    
%CHUNKERPREF class for storing preferences for chunker objects and
% sets defaults

    properties
        nchmax
        k
        dim
        nchstor
        verttol
    end
    
    methods
        function obj = chunkerpref(pref)
            %CHUNKERPREF construct chunkerpref from either a given
            % chunkerpref or struct listing properties to set

            assert(nargin <= 1,'number of arguments must be 0 or 1');
            defpref = chunkerpref.defaultprefs();
            if nargin < 1 || isempty(pref)
                pref = defpref;
            end
            assert(or(isstruct(pref),isa(pref,'chunkerpref')), ...
            'bad input, must be empty, a struct, or a chunkerpref object');
            if isa(pref,'chunkerpref')
                % do nothing if input is a chunkerpref
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
            defpref.nchmax = 10000;
            defpref.k = 16;
            defpref.dim = 2;
            defpref.nchstor = 4;
            defpref.verttol = 1e-12;
        end
    end
end