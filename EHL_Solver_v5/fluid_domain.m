classdef fluid_domain
    
    properties
        Viscosity
        De
        Wi
        beta 
        BC  = containers.Map;
    end
    
    methods
        
        function obj = fluid_domain(varargin)
            if nargin > 0
                for i = 1:2:numel(varargin)
                    propertyName = varargin{i};
                    if isprop(obj, propertyName)
                        obj.(propertyName) = varargin{i + 1};
                    else
                        error('Invalid property name: %s', propertyName);
                    end
                end
            end
            %set default BC
            obj.set_BC();
        end 
          
        function obj = set_BC(obj, varargin)
            if nargin > 1
                for i = 1:numel(varargin)
                    obj.BC(varargin) = varargin{i+1};
                end
            else % default BC. 
                obj.BC("p_i") = 0; %inlet
                obj.BC("p_o") = 0; %outlet
                obj.BC("p_c") = 0; %cavitation
            end
        end
        
        function obj = PipkinLink(reference_height)
            if isempty(De) || isempty(Wi)
                if isempty(De)
                    obj.De = obj.Wi .* reference_height
                end
                if isempty(Wi)
                    obj.Wi = obj.De ./ reference_height
                end
            end
        end
    
    end
    
end


 