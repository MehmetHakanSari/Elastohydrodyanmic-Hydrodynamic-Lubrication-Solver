classdef simulation_proporties
    properties
        relaxp = 1.95; % Sovrarelaxation rate on Reynolds pressure
        relaxD = 0.1; % Underrelaxation rate on EHL dispalacement (initial value)
        relaxL = 0.5; % Underrelaxation rate on load control
        
        maxiterR = 1e8;

        relaxDmax = 1e-2; % maximum underrelaxation rate on EHL dispalacement for Aitk. acc
        relaxDmin = 1e-3; % minimum underrelaxation rate on EHL dispalacement for Aitk. acc

        tollp = 1e-9; % Error on on Reynolds pressure
        tollD = 1e-6; % Error on EHL dispalacement
        tollL = 1e-3; % Error on EHL load control
    end

    methods
        function obj = simulation_proporties(varargin)
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
        end
    end
end