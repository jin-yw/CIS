classdef CartesianSensor < Model
    % The measurement model for a specific origin (sensor).
    
    properties
        H;                  % The measurement matrix
        StateType;          % The type of transition model.
        SensorPos;          % The position of the sensor.
    end
    
    methods
        %% Constructor
        function td = CartesianSensor(t, std, sensorPos, stateType)
            CVlike = {'CVmodel'};
            CAlike = {'CAmodel'};
            td.T = t;
            td.Std = std;
            td.StateType = stateType;
            td.SensorPos = sensorPos;
            dimension = length(std);
            switch stateType
                case CVlike
                    h = [1,0];
                case CAlike
                    h = [1,0,0];
            end
            td.H = [];
            for j = 1:dimension
                td.H = blkdiag(td.H, h);
            end
            % Compute the covariance matrix
            td.CovMatrix = diag(std.*std);
        end
        %% Measurement function
        % Obtain the measurement without noise
        function z_new = estimateNewState(this, state)
            z_new = this.H * state - this.SensorPos;
        end
        % Add noise to the measurement
        function z_new = addDisturb(this, z)
            dimension = length(this.Std);
            z_new = z + chol(this.CovMatrix) * randn(dimension,1);
        end
        % Obtain the measurement with noise
        function z_new = predictNewState(this,state)
            z = this.estimateNewState(state);
            z_new = this.addDisturb(z);
        end
        %% Obtain the linear approximation for the measurement model
        function H_linear = getLinear(this, ~)
            H_linear = this.H;
        end
        %% The likelihood function
        function prob = getProb(this, measure, state)
            z = this.estimateNewState(state);
            e = z - measure;
            prob = 1 / sqrt(det(2*pi*this.CovMatrix)) * exp(-1/2 * e' / this.CovMatrix * e);
        end
    end
    
end

