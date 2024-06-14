classdef SphereSensor < Model
    % The measurement model for a specific origin (sensor).
    
    properties
        StateType;          % The type of trasition model
        SensorPos;          % The position of the sensor.
    end
    
    methods
        function td = SphereSensor(t, std, sensorpos, statetype)
            td.T = t;
            td.Std = std;
            td.StateType = statetype;
            td.CovMatrix = diag(std.*std);
            td.SensorPos = sensorpos;
        end
        function z_new = estimateNewState(this, state)
            Dimension = length(this.Std);
            CVlike = {'CVmodel'};
            CAlike = {'CAmodel','CTmodel','SingerModel'};
            switch Dimension
                case 2
                    switch this.StateType
                        case CVlike
                            x = state(1) - this.SensorPos(1);
                            y = state(3) - this.SensorPos(2);
                            z_new = zeros(2,1);
                            z_new(1) = sqrt(x^2 + y^2);
                            z_new(2) = atan2(y,x);
                        case CAlike
                            x = state(1) - this.SensorPos(1);
                            y = state(4) - this.SensorPos(2);
                            z_new = zeros(2,1);
                            z_new(1) = sqrt(x^2 + y^2);
                            z_new(2) = atan2(y,x);
                    end
                case 3
                    switch this.StateType
                        case CVlike
                            x = state(1) - this.SensorPos(1);
                            y = state(3) - this.SensorPos(2);
                            z = state(5) - this.SensorPos(3);
                            z_new = zeros(3,1);
                            z_new(1) = sqrt(x^2 + y^2 + z^2);
                            z_new(2) = atan2(y,x);
                            z_new(3) = atan2(z, sqrt(x^2 + y^2));
                        case CAlike
                            x = state(1) - this.SensorPos(1);
                            y = state(4) - this.SensorPos(2);
                            z = state(7) - this.SensorPos(3);
                            z_new = zeros(3,1);
                            z_new(1) = sqrt(x^2 + y^2 + z^2);
                            z_new(2) = atan2(y,x);
                            z_new(3) = atan2(z, sqrt(x^2 + y^2));
                    end
            end
        end
        function z_new = addDisturb(this, z)
            dimension = length(this.Std);
            z_new = z + diag(this.Std) * randn(dimension,1);
        end
        function z_new  =predictNewState(this, state)
            temp = this.estimateNewState(state);
            z_new = this.addDisturb(temp);
        end
        function H_linear = getLinear(this, state)
            CVlike = {'CVmodel'};
            CAlike = {'CAmodel','CTmodel','SingerModel'};
            dimension = length(this.Std);
            switch dimension
                case 2
                    switch this.StateType
                        case CVlike
                            H_linear = zeros(2,4);
                            x = state(1) - this.SensorPos(1);
                            if x == 0
                                x = 1e-10;
                            end
                            y = state(3) - this.SensorPos(2);
                            H_linear(1,1) = x / sqrt(x^2 + y^2);
                            H_linear(1,3) = y / sqrt(x^2 + y^2);
                            H_linear(2,1) = 1 / (1 + (y/x)^2) * (-y/(x^2));
                            H_linear(2,3) = 1 / (1 + (y/x)^2) * (1 / x);
                        case CAlike
                            H_linear = zeros(2,6);
                            x = state(1) - this.SensorPos(1);
                            y = state(4) - this.SensorPos(2);
                            if x == 0
                                x = 1e-10;
                            end
                            H_linear(1,1) = x / sqrt(x^2 + y^2);
                            H_linear(1,4) = y / sqrt(x^2 + y^2);
                            H_linear(2,1) = 1 / (1 + (y/x)^2) * (-y/ (x^2));
                            H_linear(2,4) = 1 / (1 + (y/x)^2) * (1/x);
                    end
                case 3
                    switch this.StateType
                        case CVlike
                            H_linear = zeros(3,6);
                            x = state(1) - this.SensorPos(1);
                            y = state(3) - this.SensorPos(2);
                            z = state(5) - this.SensorPos(3);
                            if x == 0
                                x = 1e-10;
                            end
                            H_linear(1,1) = x / sqrt(x^2 + y^2 + z^2);
                            H_linear(1,3) = y / sqrt(x^2 + y^2 + z^2);
                            H_linear(1,5) = z / sqrt(x^2 + y^2 + z^2);
                            H_linear(2,1) = 1 / (1 + (y/x)^2) * (-y / (x^2));
                            H_linear(2,3) = 1 / (1 + (y/x)^2) * (1/x);
                            r = sqrt(x^2 + y^2);
                            H_linear(3,1) = 1 / (1 + (z/r)^2) * (-z * x / (r^3));
                            H_linear(3,3) = 1 / (1 + (z/r)^2) * (-z * x / (r^3));
                            H_linear(3,5) = 1 / (1 + (z/r)^2) * (1/r);
                        case CAlike
                            H_linear = zeros(3,9);
                            x = state(1) - this.SensorPos(1);
                            y = state(4) - this.SensorPos(2);
                            z = state(7) - this.SensorPos(3);
                            if x == 0
                                x = 1e-10;
                            end
                            r = sqrt(x^2 + y^2);
                            H_linear(1,1) = x / sqrt(x^2 + y^2 + z^2);
                            H_linear(1,4) = y / sqrt(x^2 + y^2 + z^2);
                            H_linear(1,7) = z / sqrt(x^2 + y^2 + z^2);
                            H_linear(2,1) = 1 / (1 + (y/x)^2) * (-y / (x^2));
                            H_linear(2,4) = 1 / (1 + (y/x)^2) * (1/x);
                            H_linear(3,1) = 1 / (1 + (z/r)^2) * (-z * x / (r^3));
                            H_linear(3,4) = 1 / (1 + (z/r)^2) * (-z * x / (r^3));
                            H_linear(3,7) = 1 / (1 + (z/r)^2) * (1/r);
                    end
            end
        end
        function prob = getProb(this, measure, state)
            z = this.estimateNewState(state);
            e = z - measure;
            dimension = length(this.Std);
            prob = 1 / sqrt((2*pi)^dimension * det(this.CovMatrix)) * exp(-1/2 * e' / this.CovMatrix * e);
        end
        function pos = getCartesianPos(this, z)
            r = z(1);
            theta = z(2);
            pos = [r * cos(theta) + this.SensorPos(1); r*sin(theta) + this.SensorPos(2)];
        end
    end
    
end

