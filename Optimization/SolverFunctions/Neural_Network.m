classdef Neural_Network
    properties
        %define structure parameters
        inputLayerSize
        outputLayerSize
        %define wieght parameters
        Wlayer1
        %define bias
        blayer1
        %define scale factor for how complex you want to allow your system
        %to be
        lambda
        %is it a classification network (generators off/on)
        classify
        %constant for node function
        nodeconst
        %statistical stuff to update with more info
        avrginputs
        stddev
    end
    methods
        function obj = Neural_Network(inputLayerSize, outputLayerSize, varargin)
            %initialization
            %inputs include UpperBound, LowerBound, Demand, quadratic
            %portion of cost, linear portion of cost, cost/kWh from grid
            %for each generator, so inputlaysize is 4*#ofgenerators+2
            if isnumeric(inputLayerSize)
                obj.inputLayerSize = inputLayerSize;
                obj.Wlayer1 = rand(obj.inputLayerSize,outputLayerSize);%give a different weight to each input's connection to each node
            end
            if isnumeric(outputLayerSize)
                obj.outputLayerSize = outputLayerSize;
                obj.blayer1 = rand(1,outputLayerSize);
            end     
            %default conditions
            obj.lambda = .0001;
            obj.classify = false;
            obj.nodeconst = 1;
            obj.avrginputs = [];
            obj.stddev = [];
            
            if length(varargin)==1%the first input is if it is a classification network, then node const, then lambda
               classify = iscellstr(varargin);
               lambda = isnumeric(varargin);
               obj.lambda = varargin(lambda);
               obj.classify = (nnz(classify)>0);
            elseif length(varargin)==2
                obj.classify = strcmp(varargin{1},'classify');
                obj.nodeconst = varargin{2};
            end
        end
        function yHat = forward(self, X)
            %forward propogate inputs, X is the inputs, this must be in a
            %genparameters x number of outputs size
            if sum(size(X')==size(self.Wlayer1))==2%if inputs directly allow multiplication
                z2 = X.*self.Wlayer1';
                yHat = sum(z2,2)+self.blayer1';
            else%if only one row of inputs, or one row of inputs per timestep
                yHat = (X*self.Wlayer1 + self.blayer1);
            end
            yHat = activationf(self, yHat); %use the activation function scaled by 1
        end
        function a2 = activationf(self,yHat)
            %apply activation function
            %if it is a classifyer use a sigmoid function
            if self.classify
                a2 = 1./(1+exp(-self.nodeconst*yHat));
                if nnz(isnan(a2))>0
                    a2(and(isnan(a2),yHat>0)) = 1;%if numbers are too big inf/inf, make a2=1
                    a2(isnan(a2)) = 0;%if numbers are to negative -inf/-inf, make a2=0
                end
            else %if it is a numeric output network don't include activation
                a2 = yHat*self.nodeconst;
            end
        end
    end
end
