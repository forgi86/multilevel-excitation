%%   MODEL STATES %%

% x(1) = m0
% x(2) = m1
% x(3) = m2
% x(4) = m3
% x(5) = m4
 
 
%% CODE %%
function continuous_reaction_model(block)
    setup(block);
end 
 
function setup(block)
 
    %% Register number of input and output ports
    block.NumInputPorts  = 1;
    block.NumOutputPorts = 3;
 
    %% Setup functional port properties to dynamically
    %% inherited.
    block.SetPreCompInpPortInfoToDynamic;
    block.SetPreCompOutPortInfoToDynamic;
 
    block.InputPort(1).Complexity   = 'Real'; 
    block.InputPort(1).DataTypeId   = 0;
    block.InputPort(1).SamplingMode = 'Sample';
    block.InputPort(1).Dimensions   = 1;
    
    block.OutputPort(1).Complexity   = 'Real'; 
    block.OutputPort(1).DataTypeId   = 0;
    block.OutputPort(1).SamplingMode = 'Sample';
    block.OutputPort(1).Dimensions   = 1;
    
    block.OutputPort(2).Complexity   = 'Real'; 
    block.OutputPort(2).DataTypeId   = 0;
    block.OutputPort(2).SamplingMode = 'Sample';
    block.OutputPort(2).Dimensions   = 1;

    
    block.OutputPort(3).Complexity   = 'Real'; 
    block.OutputPort(3).DataTypeId   = 0;
    block.OutputPort(3).SamplingMode = 'Sample';
    block.OutputPort(3).Dimensions   = 1;
    

    % Register parameters
    block.NumDialogPrms = 2;
 
    %% Set block sample time to continous
    block.SampleTimes = [0 0];
 
    %% Setup Dwork
    block.NumContStates = 3;
 
    %% Register Methods
    block.RegBlockMethod('InitializeConditions',    @InitConditions);  
    block.RegBlockMethod('Outputs',                 @Output);  
    block.RegBlockMethod('Derivatives',             @Derivatives);   
    block.RegBlockMethod('PostPropagationSetup', @DoPostPropSetup);
    block.RegBlockMethod('SetInputPortSamplingMode', @DoInputPortSamplingMode);
end
 
function InitConditions(block)
    x0 = block.DialogPrm(2).Data;
    block.ContStates.Data = x0;  
end

function DoPostPropSetup(block)
end

function DoInputPortSamplingMode(block,idx,fd)
    block.InputPort(idx).SamplingMode = fd;
    block.InputPort(idx).SamplingMode = fd;
    
    block.OutputPort(1).SamplingMode = fd;
    block.OutputPort(2).SamplingMode = fd;
    block.OutputPort(3).SamplingMode = fd;
end
 
function Derivatives(block)
    x = block.ContStates.Data;
    
    der = zeros(size(x));
    
    coeff = block.DialogPrm(1).Data;
    F = coeff.Fss;
    Ca0 = coeff.Ca0;
    E   = coeff.E;
    VR = coeff.VR;
    T0 = coeff.T0;
    TCin = coeff.TCin;
    lambda = coeff.lambda;
    k0 = coeff.k0;
    UAJ = coeff.UAJ;
    R = coeff.R;
    rho = coeff.rho;
    cp = coeff.cp;
    VJ = coeff.VJ;
    rhoJ = coeff.rhoJ;
    cJ = coeff.cJ;
    
    
    Ca = x(1);
    TR = x(2);
    TJ = x(3);
    
    FJ = block.InputPort(1).Data;
    
    k = k0*exp(-E/(R*TR));
    der(1) = F/VR*(Ca0 - Ca) - Ca*k;
    der(2) = F/VR*(T0 - TR) - lambda*Ca*k/(rho*cp) - UAJ*(TR - TJ)/(VR*rho*cp);
    der(3) = FJ/VJ*(TCin - TJ) + UAJ*(TR-TJ)/(VJ*rhoJ*cJ);
    block.Derivatives.Data = (der);
end
 
function Output(block)
    x = block.ContStates.Data; % state
    Ca = x(1);
    TR = x(2);
    TJ = x(3);
    block.OutputPort(1).Data = Ca;
    block.OutputPort(2).Data = TR;
    block.OutputPort(3).Data = TJ;    
end
