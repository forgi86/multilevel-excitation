%%   MODEL STATES %%

% x(1) = m0
% x(2) = m1
% x(3) = m2
% x(4) = m3
% x(5) = m4
 
 
%% CODE %%
function continuous_reaction_model_scaled_sensitivity(block)
    setup(block);
end 
 
function setup(block)
 
    %% Register number of input and output ports
    block.NumInputPorts  = 1;
    block.NumOutputPorts = 3 + 4;
 
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
    
    block.OutputPort(4).Complexity   = 'Real'; 
    block.OutputPort(4).DataTypeId   = 0;
    block.OutputPort(4).SamplingMode = 'Sample';
    block.OutputPort(4).Dimensions   = 3;

    block.OutputPort(5).Complexity   = 'Real'; 
    block.OutputPort(5).DataTypeId   = 0;
    block.OutputPort(5).SamplingMode = 'Sample';
    block.OutputPort(5).Dimensions   = 3;

    block.OutputPort(6).Complexity   = 'Real'; 
    block.OutputPort(6).DataTypeId   = 0;
    block.OutputPort(6).SamplingMode = 'Sample';
    block.OutputPort(6).Dimensions   = 3;
    
    block.OutputPort(7).Complexity   = 'Real'; 
    block.OutputPort(7).DataTypeId   = 0;
    block.OutputPort(7).SamplingMode = 'Sample';
    block.OutputPort(7).Dimensions   = 3;

    % Register parameters
    block.NumDialogPrms = 2;
 
    %% Set block sample time to continous
    block.SampleTimes = [0 0];
 
    %% Setup Dwork
    block.NumContStates = 3 + 4*3; % 3 states + (4 parameters)*(3 states)
 
    %% Register Methods
    block.RegBlockMethod('InitializeConditions',    @InitConditions);  
    block.RegBlockMethod('Outputs',                 @Output);  
    block.RegBlockMethod('Derivatives',             @Derivatives);   
    block.RegBlockMethod('PostPropagationSetup', @DoPostPropSetup);
    block.RegBlockMethod('SetInputPortSamplingMode', @DoInputPortSamplingMode);
end
 
function InitConditions(block)
    x0 = block.DialogPrm(2).Data;
    block.ContStates.Data(1:3) = x0;  
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
    KTR = UAJ/(VR*rho*cp);
    KTJ = UAJ/(VJ*rhoJ*cJ);
    
    
    Ca = x(1);
    TR = x(2);
    TJ = x(3);
    
    sk0 = x(4:6);
    sE  = x(7:9);
    slambda = x(10:12);
    sUAJ = x(13:15);
    
    FJ = block.InputPort(1).Data;
    FJ = FJ/1000; % scaling
    
    k = k0*exp(-E/(R*TR));
    der(1) = F/VR*(Ca0 - Ca) - Ca*k;
    der(2) = F/VR*(T0 - TR) - lambda*Ca*k/(rho*cp) - KTR*(TR - TJ);
    der(3) = FJ/VJ*(TCin - TJ) + KTJ*(TR - TJ); 

    % sk0
    der(4) = -(F*sk0(1))/VR - (exp(-E/(R*TR))*(Ca*R*TR^2 + R*TR^2*k0*sk0(1) + Ca*E*k0*sk0(2)))/(R*TR^2);
    der(5) =  -(Ca*R*TR^2*lambda*exp(-E/(R*TR)) + R*TR^2*k0*lambda*sk0(1)*exp(-E/(R*TR)) + Ca*E*k0*lambda*sk0(2)*exp(-E/(R*TR)))/(R*TR^2*cp*rho) - (R*TR^2*UAJ*sk0(2) - R*TR^2*UAJ*sk0(3) + F*R*TR^2*cp*rho*sk0(2))/(R*TR^2*VR*cp*rho);
    der(6) = (UAJ*(sk0(2) - sk0(3)))/(VJ*cJ*rhoJ) - (FJ*sk0(3))/VJ;
 
    % kE
    der(7) = - (F*sE(1))/VR - (exp(-E/(R*TR))*(R*k0*sE(1)*TR^2 - Ca*k0*TR + Ca*E*k0*sE(2)))/(R*TR^2);
    der(8) = (Ca*TR*VR*k0*lambda*exp(-E/(R*TR)) - Ca*E*VR*k0*lambda*sE(2)*exp(-E/(R*TR)))/(R*TR^2*VR*cp*rho) - (R*UAJ*sE(2) - R*UAJ*sE(3) + F*R*cp*rho*sE(2) + R*VR*k0*lambda*sE(1)*exp(-E/(R*TR)))/(R*VR*cp*rho);
    der(9) = (UAJ*(sE(2) - sE(3)))/(VJ*cJ*rhoJ) - (FJ*sE(3))/VJ;
 
    % slambda
   der(10) = - (F*slambda(1))/VR - (exp(-E/(R*TR))*(R*k0*slambda(1)*TR^2 + Ca*E*k0*slambda(2)))/(R*TR^2);
   der(11) = - (Ca*R*TR^2*k0*exp(-E/(R*TR)) + R*TR^2*k0*lambda*slambda(1)*exp(-E/(R*TR)) + Ca*E*k0*lambda*slambda(2)*exp(-E/(R*TR)))/(R*TR^2*cp*rho) - (R*TR^2*UAJ*slambda(2) - R*TR^2*UAJ*slambda(3) + F*R*TR^2*cp*rho*slambda(2))/(R*TR^2*VR*cp*rho);
   der(12) =  (UAJ*(slambda(2) - slambda(3)))/(VJ*cJ*rhoJ) - (FJ*slambda(3))/VJ;
   
   % sUAJ
    der(13) = - (F*sUAJ(1))/VR - (exp(-E/(R*TR))*(R*k0*sUAJ(1)*TR^2 + Ca*E*k0*sUAJ(2)))/(R*TR^2);
    der(14) = - (R*TR^2*k0*lambda*sUAJ(1)*exp(-E/(R*TR)) + Ca*E*k0*lambda*sUAJ(2)*exp(-E/(R*TR)))/(R*TR^2*cp*rho) - (R*TR^3 - R*TJ*TR^2 + R*TR^2*UAJ*sUAJ(2) - R*TR^2*UAJ*sUAJ(3) + F*R*TR^2*cp*rho*sUAJ(2))/(R*TR^2*VR*cp*rho);
    der(15) = - (FJ*sUAJ(3))/VJ - (TJ - TR - UAJ*sUAJ(2) + UAJ*sUAJ(3))/(VJ*cJ*rhoJ);
 
    
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
    
    block.OutputPort(4).Data = x(4:6); % sk0
    block.OutputPort(5).Data = x(7:9); % sE
    block.OutputPort(6).Data = x(10:12); % slambda
    block.OutputPort(7).Data = x(13:15); % sUAJ

end
