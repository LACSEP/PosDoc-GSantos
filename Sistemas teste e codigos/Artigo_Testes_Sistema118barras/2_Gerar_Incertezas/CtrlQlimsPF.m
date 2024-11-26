function [diffy_Qg,diffy_v,Y] = CtrlQlimsPF(sigk,sigv,sigq,Vm,Vr,Qg,Qmax,Qmin,gen,nger)

%Reactive power limits
qlimsup = Qmax - sigq;
qliminf = Qmin + sigq;
%Voltage limits
vlimsup = Vr + sigv;
vliminf = Vr - sigv;

%%
%Control equations
k = 0;
for i = 1 : 1 : nger
    j = gen(i);
    k = k + 1;
    qlimsch_ch1 = 1 / (1 + exp(-sigk*(Qg(i)-qlimsup(j)))); %Reactive power upper switch
    qlimsch_ch2 = 1 / (1 + exp(sigk*(Qg(i)-qliminf(j)))); %Reactive power lower switch
    qlimsch_ch3 = 1 / (1 + exp(sigk*(Vm(j)-vlimsup(j)))); %Voltage upper switch
    qlimsch_ch4 = 1 / (1 + exp(-sigk*(Vm(j)-vliminf(j)))); %Lower voltage switch
    
    Ynormal(k,1) = (1 - qlimsch_ch1 * qlimsch_ch3) * (1 - qlimsch_ch2 * qlimsch_ch4) * (Vm(j) - Vr(j));
    Ysuperior(k,1) =  (qlimsch_ch1 * qlimsch_ch3) * (1 - qlimsch_ch2 * qlimsch_ch4) * (Qg(i) - Qmax(j));
    Yinferior(k,1) = (1 - qlimsch_ch1 * qlimsch_ch3) * (qlimsch_ch2 * qlimsch_ch4) * (Qg(i) - Qmin(j));
    
    x1 = exp(min(100,-sigk*(Qg(i) - Qmax(j) + sigq)));
    x2 = exp(min(100,sigk*(-sigv + Vm(j) - Vr(j))));
    x3 = exp(min(100,sigk*(Qg(i) - Qmin(j) - sigv)));
    x4 = exp(min(100,-sigk*(sigv + Vm(j) - Vr(j))));
    x5 = exp(min(100,2*sigk*(Qg(i) - Qmin(j) - sigv)));
    x6 = exp(min(100,-2*sigk*(Qg(i) - Qmax(j) + sigq)));
    x7 = exp(min(100,-2*sigk*(sigv + Vm(j) - Vr(j))));
    x8 = exp(min(100,2*sigk*(-sigv + Vm(j) - Vr(j))));
    
    
    diffYnormal_Qg(k,1) = sigk*(1 - 1/((1 + x1)*(x2 + 1)))*(Vm(j) - Vr(j))*x3/((1 + x4)*(x3 + 1)^2) - sigk*(1 - 1/((1 + x4)*(x3 + 1)))*(Vm(j) - Vr(j))*x1/((1 + x1)^2*(x2 + 1));
    diffYsuperior_Qg(k,1) = sigk*(1 - 1/((1 + x4)*(x3 + 1)))*(Qg(i) - Qmax(j))*x1/((1 + x1)^2*(x2 + 1)) + sigk*(Qg(i) - Qmax(j))*x3/((1 + x1)*(1 + x4)*(x3 + 1)^2*(x2 + 1)) + (1 - 1/((1 + x4)*(x3 + 1)))/((1 + x1)*(x2 + 1));
    diffYinferior_Qg(k,1) = -sigk*(1 - 1/((1 + x1)*(x2 + 1)))*(Qg(i) - Qmin(j))*x3/((1 + x4)*(x3 + 1)^2) - sigk*(Qg(i) - Qmin(j))*x1/((1 + x1)^2*(1 + x4)*(x3 + 1)*(x2 + 1)) + (1 - 1/((1 + x1)*(x2 + 1)))/((1 + x4)*(x3 + 1));
    
    diffYnormal_v(k,1) = -sigk*(1 - 1/((1 + x1)*(x2 + 1)))*(Vm(j) - Vr(j))*x4/((1 + x4)^2*(x3 + 1)) + sigk*(1 - 1/((1 + x4)*(x3 + 1)))*(Vm(j) - Vr(j))*x2/((1 + x1)*(x2 + 1)^2) + (1 - 1/((1 + x1)*(x2 + 1)))*(1 - 1/((1 + x4)*(x3 + 1)));
    diffYsuperior_v(k,1) = -sigk*(1 - 1/((1 + x4)*(x3 + 1)))*(Qg(i) - Qmax(j))*x2/((1 + x1)*(x2 + 1)^2) - sigk*(Qg(i) - Qmax(j))*x4/((1 + x1)*(1 + x4)^2*(x3 + 1)*(x2 + 1));
    diffYinferior_v(k,1) = sigk*(1 - 1/((1 + x1)*(x2 + 1)))*(Qg(i) - Qmin(j))*x4/((1 + x4)^2*(x3 + 1)) + sigk*(Qg(i) - Qmin(j))*x2/((1 + x1)*(1 + x4)*(x3 + 1)*(x2 + 1)^2);   
end
Y = [Ynormal + Ysuperior + Yinferior];

diffy_Qg = [diffYnormal_Qg + diffYsuperior_Qg + diffYinferior_Qg];
%diffy_Qg(diffy_Qg>0.99) = 1;
%diffy_Qg(diffy_Qg<0.01) = 0;

diffy_v = [diffYnormal_v + diffYsuperior_v + diffYinferior_v];
%diffy_v(diffy_v>0.99) = 1;
%diffy_v(diffy_v<0.01) = 0;



