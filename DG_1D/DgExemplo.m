% Código em DG 1D produzido por Josh Bevan, todos os créditos ao autor. 
% GitHub de Josh Bevan: https://github.com/userjjb/22_6xx-DG-for-PDEs

close all
clear all
tau=2*pi();
format long
    
N=1;
K=12;

xNode=0:1/K:1;
elemBC = reshape(sort([xNode,xNode(2:end-1)]),2,K)';

[Qx,Qw]=gauss(N+1);
L = zeros(N+1);
for m=0:N
    temp = legendre(m,Qx); 
    L(:,m+1) = temp(1,:)'; 
end

BasisWeights = zeros(K,N+1);
map = zeros(K,N+1);
for k=1:K
    map(k,:) = (elemBC(k,2)-elemBC(k,1))*Qx/2 ...
    + (elemBC(k,1)+elemBC(k,2))/2;
    for m=0:N
        BasisWeights(k,m+1) = ...
            (m+.5)*sum(Qw.*sin(tau*map(k,:)).*L(:,m+1)');
    end
end

SelfStencil = 2*(toeplitz(mod(0:N,2),double(0:N<0)))-ones(N+1);
UpwindStencil = ones(N+1); 
UpwindStencil(2:2:N+1,:)=-1;
LHSIntegral = reshape([1:2:2*N+1]'*(1./diff(xNode)),K*(N+1),1);


A=SelfStencil;
A(N+2:2*N+2,1:N+1)=UpwindStencil;

DiagStencil = zeros(N+1,3*N+2);
for i=-(2*N+1):N
    DiagStencil(1:length(diag(A,i)),i+2*N+2) = diag(A,i);
end

DiagStencil(:,2*N+2:3*N+2) = flipud(DiagStencil(:,2*N+2:3*N+2));

A = spdiags(repmat(DiagStencil,K,1),-(2*N+1):N,K*(N+1),K*(N+1));

A(1:N+1,K*(N+1)-N:K*(N+1))=UpwindStencil;

BasisWeights = reshape(BasisWeights',K*(N+1),1);

deltaT  = 0.0001;
saveT   = 0.01; 
endT    = 5;
nsaveT  = floor(saveT/deltaT);
nT      = floor(endT/deltaT);
saved   = zeros(K*(N+1),(nT/nsaveT)+1);
saved(:,1)=BasisWeights;
i=1;
for t= 1:1:nT
    if t/nsaveT==floor(t/nsaveT)    
        i= i+1;
        saved(:,i)=BasisWeights;
    end	
    k1 = LHSIntegral.*(A*BasisWeights);
    k2 = LHSIntegral.*(A* (BasisWeights + deltaT*k1 /2) );
    k3 = LHSIntegral.*(A* (BasisWeights + deltaT*k2 /2) );
    k4 = LHSIntegral.*(A* (BasisWeights + deltaT*k3) );
    BasisWeights = BasisWeights ...
        +(deltaT/6)*(k1 + 2*k2 + 2*k3 + k4);
end

Lobatto = zeros(N+3,N+1);
for m=0:N
    temp = legendre(m,[-1; Qx; 1]);
    Lobatto(:,m+1) = temp(1,:)';
end
figure('units','normalized','outerposition',[0 0 0.5 0.7]);

for i=1:size(saved,2)
    plot([elemBC(:,1) map elemBC(:,2)]', ...
        (Lobatto*reshape(saved(:,i),N+1,K)),'LineWidth',3)
    grid on
    grid minor
    axis([0 1 -1.5 1.5])
    text(1.02,0.5,'P'+string(N))
    text(1.02,.4,'K = '+string(K));
    text(1.02,-0.4,'Erro')
    text(1.02,-.5,num2str(1/sqrt(2) ...
        -rms(reshape(L*reshape(saved(:,i),N+1,K),K*(N+1),1))));
    text(1.02,1.1,'Time(s)')
    text(1.02,1,num2str((i-1)*saveT));
    pause(.001)
end
