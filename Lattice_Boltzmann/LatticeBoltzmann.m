%  c4  c3   c2  Modelo D2Q9
%    \  |  /    
%  c5 -c9 - c1  
%    /  |  \    
%  c6  c7   c8 

clear all 
clc 

%salvar .mat para Matlab/Octave
savematt = false;
freqmatt = 10;    
    
%salvar VTK para o Paraview
vtktf = false; 
freqvtk = 10;

%Ligar ou desligar a plotagem
ploton = true;
freqplot = 10;

tmax=80000;

if (exist([pwd '\save'],'dir')==0) %criar um diretório para salvar os dados
    mkdir([pwd '\save']) 
end 


%Pesos
w9 = 4/9; 
w1 = 1/9; w3 = 1/9; w5 = 1/9; w7 = 1/9; 
w2 = 1/36; w4 = 1/36; w6 = 1/36; w8 = 1/36; 

% Ler a imagem e transforma-la no dominio da simulaçao
imdata = imread('basepronta14.bmp');
asa = double(imdata);
asa = abs(asa-1);
[ny, nx] = size(asa);

J = asa;
J(1,:) = 1;
J(ny,:) = 1;

BOUND = J'; 

alpha = 0.03; 
u0 = 0.0974; 
densidade = 1.0; 
Re = densidade*u0*ny/alpha; 
omega = 1.0/(3*alpha+0.5);

F = repmat(densidade/9,[nx ny 9]); 
FEQ = F; 
CI = [0:nx*ny:nx*ny*7]; 

ON=find(BOUND); 

refletir=[ON+CI(1) ON+CI(2) ON+CI(3) ON+CI(4) ...
            ON+CI(5) ON+CI(6) ON+CI(7) ON+CI(8)];
        
refletido= [ON+CI(5) ON+CI(6) ON+CI(7) ON+CI(8) ...
            ON+CI(1) ON+CI(2) ON+CI(3) ON+CI(4)];
 
        
% Condições inicias
F(:,:,1) = u0; 
F(:,:,2) = .0238;
F(:,:,3) = .0739;
F(:,:,4) = .0143;
F(:,:,5) = .0574;
F(:,:,6) = .0144;
F(:,:,7) = .0739;
F(:,:,8) = .0238;
F(:,:,9) = .2955;

if (ploton)
    figure('units','normalized','outerposition',[0 0 1 1]);
    hold on; 
    axis([1 +nx+1 1 ny+1]);
    axis ij
end

t=0;

while t<tmax
    %Streaming
    F(:,:,4)=F([2:nx nx],[ny 1:ny-1],4); 
    F(:,:,3)=F(:,[ny 1:ny-1],3);
    F(:,:,2)=F([1 1:nx-1],[ny 1:ny-1],2);
    F(:,:,5)=F([2:nx nx],:,5);  
    F(:,:,1)=F([1 1:nx-1],:,1);  
    F(:,:,6)=F([2:nx nx],[2:ny 1],6); 
    F(:,:,7)=F(:,[2:ny 1],7); 
    F(:,:,8)=F([1 1:nx-1],[2:ny 1],8);
    
%     if (mod(t,100)==0) 
%         F(:,1,1) = u0; %Forcar pra direita
%     end
    
       
    BOUNCEDBACK=F(refletir);     
    
    rho=sum(F,3); 
    ux=(1./rho).*(sum(F(:,:,[1 2 8]),3)-sum(F(:,:,[4 5 6]),3)); 
    uy=(1./rho).*(sum(F(:,:,[2 3 4]),3)-sum(F(:,:,[6 7 8]),3));   
    
    ux(ON)=0; uy(ON)=0; rho(ON)=0;     
    
    uxx =  ux.^2 ;
    uyy =  uy.^2 ;
    uxy = 2 .* ux .* uy ;
    uu = uxx + uyy ; 
    
    % Direções LI
    FEQ(:,:,1) = w1*rho.*(1 + 3*ux + 4.5*uxx - 1.5*uu);
    FEQ(:,:,3) = w3*rho.*(1 + 3*uy + 4.5*uyy - 1.5*uu);
    FEQ(:,:,5) = w5*rho.*(1 - 3*ux + 4.5*uxx - 1.5*uu);
    FEQ(:,:,7) = w7*rho.*(1 - 3*uy + 4.5*uyy - 1.5*uu);
    % Direções LD
    FEQ(:,:,2) = w2*rho.*(1 + 3*(ux+uy) + 4.5*(uu+uxy) - 1.5*uu);
    FEQ(:,:,4) = w4*rho.*(1 - 3*(ux-uy) + 4.5*(uu-uxy) - 1.5*uu);
    FEQ(:,:,6) = w6*rho.*(1 - 3*(ux+uy) + 4.5*(uu+uxy) - 1.5*uu);
    FEQ(:,:,8) = w8*rho.*(1 + 3*(ux-uy) + 4.5*(uu-uxy) - 1.5*uu);
    % Equilibrio
    FEQ(:,:,9) = rho - sum(FEQ(:,:,[1:8]),3); 
    
    F = (omega*FEQ+(1-omega)*F);  %colisão
    
    F(refletido) = BOUNCEDBACK;
    result = (vort(ux,uy))';
    
    t=t+1;  
    
    if (mod(t,freqvtk)==0 && vtktf)
        speed=sqrt(ux(1:nx,:).^2'+uy(1:nx,:).^2');
        vtkwrite('save/speed'+string(t)+'.vtk', 'structured_points', 'speed',speed);
        vtkwrite('save/ux'+string(t)+'.vtk', 'structured_points', 'ux',ux);
        vtkwrite('save/uy'+string(t)+'.vtk', 'structured_points', 'uy',uy);
        vtkwrite('save/vort'+string(t)+'.vtk', 'structured_points', 'vort',result);
        fprintf("Ciclo: %d.\n", t)
    end
    
    if (mod(t,freqmatt)==0 && savematt)
        speed=sqrt(ux(1:nx,:).^2'+uy(1:nx,:).^2');
        save('save/record'+string(t)+'.mat', 'speed','ux', 'uy', 'result','nx', 'ny', 'rho');
        fprintf("Ciclo: %d.\n", t)
    end
    
    if (mod(t,freqplot)==0 && ploton) 
        speed=sqrt(ux(1:nx,:).^2'+uy(1:nx,:).^2');
        clf;       
        colormap(viridis);
        imagesc(speed);
        caxis([0.01 0.17])
        hold on
        gap = 30; % Regula a desindade do campo vetorial, quanto maior o número, menos vetores e vice-versa
        quiver(2:gap:nx,1:gap:ny,ux(2:gap:nx,1:gap:ny)',uy(2:gap:nx,1:gap:ny)','k');        
        title(['Flow field steps: ',num2str(t),'; Reynolds number: ', num2str(Re),'; Omega: ', num2str(omega)] );
        colorbar
        drawnow
        fprintf("Ciclo: %d.\n", t)
    end
end

function res = vort(xvelocity, yvelocity)
res = zeros (size(xvelocity));
    for  x = 2:size(xvelocity,1)-1
        for  y = 2:size(yvelocity,2)-1
            res(x,y) = xvelocity(x,y+1)-xvelocity(x,y-1) - (yvelocity(x+1,y)-yvelocity(x-1, y));
        end
    end
end