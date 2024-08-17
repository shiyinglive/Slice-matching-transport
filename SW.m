% scaled Sliced Wasserstein distance between two images
function sw = SW(I0,I1)
I_domain = [0, 1];
Ihat_domain = [0, 1];
theta_seq = 0:4:179;
rm_edge = 1;
 I0hat = RCDT(I_domain, I0, Ihat_domain, theta_seq, rm_edge);
 I1hat = RCDT(I_domain, I1, Ihat_domain, theta_seq, rm_edge);

 sw = sqrt(sum(sum((I0hat-I1hat).^2)));

 %%% RCDT function

function rcdt_mat = RCDT(x_I, I, x_c, theta_seq, rm_edge)

eps=1e-8;
pR = radon(I,theta_seq);

for th=1:length(theta_seq)
    p = pR(:,th); p = p/sum(p);
    
    x = linspace(x_I(1),x_I(2),length(p));
    x_cdt=linspace(x_c(1),x_c(2),length(p));
    pCDT(:,1)=CDT(x, p+eps, x_cdt, rm_edge);
       
    rcdt_mat(:,th) = pCDT;
end


             
end
%%% end of RCDT function

%%%% CDT function
function [s_cdt]=CDT(x,J,x_cdt,rm_edge)

if (size(J,2) == 1) 
    J=J';
end
if (size(x,2) == 1) 
    x=x';
end
if (size(x_cdt,2) == 1) 
    x_cdt=x_cdt';
end

J=J/sum(J); 
cJ=cumsum(J);

s_cdt=interp1(cJ,x,x_cdt,'pchip');

if rm_edge
    s_cdt(1) = [];
    s_cdt(end) = [];
end
 
end
%%%% end of CDT function
end