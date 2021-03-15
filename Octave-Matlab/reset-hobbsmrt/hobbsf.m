function [f, g, r, J, H] = hobbsf(b,y)
% HOBBSF.M
% work out residual, Jacobian and Hessian of Hobbs problem at b
% assume data in y already
[m,kk] = size(y);
J=zeros(m,kk);
r=zeros(m,1);
% S is LOCAL 
S=zeros(3);
disp(S);
for i=1:m,
  if abs(0.1*b(3)*i)>50, 
       r(i) = 1e+35;
%      f = -1;
%      return;
% to indicate non-computable
      break;
  else   
    ee=exp(i*b(3)/10);
    tt=ee+10*b(2);
    qq=100*b(1)/(1+10*b(2)*exp(-0.1*b(3)*i))-y(i);
    r(i)=qq;
    J(i,1)=100*ee/tt;
    J(i,2)=-1000*b(1)*ee/tt^2;
    J(i,3)=10*i*b(1)*ee/tt-10*i*b(1)*ee^2/tt^2;
    S(1,2)=S(1,2)-qq*1000*ee/tt^2;
    S(1,3)=S(1,3)+qq*100*b(2)*i*ee/tt^2;
    S(2,2)=S(2,2)+qq*20000*b(1)*ee/tt^3;
    S(2,3)=S(2,3)+qq*(200*i*b(1)*ee^2/tt^3-100*i*b(1)*ee/tt^2);
    S(3,3)=S(3,3)+qq*(-2*i^2*b(1)*ee^2/tt^2+i^2*b(1)*ee/tt);
    S(3,3)=S(3,3)+qq*(-i^2*b(1)*ee^2/tt^2+2*i^2*b(1)*ee^3/tt^3);
  end;
end;
% disp(r); % Why doesn't this work inside fn
% if f<0, return;
% end;
S(2,1)=S(1,2);
S(3,1)=S(1,3);
S(3,2)=S(2,3);
H=2.*(J'*J+S);
fprintf("H");
disp(H);
f = r'*r;
g = 2*J'*r;
return % end of hobbsf.m
