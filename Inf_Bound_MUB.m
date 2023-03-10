d=2;%dimension
T=3;%the number of measurement bases
[A,MUB]=mub(d);MUB(T*d+1:(d+1)*d,:)=[];% T MUBs

% B=eye(d); %---T random unitaries-----------
% for k=1:T-1
%     v=2*pi*rand(1,d^2-1).^3;v=reshape(v,[1,1,d^2-1]);%expansion parameters
%     u=expm(1i*sum(times(v,A),3));B=[B;u];
% end

B=eye(d);
for k=1:T-1  %T bases rotated from MUB
    v=2*pi*rand(1,d^2-1).^11;v=reshape(v,[1,1,d^2-1]);u=expm(1i*sum(times(v,A),3));
    B=[B;u*MUB(k*d+1:(k+1)*d,:)];
end
W=conj(B)*transpose(B);W=abs(W.^2);         %overlap matrix--
eig_W=sort(eig(W),'descend');eig_g=eig_W(2);%the second eigenvalue of W=the first eigenvalue of view operator
%-------RANDOM STATE
mixity=[];inf_gain=[];
parfor j=1:10000
  state=diag(rand(1,d).^20);state=state/trace(state);
  s=trace(state^2)-1/d;
  %mixity of random state
  mixity=[mixity s];                
  %random unitary
  v=2*pi*rand(1,d^2-1).^3;v=reshape(v,[1,1,d^2-1]);u=expm(1i*sum(times(v,A),3));
  state=transpose(conj(u))*state*u;
  p=diag(conj(B)*state*transpose(B));
  %-------information gain
  p=p-1/d;
  inf_gain=[inf_gain sum(p.^2)];     
end
%-------FIGURE
 sz=30;figure,
 scatter(mixity,inf_gain,sz,[.6,1,.6],'.'), hold on;
 s=0:0.001:1-1/d;
 plot(s,eig_g*s,'r--','linewidth',1.2), hold on;  %our upper bound on information gain given in Eq.(8) of the main text
   
 grid on;xlabel('mixity');ylabel('Information gain');
 axis([0,1-1/d,0,T]);
 set(get(gca,'XLabel'),'FontSize',16);
 set(get(gca,'YLabel'),'FontSize',16);
 set(get(gca,'TITLE'),'FontSize',15);  
 set(gca,'fontsize',17);
 grid on;grid minor; 

disp(transpose(eig_W));
