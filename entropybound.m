d=2;%dimension
T=3;%the number of measurement bases
[A,MUB]=mub(d);
IU=[];%index of uncertainty
entropys=[];%optimal state-independent lower bound on Shannon entropy
entropy2=[];%optimal state-independent lower bound on Renyi-2 entropy
parfor n=1:1000
     B=eye(d);
%      for k=1:T-1 %---T random unitaries-
%          v=2*pi*rand(1,d^2-1).^3;v=reshape(v,[1,1,d^2-1]);%expansion parameters
%          u=expm(1i*sum(times(v,A),3));B=[B;u];
%      end
     for k=1:T-1  %T bases rotated from MUB
         v=2*pi*rand(1,d^2-1).^3;v=reshape(v,[1,1,d^2-1]);u=expm(1i*sum(times(v,A),3));
         B=[B;u*MUB(k*d+1:(k+1)*d,:)];
     end
     W=conj(B)*transpose(B);W=abs(W.^2);%-------overlap matrix--
     eig_W=sort(eig(W),'descend');eig_g=eig_W(2);%the first eigenvalue of view operator=the second eigenvalue of W
     IU=[IU T-eig_g];
     %-----------random pure states and upper bound on IC---------
      e1=T*log2(d);e2=T*log2(d);%temporary variables
     for m=1:5000 %minimize entropies over random pure states
         state=diag([1 zeros(1,d-1)]);
         v=2*pi*rand(1,d^2-1).^2;v=reshape(v,[1,1,d^2-1]);u=expm(1i*sum(times(v,A),3));
         state=transpose(conj(u))*state*u;
         p=diag(conj(B)*state*transpose(B));p=abs(transpose(p));%probability vector for measurement outcomes
         s=-p*transpose(log2(p+eps));%the sum of shannon entropies
         if s<e1
            e1=s;
         end
         p=reshape(p,d,T,1);
         s=sum(p.^2);s=-sum(log2(s));%the sum of Renyi-2 entropies
         if s<e2
            e2=s;
         end
     end
     entropy2=[entropy2 e2];  
     entropys=[entropys e1];
end
sz=30;figure,s=0:0.01:T-1;
scatter(IU,entropys,sz,[1,0.5,0],'.'), hold on;
%lower bound on Shannon entropy given by Eq.(16) of the main text//IC=T/d+(T-IU)(1-1/d)
plot(s,SentropyIC(T/d+(T-s)*(1-1/d),T),'b--','linewidth',1.13), hold on;
scatter(IU,entropy2,sz,[.55,0.55,.55],'.'), hold on;
%lower bound on Renyi-2 entropy given by Eq.(14) of the main text
plot(s,-T*log2((1/d+(T-s)/T/d)),'r','linewidth',1.13), hold on;

grid on;xlabel('Index of uncertainty');ylabel('Sum of entropies');
set(get(gca,'XLabel'),'FontSize',16);
set(get(gca,'YLabel'),'FontSize',16);
set(get(gca,'TITLE'),'FontSize',15);
set(gca,'fontsize',17);
grid on;grid minor;
