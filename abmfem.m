function s=abmfem(X,Ind)

% s = abmfem(X,Ind)
%
% ABMFEM Algorithm computes the posterior probabilities of two classes at
% data points in X without any required label vector. Ind is the indices of
% the individual sample, default is Ind = 1:N . N is the number of points
% in X.
%
% If you use this code in your work, please cite:
%
% Köktürk, Ba?ak Esin, and Bilge Karaçal?. "Annealing-based model-free expectation maximisation for
% multi-colour flow cytometry data clustering." International Journal of Data Mining and Bioinformatics 14,
% no. 1 (2016): 86-99.


N = length(X);
if (nargin<2)
    Ind = 1:N;
end

r=randperm(size(X,1));
m=size(X,1)/2;
r=r(1:m);
s.Y = ones(1,size(X,1));
s.Y(r) = 2;
iter = 0;
prev_Y = zeros(1,N);

D=zeros(N,N);
for k=1:N
    D((k+1):N,k)=(sum((X((k+1):N,:)-ones(N-k,1)*X(k,:)).^2,2));
    D(k,(k+1):N)=D((k+1):N,k)';
end

s.Y_vector = zeros(100,N);

s.error_vector = zeros(100,1);

s.P1_vector = zeros(100,N);


a = 0;

for n =100:-1:1
    prev_Y = zeros(1,N);
    iter = 0;
    while(sum(s.Y~=prev_Y)>N*0.01)
        a = a+1;
        iter = iter+1;
        prev_Y = s.Y;
        
        [L,P1,P2] = qsl(X,s.Y,n,D);
        
        s.Y =ones(1,N);
        ind = find(P1<0.5);
        s.Y(ind) = 2;
        
        P_C0 = length(find(s.Y==1))/N;
        P_C1 = length(find(s.Y==2))/N;
        
        P1 = (P1*P_C0)./(P1*P_C0+P2*P_C1);
        P2 = (P2*P_C1)./(P1*P_C0+P2*P_C1);
        
        s.Y =ones(1,N);
        ind = find(P1<0.5);
        s.Y(ind) = 2;
        
        if(length(find(s.Y==1))==0 || length(find(s.Y==2))==0)
            r=randperm(size(X,1));
            m=size(X,1)/2;
            r=r(1:m);
            s.Y = ones(1,size(X,1));
            s.Y(r) = 2;
        end
        if (iter>25)
            break;
        end
        s.full_y(a,:)=s.Y;
        
        s.indices(a) = n;
        
    end
    s.Y_vector(n,:)=s.Y;
    s.cost_vector(n)=4*sum(P1.*P2)+2*n;
    s.P1_vector(n,:)=P1;
end


[~,b]=min(s.cost_vector);
s.Y = s.Y_vector(b,:);
P1 = s.P1_vector(b,:);
P2=1-P1;

N1= length(find(s.Y==1));
N2 = N-N1;
C11 = sum(P1(s.Y==1))/N1;
C12 = sum(P2(s.Y==1))/N1;
C21 = sum(P1(s.Y==2))/N2;
C22 = sum(P2(s.Y==2))/N2;
s.error = sum(C21+C12)/sum(C11+C12+C21+C22);
ind = find(P1<=0.5);
s.ind1 = Ind(ind);
s.data1 = X(ind,:);
ind = find(P1>0.5);
s.ind2 = Ind(ind);
s.data2 = X(ind,:);
end
