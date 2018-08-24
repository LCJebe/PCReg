%L1PCA Polynomial-time, exact L1-principal component analysis.
%  [Qopt, Bopt]=L1PCA(X, K) produces a unitary matrix Qopt, with 
%  as many rows as X (say, D) and K columns (D>=K), and a binary
%  matrix Bopt (each entry is either +1  or -1) of as many rows as
%  the columns of X (say, N) and K columns, so that 
%  (i) Qopt=argmax_{Q'*Q=eye(K)} sum(sum(abs(X'*Q))) and
%  (ii) Bopt=argmax_{each entry of B is either +1 or -1} sum(svd((X*B)).
%  The columns of Qopt are the K L1-principal components of X.
%  This function implements the polynomial-time exact algorithm,
%  presented in the articles referenced below.
%
%   Author:
%   Dr. Panos P. Markopoulos
%   Electrical Engineering Department
%   Rochster Institute of Technology
%   Email: pxmeee@rit.edu
%   Web: https://people.rit.edu/pxmeee/
%   
%   ACKNOWLEDGEMENTS
%   Functions authored while Prof. Panos Markopoulos
%   was a PhD student at Univeristy at Buffalo, working
%   under the guidance of Prof. Dimitris Pados 
%   (now with FAU; dpados@fau.edu).
%   
%   Prof. George Karystinos (TUC, GREECE; karystinos@telecom.tuc.gr) also
%   contributed to this code.
%
%   Date: June 2013
%
%   References: 
%   1) P. P. Markopoulos, G. N. Karystinos, and D. A. Pados,"Optimal
%   algorithms for L1-subspace signal processing," IEEE Transactions on
%   Signal Processing, vol. 62, pp. 5046-5058, Oct. 2014. 
%   2) P. P. Markopoulos, G. N. Karystinos, and D. A. Pados,"Some options
%   for L1-subspace signal processing," in Proc. 10th International
%   Symposium on Wireless Communication Systems (ISWCS 2013), 
%   Ilmenau, Germany, Aug. 2013, pp. 622-626.
%   
%   **Inquiries regarding the script provided below are cordially welcome. 
%   In case you spot a bug, please let me know. 
%   If you use some piece of code for your own work, please cite the 
%   articles above.** 
 


function [Qopt, Bopt]=l1pca(X, K)

tic;
toler=1e-8;

if norm(imag(X),2)>toler
    error(['X must be a real matrix (here, norm(imag(X),2)=' num2str(norm(imag(X),2)), ').']')
end

[D, ~]=size(X);
[U, S, ~]=svd(X);
s=diag(S);
d=sum(s>toler);

if K>d
    error(['K must be less or equal to d=rank(X) (here, d=', int2str(d), ', K=', int2str(K), ').']);
end


% num_of_cands_exhaust=nchoosek(2^(N-1)+K-1,K); 
% w=zeros(1,d);
% for i=0:d-1
%     w(i+1)=nchoosek(N-1,i);
% end
% num_of_cands_poly=nchoosek(sum(w)+K-1,K);

Q=X;
if d<D
    Q=U(:,1:d)'*X;
end
[~, Bc]=compute_candidates(Q', d);


[~, a]=size(Bc);
 
muts=multisets(1:a,K);
[n_muts, ~]=size(muts);
 
mopt=-inf;
for i=1:n_muts
    
    comb=muts(i,:);
    B=Bc(:,comb);
    
    m=sum(svd(X*B,'econ'));
    if m>mopt
        mopt=m;
        Bopt=B;
    end
end
[U, ~, V]=svd(X*Bopt);
Qopt=U(:,1:K)*V';
timelapse=toc;
maxmetric=mopt;
 
disp('------------------------------');
disp(['Number of cadidates checked: ' int2str(n_muts)]);
disp(['Time elapsed (sec): ' num2str(timelapse)]);
disp(['Metric value: ' num2str(maxmetric)]);
disp('------------------------------');




function [C, B]=compute_candidates(Q, Dor)

[N, D]=size(Q);
if D>2
    ncomb=nchoosek(N,D-1);    
    combinations=nchoosek(1:N,D-1);
    B=zeros(N,ncomb);
    C=zeros(Dor,ncomb);    
    for i=1:ncomb
        I=combinations(i,:);
        Qbar=Q(I,:);       
        [~, ~, Z]=svd(Qbar);
        c=Z(:,end);    
        c=c*mysign(c(D));      
        cp=c;
        if D<Dor
            cp=[c; zeros(Dor-D,1)];
        end
        C(:,i)=cp;      
        B(:,i)=mysign(Q*c);
        for d=1:D-1           
            [~, ~, Z]=svd(Qbar([1:d-1 d+1:D-1],1:D-1));
            c=Z(:,end);           
            c=c*mysign(c(D-1));           
            B(I(d),i)=mysign(Qbar(d,1:end-1)*c);
        end   
    end  
    [Cn, Xn]=compute_candidates(Q(:,1:D-2),Dor);
    C=[C Cn];
    B=[B Xn];   
elseif D==2
    B=ones(N,N);
    C=ones(Dor,N);
    for i=1:N
        [~,~,Z]=svd(Q(i,:));
        c=Z(:,2);        
        c=c*mysign(c(2));
        cp=[c; zeros(Dor-2,1)];
        B(:,i)=mysign(Q*c);
        B(i,i)=mysign(Q(i,1));
        C(:,i)=cp;
    end
else
    B=mysign(Q);
    C=[1; zeros(Dor-1,1)];
end


function z=mysign(A)
toler=1e-8;
z=sign(A);
z(abs(z)<=toler)=1;


function M=multisets(S,K)
if min(size(S))>1
     error('S must be a N by 1 array.');
end
M=S(mulst(length(S),K));
if K<2
    M=M';
end

function x=mulst(n,k)
if k>1
    x=zeros(nchoosek(n+k-1,k),k);
    cnt=1;
    for i=1:n
        tmp=nchoosek(n-i+1+k-1-1,k-1);
        x(cnt:cnt-1+tmp,:)=[repmat(i,tmp,1) mulst(n-i+1,k-1)+i-1];
        cnt=cnt+tmp;
    end
elseif k==1
    x=(1:n)';
end