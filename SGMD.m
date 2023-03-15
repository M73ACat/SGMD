function sgc = SGMD(x,fs,tau,threshold_corr,threshold_nmse,varargin)
% Author: M73ACat
% 2023/03/15
% Basic parameters: 
% x                 list, input a signal x,
% fs                int, sample frequency of x,  
% threshold_corr    float, the threshold of the corr (optional, default 0.95),
% threshold_nmse    float, the threshold of the nmse (optional, default 0.95),

% Advanced parameters: 
% tau               int, delay time (optional, defualt 1)
% nfft              int, the number of FFT points used to calculate the PSD estimate (optional, default 0 -> 2^floor(log2(n)) ),
% PSD_function      char, functions for calculating PSD (optional ['0'(means 'periodogram'),'1'(means 'pwelch')], default '0' -> periodogram)

% Reference:
% [1] 潘海洋. 基于辛几何模态分解和支持矩阵机的机械故障诊断方法[D]. 湖南大学, 2019.
% [2] Pan H, Yang Y, Li X, et al. Symplectic geometry mode decomposition and its application to rotating machinery compound fault diagnosis[J]. Mechanical Systems and Signal Processing, 2019, 114:189–211. DOI: 10.1016/j.ymssp.2018.05.019.
% [3] https://zhuanlan.zhihu.com/p/66203573

tel=10^-3;
[n, row]=size(x);
if n<row
    x=x';
    n=row;
end
PSD_function = '0';
nfft = 2^floor(log2(n));
for i=1:length(varargin)
    if ischar(varargin{i})
        PSD_function = varargin{i};
    else
        if varargin{i} > 0; nfft = varargin{i}; end
    end
end
if PSD_function == '0'
    [a,f]=periodogram(x,[],'twosided',nfft,fs);
else
    window = hanning(nfft);
    [a,f] = pwelch(x,window,nfft*(3/4),nfft,fs,'twoside');
end
[~,ra]=max(a(1:end/2));
if f(ra)/fs<tel
    d=floor(n/3);
else
    d=floor(1.2*fs/f(ra));
end
[X,m]= Trajectorymatrix(x,d,tau);
A=X'*X;

[Q,~] = eig(A*A);

Q = real(Q);
sgc=zeros(n,d);
for nn=1:d
    z=Q(:,nn)*Q(:,nn)'*X';
    y=zeros(n,1);

    newd=min(size(z));
    newm=max(size(z));
    if m>d
        z1=rot90(z',3);
    else
        z1=rot90(z,3);   
    end
    for k=1:n
        if k<=newd
            y(k)=sum(diag(z1,newm-k))/k;
        elseif k>newd&&k<=newm
            y(k)=sum(diag(z1,newm-k))/newd;
        else
            y(k)=sum(diag(z1,newm-k))/(n-newm+1);
        end
    end
    sgc(:,nn)=y;
end

if std(sgc(:,end)) > std(sgc(:,1))
    y = fliplr(sgc);
else
    y = sgc;
end
indexs = linspace(1,d,d);
flags = logical(indexs);
x_e = sum((x-mean(x)).^2);
sgc = [];
g_h = 0;
g_h_e = 0;
while ~isempty(indexs(flags))
    temp_index = indexs(flags);
    source = y(:,temp_index(1));
    flags(temp_index(1)) = 0;
    temp_index = temp_index(2:end);
    temp_flag = [];
    for i = 1: length(temp_index)
        corrs = corrcoef(source,y(:,temp_index(i)));
        corrs = corrs(1,2);
        if corrs >= threshold_corr
            temp_flag(end+1) = i;
        end
    end
    flags(temp_index(temp_flag)) = 0;
    sgc(:,end+1) = source+sum(y(:,temp_index(temp_flag)),2);
    g_h = sum(sgc,2);
    g_h_e = sum((x-g_h).^2);
    if g_h_e / x_e < threshold_nmse
        break
    end
end
sgc(:,end+1) = x - g_h;

end


function [t,m] = Trajectorymatrix(x,d,tau)
%根据Takens嵌入定理构建的轨迹矩阵
%input
%x          input a signal x,
%d          embedding dimension
%tau        delay time
%output
%t m*n Trajectory matrix
%20180611 PhD. Jin Hang   in Chengdu
[n ,c]=size(x);
if n==1 || c==1
    if c>n
        x=x';
        [n ,~]=size(x);
    end
else
    error('x isn`t  Column vector')
end
m=n-(d-1)*tau;

if m<=0
    error('Index out of matrix,Replace d or Tau')
end

t=zeros(m,d);

for kk=1:d
    t(1:m,kk)=x(1+(kk-1)*tau:m+(kk-1)*tau);
end

end