""" 
Author: M73ACat
2023/02/06

Reference: 
[1] 潘海洋. 基于辛几何模态分解和支持矩阵机的机械故障诊断方法[D]. 湖南大学, 2019.
[2] Pan H, Yang Y, Li X, et al. Symplectic geometry mode decomposition and its application to rotating machinery compound fault diagnosis[J]. Mechanical Systems and Signal Processing, 2019, 114:189–211. DOI: 10.1016/j.ymssp.2018.05.019.
[3] https://zhuanlan.zhihu.com/p/66203573
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.linalg as lg

class SGMD:
    def __init__(self,sig,fs,mode='eig',delay_time=1,nfft=256,is_plot=False,threshold_corr=0.95,threshold_nmse=0.25, end_indicator='nmse') -> None:
        """ 
        Basic parameters: 
            sig: ndarry, the target signal,
            fs: int, sample frequency of sig,
            threshold_corr: float, the threshold of the corr (optional, default 0.95),
            threshold_nmse: float, the threshold of the nmse (optional, default 0.5).
        Advanced parameters: 
            mode: str, 'eig', 'schur' or 'qr' (optional, default eig),
            delay_time: int, one parameter for SGMD.trajectory_matrix (optional, default 1),
            nfft: int, window length for SGMD.trajectory_matrix (optional, default 256),
            is_plot: bool, whether to plot the PSD figure (optional, default False),
            end_indicator: str, metrics for deciding whether to stop iterating, 'nmse', 'std' or 'var' (optional, default nmse).
        """
        self.sig = sig
        self.n = len(sig)
        self.fs = fs
        self.nfft = nfft
        self.threshold_corr = threshold_corr
        self.threshold_nmse = threshold_nmse
        self.mode = mode
        self.end_indicator = end_indicator
        self.trajectory_matrix(is_plot=is_plot,delay_time=delay_time)

    def trajectory_matrix(self,is_plot=False,delay_time=1):
        if 'X' in dir(self):
            if is_plot:
                plt.plot(self.f1,self.psd)
                plt.xscale('log')
                plt.yscale('log')
                plt.xlabel('f/Hz',fontsize=12)
                plt.ylabel('PSD', fontsize=12)
                plt.grid(False)
                plt.plot(self.f1[self.max_index[0]],self.psd[self.max_index[0]],'*r')
                plt.show()
            return self.X
        else:
            psd,f1 = plt.psd(self.sig,
                                NFFT=self.nfft,
                                Fs=self.fs,
                                detrend='mean',
                                window=np.hanning(self.nfft),
                                noverlap=int(self.nfft*3/4),
                                sides='twosided')
            plt.close()
            self.f1 = f1[self.nfft//2:]
            self.psd = psd[self.nfft//2:]
            self.max_index = np.where(self.psd[:]==np.max(self.psd))
            f_max = self.f1[self.max_index[0]][0]
            norm_val = f_max / self.fs
            self.d = 1.2*(self.fs/f_max) if norm_val >= 1e-3 else self.n/3
            self.d = np.int16(np.round(self.d))
            self.m = int(self.n - (self.d-1)*delay_time)
            self.X = np.zeros((self.m,int(1+(self.d-1)*delay_time)))
            if self.d > self.m:
                for i in range(self.m):
                    self.X[i,:] = self.sig[i:i+(self.d-1)*delay_time+1]
            else:
                for i in range(self.d):
                    self.X[:,i] = self.sig[i*delay_time:self.m+i*delay_time]

    def sgmd(self):
        matrix_A = np.dot(self.X.T,self.X)
        if self.mode == 'eig':
            _, matrix_Q = lg.eig(np.dot(matrix_A,matrix_A))
        elif self.mode == 'schur':
            _, matrix_Q = lg.schur(np.dot(matrix_A,matrix_A),'complex')
        elif self.mode == 'qr':
            matrix_Q, _ = lg.qr(np.dot(matrix_A,matrix_A))
        matrix_Q = np.real(matrix_Q)
        matrix_Y = np.zeros((self.d,self.n))
        self.d, self.m = min(self.d,self.m), max(self.d,self.m)
        for i in range(self.d):
            matrix_Z = np.dot(np.dot(np.expand_dims(matrix_Q[:,i],-1),np.expand_dims(matrix_Q[:,i],-1).T),self.X.T)
            matrix_Z = matrix_Z.T if self.m < self.d else matrix_Z
            matrix_Y[i,:self.d] = [np.mean(np.diag(np.flip(matrix_Z[:j+1,:j+1],axis=1))) for j in range(self.d)]
            matrix_Y[i,self.d:self.m] = [np.mean(np.diag(np.flip(matrix_Z[:,j+1-self.d:j+1],axis=1))) for j in range(self.d,self.m)]
            matrix_Y[i,self.m:] = [np.mean(np.diag(np.flip(matrix_Z[j+1-self.m:,j+1-self.d:],axis=1))) for j in range(self.m,self.n)]
                
        index = np.arange(self.d)
        flags = np.array([True]*(self.d))

        if self.end_indicator == 'nmse':
            x_e = np.sum((self.sig-np.mean(self.sig))**2)
        else: 
            x_e = np.std(self.sig) if self.end_indicator == 'std' else np.var(self.sig)
        while len(index[flags]):
            source = matrix_Y[index[flags][0]]
            flags[index[flags][0]] = False
            tmp_flag = [i for i in index[flags] if np.corrcoef(source,matrix_Y[i])[0,1] >= self.threshold_corr]
            flags[tmp_flag] = False
            SGCs = np.vstack((SGCs,source+np.sum(matrix_Y[tmp_flag],axis=0))) if 'SGCs' in dir() else source+np.sum(matrix_Y[tmp_flag],axis=0)
            if self.end_indicator == 'nmse':
                g_h = np.sum(SGCs,axis=0)
                g_h_e = np.sum((self.sig-g_h)**2)
            else: 
                g_h = np.sum(matrix_Y[index[flags]],axis=0)
                g_h_e = np.std(g_h) if self.end_indicator == 'std' else np.var(g_h)
            if g_h_e / x_e < self.threshold_nmse:
                break
        if self.end_indicator == 'nmse':
            SGCs = np.vstack((SGCs,self.sig-g_h))
        else:
            SGCs = np.vstack((SGCs,g_h))
        return SGCs


if __name__ == '__main__':
    fs = 5120
    time = 1
    t = np.arange(0,time,1/fs)
    x1_t = 2*(1+0.5*(np.sin(2*np.pi*t)))*np.sin(60*np.pi*t)
    x2_t = np.sin(120*np.pi*t)
    x3_t = 0.5*np.cos(10*np.pi*t)
    sig = x1_t + x2_t + x3_t

    nfft = 256
    sgmd = SGMD(sig,fs,nfft=nfft,threshold_corr=0.8,threshold_nmse=0.001,mode='eig')

    SGCs = sgmd.sgmd()

    plt.rc('font',family='Times New Roman',size='12')

    plt.figure()
    plt.title('signals')
    plt.subplot(4,1,1)
    plt.plot(sig)
    plt.xlim(0,len(sig))
    plt.xlim([0,fs*t[-1]])
    plt.subplot(4,1,2)
    plt.plot(x1_t)
    plt.xlim([0,fs*t[-1]])
    plt.subplot(4,1,3)
    plt.plot(x2_t)
    plt.xlim([0,fs*t[-1]])
    plt.subplot(4,1,4)
    plt.plot(x3_t)
    plt.xlim([0,fs*t[-1]])

    plt.figure()
    plt.title('SGCs From SGMD')
    for i in range(SGCs.shape[0]):
        plt.subplot(int(SGCs.shape[0]),1,i+1)
        plt.plot(SGCs[i])
        plt.xlim([0,fs*t[-1]])
    plt.show()