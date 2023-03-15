# Python and Matlab program for SGMD

A simple program to implement the Symplectic geometry mode decomposition (SGMD), including python and matlab versions.

# References

[1] Pan H, Yang Y, Li X, et al. Symplectic geometry mode decomposition and its application to rotating machinery compound fault diagnosis[J]. Mechanical Systems and Signal Processing, 2019, 114:189–211. DOI: 10.1016/j.ymssp.2018.05.019.

[2] 潘海洋. 基于辛几何模态分解和支持矩阵机的机械故障诊断方法[D]. 湖南大学, 2019.

[3] https://zhuanlan.zhihu.com/p/66203573

# Additional information

1. More information about the program can be found at https://zhuanlan.zhihu.com/p/603813618.
2. The existing program certainly has unreasonable, negligent or even wrong places, please use discretion. And the provider is not responsible for any consequences caused by this program.
3. If there is something that can be improved or wrong, please do not hesitate to point out.
4. In addition, it should be noted that the program does not have the same result in its Matlab and python versions, and the matlab version offers fewer expansion parameters and functions than python.

# Example
![1](https://latex.codecogs.com/svg.image?%5Cbegin%7Bcases%7Dx(t)%20=%20%7Bx_1%7D(t)%20&plus;%20%7Bx_2%7D(t)%20&plus;%20%7Bx_3%7D(t)%5C%5C%7Bx_1%7D(t)%20=%202(1%20&plus;%200.5%5Csin%20(2%5Cpi%20t))%5Csin%20(60%5Cpi%20t)%5C%5C%7Bx_2%7D(t)%20=%20%5Csin%20(120%5Cpi%20t)%5C%5C%7Bx_3%7D(t)%20=%200.5%5Ccos%20(10%5Cpi%20t)%5Cend%7Bcases%7D)
fs = 5120 Hz, time = 1 s

The correlation coefficient threshold and NMSE threshold were 0.95 and 0.01 in Matlab, and 0.8, 0.001 in Python.

<div align=center>
<img src=https://user-images.githubusercontent.com/72395068/216860430-8d09f070-c548-4538-a1e2-d514c382a60b.png width=640/>
<div>Fig.1 The same simulation signal as in literature 1.</div>
</div>
<div align=center>
<img src=https://user-images.githubusercontent.com/72395068/216860595-a8c12df0-b32f-4528-a40a-175acecc5f1d.png width=640/>
<div>Fig.2 The decomposition results given in the SGMD paper.</div>
</div>
<div align=center>
<img src=https://user-images.githubusercontent.com/72395068/216865979-4f0ed357-5822-4f6d-9605-be5e7e548ed8.png width=640/>
<div>Fig.3 Decomposition results of this program (Python)</div>
</div>
<div align=center>
<img src=https://user-images.githubusercontent.com/72395068/216860769-21d85100-af65-4a48-8bb9-0254668163d4.png width=640/>
<div>Fig.4 Decomposition results of this program (Matlab)</div>
</div>
