WMAN
==
Description                
--
>WNANet is an erratic noise attenuation method using DNN waveform attention.  

Reference
==
If you find this package useful, please do not forget to cite the following paper.
>[1]C. Fu, Y. Huo, G. Li and Y. Chen, "Unsupervised Learning With Waveform Multibranch Attention Mechanism for Erratic Noise Attenuation," in IEEE Transactions on Geoscience and Remote Sensing, vol. 62, pp. 1-11, 2024, Art no. 5933711, doi: 10.1109/TGRS.2024.3487306.<br>
[2]Y. Chen, 2020,“Fast dictionary learning for noise attenuation of multidimensional seismicdata,”Geophysical Journal International, vol. 221, no. 3, pp. 1717–1727, 


BibTeX:
==
>[1]@ARTICLE{CWMANet,
  author={Fu, Chao and Huo, Yupeng and Li, Guodong and Chen, Yangkang},
  journal={IEEE Transactions on Geoscience and Remote Sensing}, 
  title={Unsupervised Learning With Waveform Multibranch Attention Mechanism for Erratic Noise Attenuation}, 
  year={2024},
  volume={62},
  pages={1-11},
  doi={10.1109/TGRS.2024.3487306}
  }<br>
  [2]@article{chen2020fast,
  title={Fast dictionary learning for noise attenuation of multidimensional seismic data},
  author={Chen, Yangkang},
  journal={Geophysical Journal International},
  volume={222},
  number={3},
  pages={1717--1727},
  year={2020},
  publisher={Oxford University Press}
}

The accepted version of the paper can be downloaded from:
--
>https://ieeexplore.ieee.org/document/10737129

Copyright
==
>Developers of the WMAN package, 2024-present

# Environment Configuration
This repository contains code that uses several Python libraries. Below are the specific versions of the dependencies used in this project:
## Dependencies
| Library | Version |
|---------|---------|
| MatplotLib | 3.6.2 |
| NumPy | 1.23.5 |
| Pandas | 1.5.2 |
| PyTorch | 1.13.0+CU117 |
| Scikit-learn | 1.2.0 |
| TQDM | 4.64.1 |
| SciPy | 1.9.3 |


# Install
==
## This environment uses CUDA 11.7 with PyTorch. Make sure your system has compatible NVIDIA drivers installed.
conda create -n wman python=3.11.7 <br>
conda activate wman <br>
conda install ipython notebook <br>
pip install matplotlib==3.6.2 <br>
pip install numpy==1.23.5 <br>
pip install pandas==1.5.2 <br>
pip install torch==1.13.0+cu117 <br>
pip install scikit-learn==1.2.0 <br>
pip install tqdm==4.64.1 <br>
pip install scipy==1.9.3 <br>


Then install WNAN using the latest version
>git clone https://github.com/furser1/WMAN <br>
cd wman <br>
pip install -v -e <br>

Code running 
==
Data description: The uploaded file is divided into the following parts:
--
> 1.1 Data part: corresponding to synthesized data and actual data respectively
* syn2d.mat<br>
* syn3d.mat<br>
* real2d.mat<br>
* real3d.mat<br>
> 1.2 Pytorch network part: corresponds to WMANet network and PATCHUNET network respectively
* model_wma_2dsyn1.ipynb<br>
* model_wma_3dsyn2.ipynb<br>
* model_wma_2dreal1.ipynb<br>
* model_wma_3dreal2.ipynb<br>
* model_patch_2dsyn1.ipynb<br>
* model_patch_3dsyn2.ipynb<br>
* model_patch_2dreal1.ipynb<br>
* model_patch_3dreal2.ipynb<br>
> 1.3 Data presentation part: used for data processing and presentation, corresponding to
* WMANet_syn1_git.m 
* WMANet_syn2_git.m
* WMANet_real1d_git.m
* WMANet_real3d_git.m
  
 Usage
---
>  First, You can use 1.3 to load the data from 1.1 and Patching processing<br>
>  Then, You can use 1.2's network for denoising<br>
>   Finally, You can use 1.3 data Unpatching part to recover the data and display.<br>


