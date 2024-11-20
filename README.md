CWMANet
==
Description                
--
>CWNANet is an erratic noise attenuation method using DNN waveform attention.  

Reference
==
If you find this package useful, please do not forget to cite the following paper.
>C. Fu, Y. Huo, G. Li and Y. Chen, "Unsupervised Learning With Waveform Multibranch Attention Mechanism for Erratic Noise Attenuation," in IEEE Transactions on Geoscience and Remote Sensing, vol. 62, pp. 1-11, 2024, Art no. 5933711, doi: 10.1109/TGRS.2024.3487306.

BibTeX:
==
>@ARTICLE{CWMANet,
  author={Fu, Chao and Huo, Yupeng and Li, Guodong and Chen, Yangkang},
  journal={IEEE Transactions on Geoscience and Remote Sensing}, 
  title={Unsupervised Learning With Waveform Multibranch Attention Mechanism for Erratic Noise Attenuation}, 
  year={2024},
  volume={62},
  pages={1-11},
  doi={10.1109/TGRS.2024.3487306}
  }

The accepted version of the paper can be downloaded from:
--
>https://ieeexplore.ieee.org/document/10737129

Copyright
==
>Developers of the CWMAN package, 2024-present

Install
==
conda create -n cwman python=3.11.7 <br>
conda activate cwman <br>
conda install ipython notebook <br>
pip install torch==2.2.1 <br>

Then install CWNAN using the latest version
>git clone https://github.com/furser1/cwman <br>
cd cwman <br>
pip install -v -e <br>
