# FluidSound

Bubble-based water sound synthesis code based on the papers:
    
> [Improved Water Sound Synthesis using Coupled Bubbles](https://graphics.stanford.edu/papers/coupledbubbles/). Kangrui Xue, Ryan M. Aronson, Jui-Hsien Wang, Timothy R. Langlois, Doug L. James. *ACM Transactions on Graphics (SIGGRAPH North America 2023)*.  

> [Toward Animating Water with Complex Acoustic Bubbles](https://www.cs.cornell.edu/projects/Sound/bubbles/). Timothy R. Langlois, Changxi Zheng, Doug L. James. *ACM Transactions on Graphics (SIGGRAPH North America 2016)*. 

## Build Instructions

**Dependencies:** C++11, Eigen 3.4

Building is handled by CMake. For example, to build from source on Mac & Linux:

    git clone https://github.com/kangruix/FluidSound
    cd FluidSound
    mkdir build && cd build && cmake -DCMAKE_BUILD_TYPE=Release ..
    make -j4

We provide an example scene in Scenes/GlassPour/. To run the code:

    ./runFluidSound ../Scenes/GlassPour/trackedBubInfo.txt 48000 1
    python ../scripts/write_wav.py output.txt 48000 

(where the 1 indicates the scheme: 0 - uncoupled, 1 - coupled). Afterwards, the simulated audio will be written to 'output.wav'. More scenes are available [here](https://graphics.stanford.edu/projects/waveblender/dataset/index.html).

### Miscellaneous

Documentation can be built using Doxygen: `doxygen Doxyfile`

### Citation

```bibtex
@article{Xue:2023:CoupledBubbles,
  author = {Xue, Kangrui and Aronson, Ryan M. and Wang, Jui-Hsien and Langlois, Timothy R. and James, Doug L.},
  title = {Improved Water Sound Synthesis using Coupled Bubbles},
  journal = {ACM Transactions on Graphics (Proceedings of SIGGRAPH 2023)},
  year = {2023}, volume = {42}, number = {4}, month = jul
}
@article{Langlois:2016:Bubbles,
  author = {Langlois, Timothy R. and Zheng, Changxi and James, Doug L.},
  title = {Toward Animating Water with Complex Acoustic Bubbles},
  journal = {ACM Transactions on Graphics (Proceedings of SIGGRAPH 2016)},
  year = {2016}, volume = {35}, number = {4}, month  = jul
}
```
