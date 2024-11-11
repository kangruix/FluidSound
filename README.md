# FluidSound

Bubble-based water sound synthesis code based on the paper:
    
> [Improved Water Sound Synthesis using Coupled Bubbles](https://graphics.stanford.edu/papers/coupledbubbles/). Kangrui Xue, Ryan M. Aronson, Jui-Hsien Wang, Timothy R. Langlois, Doug L. James. *ACM Transactions on Graphics (SIGGRAPH North America 2023)*. 

## Build Instructions

**Dependencies:** C++11, Eigen 3.4.0

Building is handled by CMake. For example, to build from source on Mac & Linux:

    git clone https://github.com/kangruix/FluidSound
    cd FluidSound
    mkdir build && cd build && cmake ..
    make -j4

We provide an example scene in Scenes/GlassPour/. To run the code:

    ./runFluidSound ../Scenes/GlassPour/trackedBubInfo.txt 48000 1
    python ../scripts/write_wav.py output.txt 48000 

Afterwards, the simulated audio will be written to 'output.wav'. More scenes are available [here](https://graphics.stanford.edu/projects/waveblender/dataset/index.html).

### Citation

```bibtex
@article{xue2023coupledbubbles,
  author = {Xue, Kangrui and Aronson, Ryan M. and Wang, Jui-Hsien and Langlois, Timothy R. and James, Doug L.},
  title = {Improved Water Sound Synthesis using Coupled Bubbles},
  year = {2023},
  journal = {ACM Trans. Graph.},
}
```
