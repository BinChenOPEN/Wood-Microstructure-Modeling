# Wood-Microstructure-Modeling
## Summary
This project generates realistic wood (Birch, Spruce, etc.) microstructures from a set of given parameters 
This project includes two parts：
1. **Birch microstructure generation;**
2. **Spruce microstructure generation;**

## Key features in the generated microstructure
- Ray cells are included in the wood structure. For birch, each group includes two columns of ray cells. Regarding to spruce, each group includes a signle column of ray cells.  
- The ray cells will distort the fibers nearby, to mimic the real wood growth.  
- The impact of the ray cells in contraint into a small region.  
- Randomness is applied on fibers, vessels and ray cells.  
- The fibers, vessels are distorted by distortion map, which influence only the fibers nearby.
- In spruce, the difference in early wood and late wood are simulated. The late fibers are compressed more, with thicker cell wall.

## Note  
The generation is time consumming for large wood structure generation. Parallel computation uisng multi-core CPU is used to speed up the computation. We tried to reduce the requirement of the memory size. In our application, 3000^3 voxels volume image can be generated with 32GB RAM. Larger structure needs more memory. The computaion also requires saving and loading large amount of data. Therefore, SSD is recommended.

## Usage
Download the project and save it into a selected folder. 
Read the notation of the parameters in the Matlab scripts. 
If you want to quickly run the projects. **Here are some notes**:  
`cellWallThick` is should not be too small. We recommend to use value larger than 3  
`sizeVolume` should not be too small in order to include the ray cell features into the volume. I recommend to start with 1500 $\times$ 1500 $\times$ 750 as a start.  
`cellR` is an important variable. It is roughly the radius of the fiber.  
`extraSZ` works good in most of cases with the default values. If the obtained final structure in the final subfolder  `FinalVolumeSlice` includes some white boundary, you can increase this value a bit. 
`isExistRayCell` is used to control whether there is any ray cell in the structure.
`isExistVessel` is used to control whether there is any vessel in the structure.

## Results
Here it is a generated microstructure of birch  
![Birch](https://github.com/cbbuaa/Wood-Microstructure-Modeling/assets/23333414/5cc4806e-52ec-46e0-888a-ccaa9a74a855)  
Here it is a generated microstructure of Sprucce  
![image](https://github.com/cbbuaa/Wood-Microstructure-Modeling/assets/23333414/15b9ab1c-c080-4141-93d7-3b4862555be2)  


## Citation
[Chen, Bin, et al. "A distortion-map-based method for morphology generation in multi-phase materials-application to wood." Composites Science and Technology 244 (2023): 110262.](https://www.sciencedirect.com/science/article/pii/S0266353823003561)
