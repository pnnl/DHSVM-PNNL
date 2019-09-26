<h1>Distributed Hydrology Soil Vegetation Model (DHSVM) </h1>

This repository serves as the public source code repository of the Distributed Hydrology Soil Vegetation Model (DHSVM). You can find DHSVM documentation, and selected past and ongoing DHSVM-based research & projects on the <a href="http://dhsvmdev.pnl.gov///">DHSVM website </a>.

DHSVM (<a href="http://onlinelibrary.wiley.com/doi/10.1029/94WR00436/abstract">Wigmosta et al., 1994</a>) numerically represents with high spatial resolution (typically on the order of 100 m) the effects of local weather, topography, soil type, and vegetation on the hydrology within watersheds. The model is used to study the impacts of climate change, land use change, forest management practices, flooding, glacier dynamics, stream temperature and stream quality.

<strong>DHSVM is a research model that does not come with any warrantee or guarantee</strong>. Please be advised that no technical support is available other than the model web page. Because the model is under continous development, there is no guarantee that the newly developed modules or options are exhaustively tested or work properly. 

If you decide to use DHSVM, please acknowledge <a href="http://onlinelibrary.wiley.com/doi/10.1029/94WR00436/abstract">Wigmosta et al. [1994]</a> and any other relevant publications. We are very interested in receiving a copy of any manuscripts of studies in which the model is used. Finally, if you do find bugs in the model or if you have improvements to the model code, we are interested in incorporating your suggestions and/or contributions. 


## DHSVM v 3.2 ##
### _Release date: February 23, 2018_ ###

This is a major release from DHSVM 3.1.2. It includes several new features, function enhancements and bug fixes.<br />
The tutorial and sample data to run DHSVM v 3.2 will be made available on the <a href="http://dhsvmdev.pnl.gov//tutorials.stm">DHSVM website </a>.

__New Capabilities__
  * Variable radiation transmittance (with solar position and tree characteristics) 
  * Canopy gap (<a href="https://onlinelibrary.wiley.com/doi/abs/10.1002/hyp.13150">Sun et al., 2018</a>)
  * Snow sliding 
  * Python scripts to create stream network
  * Support of gridded meteorological data input

__Enhancement & Fixes__
  * Negative soil moisture 
  * Configuration and Build with CMake
<br />

## Configuration and Build with CMake ##

DHSVM and related utilities can be configured and built using
[CMake](https://cmake.org).  This provides an automated,
cross-platform way to locate and use system libraries (X11,
[NetCDF](http://www.unidata.ucar.edu/software/netcdf/), etc.) and
select optional features.  Here are some terse instructions: 

  * In the top DHSVM (where `CMakeLists.txt` is located), make a
    directory for the build, called `build` maybe.
    
  * In the `build` directory, run [CMake](https://cmake.org) with
    appropriate options, for example,
    
        cmake -D CMAKE_BUILD_TYPE:STRING=Release ..

    Look at `example_configuration.sh` for configurations used on
    several developers' systems. Alternatively, just use the script:

        sh ../example_configuration.sh

    This provides a vanilla configuration without X11,
    [NetCDF](http://www.unidata.ucar.edu/software/netcdf/), or RBM
    using the default C compiler.
    
  * If successful, build DHSVM and related programs, using

        cmake --build .

    The resulting executable programs will be in the build directory
    in a tree mirroring the source tree.  For example, DHSVM is
    `build/DHSVM/sourcecode/DHSVM`. 
    
The original Makefiles are in the source tree and can still be used as
described in the tutorial if preferred.

### Experimental Surface/Subsurface Mode ###

The normal DHSVM surface/subsurface routing scheme uses 4 neighbors
and directs cell outflow to all down-gradient neighbors ("D4").  Another
method, which uses 8 neighbors and directs outflow only to the
neighbor with the steepest gradient ("D8").  The latter has proved useful
in very low-gradient areas that span several cells.  To enable D8
routing add this option

    -D DHSVM_D8:BOOL=ON
    
to the configuration. D4 is the default.  

### Snow-only mode ###

If DHSVM is configured with this option,

   -D DHSVM_SNOW_ONLY:BOOL=ON
   
an additional executable is built, `DHSVM_SNOW`, which operates in
snow-only mode. 

