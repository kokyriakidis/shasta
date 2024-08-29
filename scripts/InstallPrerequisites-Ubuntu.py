#!/usr/bin/python3

import os
import subprocess
import tempfile

def isArm():
    return subprocess.getoutput("uname -p") == "aarch64"

def runCommand(command):
    if(os.system(command)):
        raise Exception("Error running command: " + command)
        
def installPackage(package):
    runCommand("sudo apt-get install --assume-yes " + package)

def installAptPackages():
    packages = [
    "git",
    "g++",
    "make",
    "cmake",
    "libboost-system-dev", 
    "libboost-program-options-dev",
    "libboost-graph-dev",
    "libboost-chrono-dev",
    "libpng-dev", 
    "libblas-dev", 
    "liblapack-dev",
    "gfortran",
    "libseqan2-dev",
    "ncbi-blast+",
    "graphviz",
    "gnuplot",
    "python3-dev", 
    ]
    runCommand("sudo apt-get install --assume-yes " + " ".join(packages))
    


def installPybind11():
    try:
        # This works for Ubuntu 22.04 and older
        print("Attempting pybind11 installation using pip3")
        installPackage("python3-pip")
        runCommand("sudo pip3 install pybind11")
    except:
        # This works for Ubuntu 24.04 (and newer, presumably).
        print("Pybind11 installation using pip3 did not work, installing using apt-get.")
        installPackage("python3-pybind11")


        
def installSpoa():
    # The spoa library is available in the stable Ubuntu repository, but
    # without the static version.
    # So we have to build it from source.
    
    with tempfile.TemporaryDirectory() as temporaryDirectory:
        print("Building spoa library using temporary directory", temporaryDirectory)
        
        # Change to the temporary directory.
        oldDirectory = os.getcwd()
        os.chdir(temporaryDirectory)
                
        # Get the code.
        runCommand("sudo apt-get install curl")
        runCommand("curl -L https://github.com/rvaser/spoa/archive/refs/tags/4.0.8.tar.gz -o 4.0.8.tar.gz")
        runCommand("tar -xvf 4.0.8.tar.gz")
    
        # Set spoa build flags.
        if isArm():
            spoaBuildFlags="-DCMAKE_BUILD_TYPE=Release -Dspoa_build_tests=OFF"
        else:
            # The spoa dispatcher feature selects code at run time based on available hardware features,
            # which can improve performance.
            # However, in spoa v4.0.8 it introduces two additional dependencies:
            # - USCiLab/cereal
            # - google/cpu_features
            # To avoid these additional dependencies, we turn off the dispatcher feature for now.
            # We could turn it back on if we see significant performance degradation in this area.
            spoaBuildFlags = "-DCMAKE_BUILD_TYPE=Release -Dspoa_optimize_for_portability=ON -Dspoa_build_tests=OFF"

        # Build the shared library.
        print("\n******** Building spoa shared library")
        os.mkdir("build")
        os.chdir("build")
        runCommand("cmake ../spoa-4.0.8 -DBUILD_SHARED_LIBS=ON " + spoaBuildFlags)
        runCommand("make -j all")
        runCommand("make install")
        
        # Build the static library.
        os.chdir("..")
        print("\n******** Building spoa static library")
        os.mkdir("build-static")
        os.chdir("build-static")
        runCommand("cmake ../spoa-4.0.8 -DBUILD_SHARED_LIBS=OFF " + spoaBuildFlags)
        runCommand("make -j all")
        runCommand("make install")
        
        # Change back to the original directory.
        os.chdir(oldDirectory)


    
installAptPackages() 
installPybind11() 
installSpoa()
  
# Make sure the newly created libraries are immediately visible to the loader.
runCommand("ldconfig")

print("Installation of Shasta prerequisies completed successfully.")
