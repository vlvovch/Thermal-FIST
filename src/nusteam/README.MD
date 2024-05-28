# Calculating thermodynamics in the van der Waals HRG model

This example calculates thermodynamics of the HRG model at zero chemical potential with van der Waals interactions as implemented in paper 
- V. Vovchenko, M.I. Gorenstein, H. Stoecker, *van der Waals Interactions in Hadron Resonance Gas: From Nuclear Matter to Lattice QCD*, [Phys.Rev.Lett. **118** (2017) 182301](https://doi.org/10.1103/PhysRevLett.118.182301), [arXiv:1609.03975 [hep-ph]](https://arxiv.org/abs/1609.03975)

The file [NuSTEAM-HRG-EoS.cpp](NuSTEAM-HRG-EoS.cpp) performs the calculation by linking to Thermal-FIST library. Calculation in the ideal HRG limit is readily implemented. It is users' task to implement ``SetVDWInteractionParameters`` function to enable vdW interaction calculation.
If you are stuck, see [NuSTEAM-HRG-EoS-completed.cpp](NuSTEAM-HRG-EoS-completed.cpp) for solution.

## Usage

1. Clone the ``nusteam`` branch from GitHub.
```bash
git clone -b nusteam https://github.com/vlvovch/Thermal-FIST.git
cd Thermal-FIST
```

2. Use CMake to build library and the ``NuSTEAM-HRG-EoS`` executable. Adjust the commands for your operating system as needed.
```bash
mkdir build
cd build
cmake ../
make NuSTEAM-HRG-EoS -j
```

3. Run the ``NuSTEAM-HRG-EoS`` executable in the ideal gas limit
```bash
./bin/nusteam/NuSTEAM-HRG-EoS 0. 0.
```
Compate the output with the dashes lines in Figures 1a and 2b of the paper [https://arxiv.org/pdf/1609.03975](https://arxiv.org/pdf/1609.03975).

4. Try to run the script with non-zero interaction parameters. The program will throw an error. You will need to implement the function below to set the interaction for different hadron pairs
```cpp
// Function which sets the a and b parameters of the van der Waals HRG model 
// for baryon-baryon and antibaryon-antibaryon pairs
void SetVDWInteractionParameters(ThermalModelBase *model, double a, double b)
{
  // Iterate over all hadron pairs
  for (int i1 = 0; i1 < model->TPS()->Particles().size(); ++i1) {
    for (int i2 = 0; i2 < model->TPS()->Particles().size(); ++i2) {
      // Baryon charge of first species
      int B1 = model->TPS()->Particles()[i1].BaryonCharge();
      // Baryon charge of second species
      int B2 = model->TPS()->Particles()[i2].BaryonCharge();

      // TODO part:
      // Comment the lines below and write code to
      // set the vdW parameters for baryon-baryon and antibaryon-antibaryon pairs
      // based on the values (signs) of B1 and B2
      // Use model->SetAttraction(i1, i2, a) for attraction and
      // model->SetRepulsion(i1, i2, b) for repulsion
      cout << "Implement setting the vdW parameters for baryon-baryon and antibaryon-antibaryon interaction" << endl;
      exit(1);
    }
  }
}
```
If you are stuck, see [NuSTEAM-HRG-EoS-completed.cpp](NuSTEAM-HRG-EoS-completed.cpp) for solution.

5. Perform the calculation within van der Waals HRG with parameters from [https://arxiv.org/pdf/1609.03975](https://arxiv.org/pdf/1609.03975).
```bash
./bin/nusteam/NuSTEAM-HRG-EoS 0.329 3.42
```
Compare the results with solid blue lines in Figures 1a and 2b of the paper [https://arxiv.org/pdf/1609.03975](https://arxiv.org/pdf/1609.03975).

*Copyright (C) 2024  Volodymyr Vovchenko*