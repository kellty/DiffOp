# Implementation of Differential Operator Approach to Function-on-function Regression

## Data
The ERA5 data of TX1day temperature from [https://doi.org/10.3929/ethz-b-000571107](https://doi.org/10.3929/ethz-b-000571107) was used for this work. The file "***TX1.RData***" was obtained after pre-processing by "***real_prep.R***".

## Code
- "***func.R***": Build functions to implement the regression analysis for differential operators.
- "***simu.R***" and "***simu_np.R***": Perform simulation studies under different settings.
- "***real.R***": Perform real data analysis.

## Workflow
First of all, make sure that "***func.R***" is in the working directory.
- Simulation: Run "***simu.R***" and "***simu_np.R***".
- Real data example: Run "***real.R***".
