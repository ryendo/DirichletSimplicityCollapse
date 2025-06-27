# Rigorous Verification of the Simplicity of the Second Dirichlet Eigenvalue

## Overview

This project provides the MATLAB implementation for the computer-assisted proof presented in the paper: **"The second Dirichlet eigenvalue is simple on every non-equilateral triangle"** by Ryoki Endo and Xuefeng Liu.


## Background

The relationship between the shape of a domain and its Laplacian eigenvalues is a central topic in spectral geometry. For a triangular domain, it is known that the second and third Dirichlet eigenvalues are equal (`λ₂ = λ₃`) for an equilateral triangle due to its symmetry. The conjecture, posed by R. Laugesen and B. Siudeja, asserts that this degeneracy is immediately broken for any non-equilateral triangle.

Previous work by the authors had already confirmed the conjecture for triangles whose minimum normalized height is greater than or equal to `tan(π/60)/2`. This project provides the final piece of the puzzle by proving the conjecture for triangles with a minimum normalized height **less than or equal to `tan(π/60)/2`**.

## Methodology

The proof is computationally intensive and relies on a clever transformation of the original problem.

#### 1\. Problem Transformation

Instead of directly comparing `λ₂(s, t)` and `λ₃(s, t)`, the problem is reformulated by analyzing a related quantity, `μₖᵗ`, defined as:

$$\mu_k^t = t^{4/3} \left( \lambda_k^t - \frac{\pi^2}{t^2} \right)$$

where `t` is the height of the triangle. Proving `λ₂ < λ₃` is equivalent to proving `μ₂ᵗ < μ₃ᵗ`. This transformation is crucial because `μₖᵗ` converges to a finite value as the triangle collapses (`t → 0`), making it suitable for numerical analysis.

#### 2\. Bounding `μₖᵗ` with an Inequality

The core of the proof lies in establishing rigorous upper and lower bounds for `μₖᵗ`. This is achieved by relating the original problem to two simpler, one-dimensional Schrödinger operator eigenvalue problems (Problem 2 and Problem 3). This relationship yields the central inequality for `t ∈ (0, t₀]`, where `t₀ = tan(π/60)/2`:

$$\frac{\overline{\mu}_k(s)}{1 + t_0^{2/3}/(3\pi^2) \cdot \overline{\mu}_k(s)} \le \mu_k^t(s) \le \hat{\mu}_k^{t_0}(s)$$

Here, `μ̂ₖᵗ⁰` and `μ̄ₖ` are eigenvalues of the related 1D problems that can be bounded rigorously.

#### 3\. Computer-Assisted Proof

The proof proceeds in two main computational steps:

1.  **Compute a rigorous upper bound for `μ₂ᵗ(s)`** by calculating an upper bound for `μ̂₂ᵗ⁰(s)`.
2.  **Compute a rigorous lower bound for `μ₃ᵗ(s)`** by first finding a lower bound for `μ̄₃(s)`.

By showing that the upper bound for `μ₂ᵗ` is strictly less than the lower bound for `μ₃ᵗ`, the conjecture is proven for the case of collapsing triangles.

## Project Structure

```
.
├── README.md
├── functions/
│   ├── func_left_hand_side.m
│   └── ...
├── mode_swith_interface/
│   ├── I_intval.m
│   └── ...
├── my_intlab_mode_config.m
├── results/
│   └── (Output files will be generated here)
├── step1_compute_mu_hat.m
├── step2_1_verification_of_mu3_s_at_0.m
├── step2_2_investigation_of_the_branch_mu3_s_for_s_in_0_1.m
├── upper_bound_matrix/
│   └── (Symbolic computation cache)
└── verified_eig_estimation/
    └── ...
```

  * `my_intlab_mode_config.m`: Initializes INTLAB and sets up the project paths.
  * `step1_compute_mu_hat.m`: **(Step 1)** Implements the proof described in Section 4.1 of the paper to compute the upper bound of `μ₂ᵗ`.
  * `step2_1_verification_of_mu3_s_at_0.m`: **(Step 2a)** Verifies the lower bound of `μ̄₃(s)` at `s=0`, as described in Section 4.2.
  * `step2_2_investigation_of_the_branch_mu3_s_for_s_in_0_1.m`: **(Step 2b)** Confirms that `μ̄₃(s)` remains above the critical value for all `s ∈ [0, 1)`, completing the proof for the lower bound.
  * `functions/`: Contains core mathematical helper functions, such as the implementation of the implicit equation from the paper.
  * `mode_swith_interface/`: A set of wrapper functions for INTLAB to facilitate interval arithmetic.
  * `upper_bound_matrix/`: Caches symbolic computations to speed up subsequent runs of `step1_compute_mu_hat.m`.
  * `results/`: Default directory for storing output logs and data files.
  * `verified_eig_estimation/`: Helper functions for rigorous eigenvalue estimations.

## Prerequisites

  * **MATLAB** (R2019b or later recommended)
  * **MATLAB Symbolic Math Toolbox**
  * **INTLAB - The INterval LABoratory** (Version 12 is used in this project). It must be downloaded separately.

## Installation and Setup

1.  **Clone the repository:**
    ```bash
    git clone <repository-url>
    ```
2.  **Install INTLAB:**
      * Download INTLAB from the [INTLAB website](https://www.tuhh.de/ti3/rump/intlab/).
      * Unpack the toolbox into the project's root directory. The configuration script expects it to be in a folder named `Intlab_V12`.
3.  **Configure the project:**
      * Open MATLAB and navigate to the project directory.
      * Run the configuration script:
        ```matlab
        my_intlab_mode_config
        ```
        This will add all necessary folders to the MATLAB path and initialize INTLAB.

## How to Run the Proof

Execute the following scripts in order. **Note:** These computations are intensive and may take a significant amount of time.

#### Step 1: Compute the Upper Bound for `μ₂ᵗ`

Run the script for Step 1. This computes the rigorous upper bound for `μ₂ᵗ` by evaluating `μ̂₂ᵗ⁰(s)`.

```matlab
step1_compute_mu_hat
```

  * **What it does:** This script first performs a one-time symbolic integration to define the matrices for the Rayleigh-Ritz method. It then iterates through 100 subintervals of `s ∈ [0, 1]` to compute a verified upper bound for the eigenvalue.
  * **Output:** The results are saved in `results/step1_1_results_N17_Ns100.csv`. Upon completion, the code will have rigorously shown that **`μ₂ᵗ(s) ≤ 21.091`** for all relevant `s` and `t`.

#### Step 2a: Verify the Lower Bound for `μ̄₃(s)` at s=0

Run the script to establish the starting point for the lower bound of `μ₃ᵗ`.

```matlab
step2_1_verification_of_mu3_s_at_0
```

  * **What it does:** This script uses the `verifynlssall` solver from INTLAB to find the rigorous roots of the implicit equation defining `κ₃(0)`.
  * **Output:** The console will display the verified enclosure for the third root (`κ₃`) and the resulting lower bound for `μ̄₃(0)`, confirming it is **greater than 23.5**.

#### Step 2b: Investigate the `μ̄₃(s)` Branch

Run the final script to ensure the lower bound holds for all `s ∈ [0, 1)`.

```matlab
step2_2_investigation_of_the_branch_mu3_s_for_s_in_0_1
```

  * **What it does:** This script performs a fine-grained sweep over the parameter `s`, rigorously proving that the function `f_s(κ)` (where `κ` corresponds to `μ̄₃ = 23.5`) never crosses zero. This confirms that `μ̄₃(s)` remains strictly above 23.5.
  * **Output:** A log file is created in the `results/` directory, documenting the successful verification for each subinterval of `s`.

The successful execution of these three steps provides a complete, rigorous, and verifiable proof that `μ₂ᵗ < μ₃ᵗ`, which in turn proves `λ₂ < λ₃` for all nearly degenerate, non-equilateral triangles.

## Citation

If you use this work, please cite the original paper:

R. Endo, X. Liu, "The second Dirichlet eigenvalue is simple on every non-equilateral triangle," arXiv preprint arXiv:2503.06786v2 (2025).

## Authors of the Paper and Code

  * Ryoki Endo (Niigata University)
  * Xuefeng Liu (Tokyo Woman's Christian University)