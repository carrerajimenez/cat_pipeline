
# Contributing to the Ca II Triplet EW Pipeline

First off, thanks for taking the time to contribute\!

The goal of this project is to provide a robust, flexible, and easy-to-use pipeline for analyzing Ca II Triplet spectra. Whether you are fixing a bug, adding a new line profile model, or improving the documentation, your help is appreciated.

## Table of Contents

1.  [Reporting Bugs]([https://github.com/carrerajimenez/cat_pipeline/blob/main/CONTRIBUTING.md/search?q=%23reporting-bugs)
2.  [Suggesting Enhancements]([https://github.com/carrerajimenez/cat_pipeline/blob/main/CONTRIBUTING.md/search?q=%23suggesting-enhancements)
3.  [Contributing Code](https://github.com/carrerajimenez/cat_pipeline/blob/main/CONTRIBUTING.md/search?q=%23contributing-code)
      * [General Workflow](https://github.com/carrerajimenez/cat_pipeline/blob/main/CONTRIBUTING.md/search?q=%23general-workflow)
      * [Adding a New Line Profile](https://github.com/carrerajimenez/cat_pipeline/blob/main/CONTRIBUTING.md/search?q=%23adding-a-new-line-profile)
      * [Adding a Pre-Normalization Method](https://github.com/carrerajimenez/cat_pipeline/blob/main/CONTRIBUTING.md/search?q=%23adding-a-pre-normalization-method)
4.  [Documentation]([https://github.com/carrerajimenez/cat_pipeline/blob/main/CONTRIBUTING.md/search?q=%23documentation)

-----

## Reporting Bugs

If you are using the code and reading the documentation you are already contributing to its development. If you find a bug in the code or an error in the results, documentation, please open an issue here on GitHub and we will try to help you out. To help us solve it quickly, please include:

  * **A clear title and description.**
  * **Steps to reproduce:** Provide the specific configuration (values in your `config.yaml`) that caused the error.
  * **Error Logs:** Paste the full error trace from the terminal.
  * **System Info:** Your operating system and Python version.

> **Tip:** If possible, try to reproduce the error using a small subset of your data or standard example files to isolate the issue.

## Suggesting Enhancements

If you have an idea for a new feature (e.g., a new minimization algorithm, a different way to handle SNR, or a new plot style):

  * Open an issue to discuss it first.
  * Explain *why* this enhancement would be useful.
  * If you have a specific implementation in mind, feel free to describe it.

## Contributing Code

We welcome pull requests (PRs) for bug fixes and new features\!

### General Workflow

1.  **Fork** the repository to your own GitHub account.
2.  **Clone** the project to your local machine.
3.  Create a new **branch** for your feature or fix.
4.  Make your changes.
5.  **Test** your changes by running the pipeline on a dataset to ensure it completes successfully and produces valid output.
    ```bash
    python cat_pipeline.py -c config.yaml
    ```
6.  **Commit** your changes with descriptive messages.
7.  Push to your fork and submit a **Pull Request**.

### Adding a New Line Profile

One of the most powerful features of this pipeline is the ability to add custom line profiles.

1.  Open `line_fitters.py`.
2.  Define your analytical function (e.g., `my_custom_profile`) and its residual function.
3.  Create a parameter creator function (e.g., `create_params_my_custom`) to set initial guesses and bounds.
4.  Create an EW calculation function (e.g., `calculate_ew_my_custom`).
5.  **Register it:** Add your new functions to the `PROFILE_REGISTRY` dictionary at the bottom of `line_fitters.py`.
6.  Update `README.md` to list the new option.

### Adding a Pre-Normalization Method

If you want to implement a new continuum fitting algorithm (e.g., a spline fit or a different sigma-clipping approach):

1.  Open `pre_normalizer.py`.
2.  Add your new method as a private function inside the `SpectrumNormalizer` class (e.g., `_fit_spline`).
3.  Update the `normalize` dispatcher method to call your function when a specific `order` or config flag is set.
4.  Update `config.yaml` comments to explain how to trigger this new method.

## Documentation

Documentation is just as important as code.

  * If you change how the pipeline is run or configured, please update `README.md`.
  * If you add new parameters to `config.yaml`, make sure to add comments explaining them in the example file.
  * Ensure your code is commented, especially complex mathematical operations in `line_fitters.py`.
