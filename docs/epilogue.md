
# Support

In case of bugs and feature improvements, you are welcome to create a
new issue on the [github repo](https://github.com/pyscal/pyscal3). You
are also welcome to fix a bug or implement a feature. 

Any other questions or suggestions are welcome, please contact
[us](mailto:sarath.menon@pyscal.org).

`pyscal3` welcomes and appreciates contribution and extension to the
module. Rather than local modifications, we request that the
modifications be submitted through a pull request, so that the module
can be continuously improved.

## Reporting and fixing bugs

In case a bug is found in the module, it can be reported on the [issues
page of the repository](https://github.com/pyscal/pyscal3/issues). Once a bug is reported, the status can once again monitored on
the issues page. Additionally, you are of course very welcome to fix any
existing bugs.

## New features

If you have an idea for new feature, you can submit a feature idea
through the [issues page of the
repository](https://github.com/pyscal/pyscal3/issues). As much as
information as you can provide about the new feauture would be greatly
helpful. Additionally, you could also work on feature requests already
on the issues page. The following instructions will help you get started
with local feature development.

### Setting up local environment

1.  The first step is to fork `pyscal3`. A detailed tutorial on forking can
    be found [here](https://help.github.com/en/articles/fork-a-repo).
    After forking, clone the repository to your local machine.
2.  We recommend creating a virtual environment to test new features or
    improvements to features. See this
    [link](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html)
    for help on managing environments.
3.  Once the environment is set up, you can create a new branch for your
    feature by `git checkout -b new_feauture`.
4.  Now implement the necessary feature.
5.  Once done, you can reinstall `pyscal3` by `pip install .`.
    After that please make sure that the existing tests work by running
    `pytest tests/` from the main module folder.
6.  If the tests work, you are almost done! If the new feature is not
    covered in existing tests, you can to write a new test in the tests
    folder. `pyscal3` uses pytest for tests. [This
    link](http://doc.pytest.org/en/latest/getting-started.html) will
    help you get started.
7.  Add the necessary docstrings for the new functions implemented.
    `pyscal3` uses the [numpy docstring
    format](https://numpydoc.readthedocs.io/en/latest/format.html) for
    documentation.
8.  Bonus task: Set up few examples that document how the feature works
    in the `docs` folder and link it to the examples section.
9.  Final step - Submit a pull request through github. Before you
    submit, please make sure that the new feature is documented and has
    tests. Once the request is submitted, automated tests would be done.
    If all tests are successful, your feauture will be incorporated to calphy and your contributions
    will be credited.

If you have trouble with any of the steps, or you need help, please
[send an email](mailto:sarath.menon@pyscal.org) and we will be happy to
help! 

# Acknowledgements

## Developer

-   [Sarath Menon](http://sarathmenon.me)  
    sarath.menon@pyscal.org


## Contributers

Please see the complete list of contributers [here](https://github.com/pyscal/pyscal3/graphs/contributors).


## Acknowledgements

- [Bond order analysis](https://github.com/WolfgangLechner/StructureAnalysis) code for the inspiration and the base for what later grew to be `pyscal`. 
- [Voro++](math.lbl.gov/voro++/) and [pybind11](https://pybind11.readthedocs.io/en/stable/) for developing the great tools that we could use in `pyscal3`.  
- [E-CAM High Throughput Computing ESDW](https://www.e-cam2020.eu/event/4424/?instance_id=71) held in [Turin](https://www.polito.it/?lang=en) in 2018 and 2019 for programming help, especially David W.H. Swenson and Alan O'Cais. 
- Scholarship from the International Max Planck Research School for Interface Controlled Materials for Energy Conversion for funding the initial stages of this work.
- [Interdisciplinary Centre for Advanced Materials Simulation](http://www.icams.de/content), at the [Ruhr University Bochum](https://www.ruhr-uni-bochum.de/en), Germany for the resources.
- Funding by the NFDI consortium [NFDI-MatWerk](https://nfdi-matwerk.de/) in the context of the work of the association German National Research Data Infrastructure (NFDI) e.V. NFDI is financed by the Federal Republic of Germany and the 16 federal states and funded by the Federal Ministry of Education and Research (BMBF) - funding code M532701 / the Deutsche Forschungsgemeinschaft (DFG, German Research Foundation) - [project number 460247524](https://gepris.dfg.de/gepris/projekt/460247524?language=en).
- Abril Azócar Guzmán for the pyscal logo.

# LICENSE

BSD 3-Clause License

Copyright (c) 2024, Sarath Menon
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.