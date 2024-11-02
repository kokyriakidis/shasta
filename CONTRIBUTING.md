# Contributing to shasta

Contributions of code, ideas, computational experiments,
or documentation are welcome. 
To contribute, please use the standard GitHub Pull Request process
and take into consideration the coding suggestions below.

#### Contributing C++ code

* Add at least a minimal amount
of comments to explain or clarify what the code does. 
* Use blank lines generously to separate code blocks.
* Set your editor to convert tabs to 4 spaces.
* To minimize name pollution of the global namespace,
all code should be in the `shasta` namespace.
Do not use `using` directives at global scope.
* If you need to add a large amount of code to an existing
file, consider creating a new file instead.
* Keep in mind that the code is currently compiled
using the C++20 standard. Use features in that standard
liberally. However, you cannot use C++ features
that are being considered for the C++23 standard.
* Build and test your code using the current long term support Ubuntu release.
* For functionality not provided
by the standard libraries, you can also use the Boost libraries.
* If possible, avoid adding dependencies on other packages.

#### Contributing ideas, computational experiments, comments, or criticism
* Please use the GitHub
[Issues](https://github.com/paoloshasta/shasta/issues)
section of the GitHub repository as appropriate.

#### Contributing documentation

* Put your contribution in the `docs`
directory, so it becomes automatically visible in GitHub Pages.
* Feel free to create new `.html` files, if appropriate.
* Attempt to harmonize your new materials with the existing
documentation, if possible.

